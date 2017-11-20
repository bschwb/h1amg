#include <cassert>

#include <ngstd.hpp>
using namespace ngstd;

#include "h1_helpers.hpp"

#include "h1.hpp"
#include "build_h1.hpp"

namespace h1amg
{

void ComputeCollapseWeights(
    const ngstd::Array<INT<2>>& edge_to_vertices, const ngstd::Array<double>& weights_edges,
    const ngstd::Array<double>& weights_vertices, ngstd::Array<double>& vertex_strength,
    ngstd::Array<double>& edge_collapse_weight, ngstd::Array<double>& vertex_collapse_weight)
{
  static Timer Tcoll_weights("H1-AMG::ComputeCollapseWeights");
  RegionTimer Rcoll_weights(Tcoll_weights);

  assert(edge_to_vertices.Size() == weights_edges.Size());
  int nr_edges = edge_to_vertices.Size();
  int nr_vertices = weights_vertices.Size();

  vertex_strength.SetSize(nr_vertices);
  vertex_strength = 0.0;
  edge_collapse_weight.SetSize(nr_edges);
  edge_collapse_weight = 0.0;
  vertex_collapse_weight.SetSize(nr_vertices);
  vertex_collapse_weight = 0.0;

  // TODO: Try switching loops to avoid branch misprediction?
  static Timer Tcweights_vertstr("H1-AMG::ComputeCollapseWeights::VertStrength");
  Tcweights_vertstr.Start();
  ParallelFor(nr_edges, [&] (int i) {
    for (int j = 0; j < 2; ++j) {
      AsAtomic(vertex_strength[edge_to_vertices[i][j]]) += weights_edges[i];
    }
  });
  Tcweights_vertstr.Stop();

  static Timer Tcweights_vcollweight("H1-AMG::ComputeCollapseWeights::VertCollWeight");
  Tcweights_vcollweight.Start();
  ParallelFor(nr_vertices, [&] (int i) {
    double current_weight = weights_vertices[i];
    vertex_strength[i] += current_weight;

    // when vertex weight is not 0.0, then also vertex_strength of that vertex
    // can't be 0.0
    if (current_weight != 0.0) {
      vertex_collapse_weight[i] = current_weight / vertex_strength[i];
    }
  });
  Tcweights_vcollweight.Stop();

  static Timer Tcweights_ecollweight("H1-AMG::ComputeCollapseWeights::EdgeCollWeight");
  Tcweights_ecollweight.Start();
  ParallelFor(nr_edges, [&] (int i) {
    double vstr1 = vertex_strength[edge_to_vertices[i][0]];
    double vstr2 = vertex_strength[edge_to_vertices[i][1]];

    // share of the edge weight to the vertex strength
    // same as: weights_edges[i] / vstr1 + weights_edges[i] / vstr2
    edge_collapse_weight[i] = weights_edges[i] * (vstr1+vstr2) / (vstr1 * vstr2);
  });
  Tcweights_ecollweight.Stop();
}


int ComputeFineToCoarseVertex(
  const ngstd::Array<INT<2>>& edge_to_vertices, int nv,
  const ngstd::Array<bool>& edge_collapse, const ngstd::Array<bool>& vertex_collapse,
  ngstd::Array<int>& vertex_coarse )
{
  static Timer Tf2c_verts("H1-AMG::ComputeFineToCoarseVertex");
  RegionTimer Rf2c_verts(Tf2c_verts);

  int nr_edges = edge_to_vertices.Size();
  int nr_coarse_vertices = 0;

  ngstd::Array<int> connected(nv);
  vertex_coarse.SetSize(nv);

  vertex_coarse = -4;

  static Timer Tf2cv_sc("H1-AMG::ComputeFineToCoarseVertex::SelfConnect");
  Tf2cv_sc.Start();
  ParallelFor(nv, [&connected] (int vertex) {
    connected[vertex] = vertex;
  });
  Tf2cv_sc.Stop();

  static Timer Tf2cv_cc("H1-AMG::ComputeFineToCoarseVertex::CoarseConnection");
  Tf2cv_cc.Start();
  // TODO: not sure if we can parallize this
  // Is it possible for more than 1 edge of a vertex to collapse?
  ParallelFor(nr_edges, [&](int edge) {
    if (edge_collapse[edge])
    {
      int vertex1 = edge_to_vertices[edge][0];
      int vertex2 = edge_to_vertices[edge][1];
      if (vertex2>vertex1) {
        AsAtomic(connected[vertex2]) = vertex1;
      }
      else {
        AsAtomic(connected[vertex1]) = vertex2;
      }
    }
  });
  Tf2cv_cc.Stop();

  static Timer Tf2cv_cntcoarse("H1-AMG::ComputeFineToCoarseVertex::CountCoarse");
  Tf2cv_cntcoarse.Start();
  for (int vertex = 0; vertex < nv; ++vertex)
  {
    if (connected[vertex] == vertex)
    {
      if (vertex_collapse[vertex]) {
        vertex_coarse[vertex] = -1;
      }
      else {
        vertex_coarse[vertex] = nr_coarse_vertices++;
      }
    }
  }
  Tf2cv_cntcoarse.Stop();

  // *testout << "vertex_coarse before | after fillup:" << endl;
  static Timer Tf2cv_mapping("H1-AMG::ComputeFineToCoarseVertex::Mapping");
  Tf2cv_mapping.Start();
  ParallelFor(nv, [&connected, &vertex_coarse] (int vertex) {
    // *testout << vertex << ":  " << vertex_coarse[vertex] << " | ";
    if (connected[vertex] != vertex) {
      vertex_coarse[vertex] = vertex_coarse[connected[vertex]];
    }
    // *testout << vertex_coarse[vertex] << endl;
  });
  Tf2cv_mapping.Stop();

  return nr_coarse_vertices;
}


void ComputeFineToCoarseEdge(
  const ngstd::Array<INT<2>>& edge_to_vertices, const ngstd::Array<int>& vertex_coarse,
  ngstd::Array<int>& edge_coarse, ngstd::Array<INT<2>>& coarse_edge_to_vertices)
{
  static Timer t("ComputeFine2CoarseEdge"); RegionTimer reg(t);

  static Timer t1("ComputeFineToCoarseEdge 1");
  static Timer t2("ComputeFineToCoarseEdge 2");
  static Timer t2a("ComputeFineToCoarseEdge 2a");
  static Timer t3("ComputeFineToCoarseEdge 3");

  t1.Start();
  int nr_edges = edge_to_vertices.Size();
  edge_coarse.SetSize(nr_edges);
  /*
  ngstd::ClosedHashTable<INT<2>, int> edge_coarse_table(nr_edges);

  // compute fine edge to coarse edge map (edge_coarse)
  for (int edge = 0; edge < nr_edges; ++edge)
  {
    auto verts = edge_to_vertices[edge];
    int vertex1 = vertex_coarse[verts[0]];
    int vertex2 = vertex_coarse[verts[1]];

    // only edges where both coarse vertices are different and don't
    // collapse to ground will be coarse edges
    if (vertex1 != -1 && vertex2 != -1 && vertex1 != vertex2) {
      edge_coarse_table.Set(INT<2>(vertex1, vertex2).Sort(), -1);
    }
  }
  */

  ngstd::ParallelHashTable<INT<2>, int> edge_coarse_table;
  // compute fine edge to coarse edge map (edge_coarse)
  // for (int edge = 0; edge < nr_edges; ++edge)
  ParallelFor (nr_edges, [&] (int edge)
  {
    auto verts = edge_to_vertices[edge];
    int vertex1 = vertex_coarse[verts[0]];
    int vertex2 = vertex_coarse[verts[1]];

    // only edges where both coarse vertices are different and don't
    // collapse to ground will be coarse edges
    if (vertex1 != -1 && vertex2 != -1 && vertex1 != vertex2) {
      edge_coarse_table.Do(INT<2>(vertex1, vertex2).Sort(), [](auto & val) { val = -1; });
    }
  });

  t1.Stop();
  t2a.Start();

  /*
  int cnt = 0;
  for (int i = 0; i < edge_coarse_table.Size(); ++i)
    if (edge_coarse_table.UsedPos(i))
      edge_coarse_table.SetData (i, cnt++);
  */

  /*
  int cnt = 0;
  edge_coarse_table.Iterate
    ( [&cnt] (INT<2> key, int & val)
      {
        val = cnt++;
      });
  */
  Array<int> prefixsums(edge_coarse_table.NumBuckets());
  size_t sum = 0;
  for (size_t i = 0; i < edge_coarse_table.NumBuckets(); i++)
    {
      prefixsums[i] = sum;
      sum += edge_coarse_table.Used(i);
    }
  coarse_edge_to_vertices.SetSize(sum);
  ParallelFor (edge_coarse_table.NumBuckets(),
               [&] (size_t nr)
               {
                 int cnt = prefixsums[nr];
                 edge_coarse_table.Iterate
                   (nr, [&cnt] (INT<2> key, int & val)
                    {
                      val = cnt++;
                    });
               });

  t2a.Stop();


  // compute coarse edge to vertex mapping to user for next recursive
  // coarsening
  t2.Start();

  /*
  ParallelFor(edge_coarse_table.Size(), [&] (int i)
  {
    if (edge_coarse_table.UsedPos(i)) {
      auto both =  edge_coarse_table.GetBoth(i);
      coarse_edge_to_vertices[both.second] = both.first;
    }
  });
  */
  ParallelFor (edge_coarse_table.NumBuckets(),
               [&] (size_t nr)
               {               
                 edge_coarse_table.Iterate
                   (nr, [&coarse_edge_to_vertices] (INT<2> key, int val)
                    {
                      coarse_edge_to_vertices[val] = key;
                    });
               });
  
  t2.Stop();

  t3.Start();
  /*
  ParallelFor(nr_edges, [&] (int edge)
  {
    int vertex1 = vertex_coarse[ edge_to_vertices[edge][0] ];
    int vertex2 = vertex_coarse[ edge_to_vertices[edge][1] ];

    // only edges where both coarse vertices are different and don't
    // collapse to ground will be coarse edges
    if (vertex1 != -1 && vertex2 != -1 && vertex1 != vertex2) {
      edge_coarse[edge] = edge_coarse_table.Get (INT<2>(vertex1, vertex2).Sort());
    }
    else {
      edge_coarse[edge] = -1;
    }
  });
  */

  ParallelFor(nr_edges, [&] (int edge)
  {
    int vertex1 = vertex_coarse[ edge_to_vertices[edge][0] ];
    int vertex2 = vertex_coarse[ edge_to_vertices[edge][1] ];

    // only edges where both coarse vertices are different and don't
    // collapse to ground will be coarse edges
    if (vertex1 != -1 && vertex2 != -1 && vertex1 != vertex2) {
      // edge_coarse[edge] = edge_coarse_table.Get (INT<2>(vertex1, vertex2).Sort());
      /*
      int getval = 0;
      edge_coarse_table.Do(INT<2>(vertex1, vertex2).Sort(), [&getval] (int val) { getval = val; });
      edge_coarse[edge] = getval;
      */
      /*
      edge_coarse[edge] =
        edge_coarse_table.Do(INT<2>(vertex1, vertex2).Sort(), [] (int val) { return val; });
      */
      edge_coarse[edge] = edge_coarse_table.Get(INT<2>(vertex1, vertex2).Sort());
    }
    else {
      edge_coarse[edge] = -1;
    }
  });
  // cout << "edge_coarse = " << edge_coarse << endl;
  t3.Stop();

  // find costs:
  /*
  for (size_t i = 0; i < edge_coarse_table.NumBuckets(); i++)
    cout << i << ": " << edge_coarse_table.Used(i) << "/" <<  edge_coarse_table.BucketSize(i) << endl;
  
  size_t costs = 0;
  size_t used = 0;
  for (size_t edge = 0; edge < nr_edges; edge++)
  {
    int vertex1 = vertex_coarse[ edge_to_vertices[edge][0] ];
    int vertex2 = vertex_coarse[ edge_to_vertices[edge][1] ];
    if (vertex1 != -1 && vertex2 != -1 && vertex1 != vertex2) {
      {
        costs += edge_coarse_table.GetCosts(INT<2>(vertex1, vertex2).Sort());
        used++;
      }
    }
  }
  cout << "avg hashing costs = " << double(costs)/used << endl;
  */
}

void ComputeCoarseWeightsEdges(
  const ngstd::Array<INT<2>>& edge_to_vertices,
  const ngstd::Array<INT<2>>& coarse_edge_to_vertices, const ngstd::Array<int>& edge_coarse,
  const ngstd::Array<double>& weights_edges, ngstd::Array<double>& weights_edges_coarse)
{
  static Timer Tcoarse_eweights("H1-AMG::ComputeCoarseWeightsEdges");
  RegionTimer Rcoarse_eweights(Tcoarse_eweights);

  weights_edges_coarse.SetSize(coarse_edge_to_vertices.Size());
  weights_edges_coarse = 0;

  ParallelFor(edge_to_vertices.Size(), [&] (int i) {
    if (edge_coarse[i] != -1) {
      AsAtomic(weights_edges_coarse[edge_coarse[i]]) += weights_edges[i];
    }
  });
}


void ComputeCoarseWeightsVertices(
  const ngstd::Array<INT<2>>& edge_to_vertices, const ngstd::Array<int>& vertex_coarse,
  const int nr_coarse_vertices, const ngstd::Array<double>& weights_edges,
  const ngstd::Array<double>& weights_vertices, ngstd::Array<double>& weights_vertices_coarse)
{
  static Timer Tcoarse_vweights("H1-AMG::ComputeCoarseWeightsVertices");
  RegionTimer Rcoarse_vweights(Tcoarse_vweights);

  int nr_vertices = weights_vertices.Size();
  int nr_edges = edge_to_vertices.Size();

  *testout << "nrv fine | coarse: " << nr_vertices << " " << nr_coarse_vertices << endl;
  *testout << "nr_edges: " << nr_edges << endl;
  weights_vertices_coarse.SetSize(nr_coarse_vertices);
  weights_vertices_coarse = 0;

  static Timer Tcvweights_fvweight("H1-AMG::ComputeCoarseWeightsVertices::AddFVWeight");
  Tcvweights_fvweight.Start();
  // for (int fine_vertex = 0; fine_vertex < nr_vertices; ++fine_vertex) {
  ParallelFor(nr_vertices, [&] (int fine_vertex) {
    int coarse_vertex = vertex_coarse[fine_vertex];
    if (coarse_vertex != -1) {
      AsAtomic(weights_vertices_coarse[coarse_vertex]) += weights_vertices[fine_vertex];
    }
  });
  Tcvweights_fvweight.Stop();

  static Timer Tcvweights_feweight("H1-AMG::ComputeCoarseWeightsVertices::AddFEWeight");
  Tcvweights_feweight.Start();
  ParallelFor(nr_edges, [&] (int fine_edge) {
    for (int i = 0; i < 2; ++i) {
      int cvertex1 = vertex_coarse[ edge_to_vertices[fine_edge][i] ];
      int cvertex2 = vertex_coarse[ edge_to_vertices[fine_edge][1-i] ];
      if (cvertex1 == -1 && cvertex2 != -1) {
        // *testout << "edge " << fine_edge << " between cvert " << cvertex1
        //          << " and " << cvertex2 << endl;
        AsAtomic(weights_vertices_coarse[cvertex2]) += weights_edges[fine_edge];
      }
    }
  });
  Tcvweights_feweight.Stop();
}

}  // h1amg
