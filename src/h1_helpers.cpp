#include <cassert>

#include <ngstd.hpp>

#include "h1_helpers.hpp"

#include "h1.hpp"
#include "build_h1.hpp"

namespace h1amg
{

void IterativeCollapse(
    const ngstd::Array<INT<2>>& edge_to_vertices, const ngstd::Array<double>& edge_collapse_weight,
    const ngstd::Array<double>& vertex_collapse_weight, shared_ptr<ngstd::BitArray> free_dofs,
    ngstd::Array<bool>& edge_collapse, ngstd::Array<bool>& vertex_collapse)
{
  static Timer Titer_coll("H1-AMG::IterativeCollapse");
  RegionTimer Riter_coll(Titer_coll);

  assert(edge_to_vertices.Size() == edge_collapse_weight.Size());
  int nr_vertices = vertex_collapse_weight.Size();
  int nr_edges = edge_collapse_weight.Size();

  ngstd::Array<int> vertex_to_edge(nr_vertices);
  vertex_to_edge = -1;

  edge_collapse.SetSize(nr_edges);
  edge_collapse = false;
  vertex_collapse.SetSize(nr_vertices);
  vertex_collapse = false;

  bool changed;
  int unused = -1;
  int ground = -2;

  do {
    changed = false;

    for (int edge = 0; edge < nr_edges; ++edge) {
      if (edge_collapse_weight[edge] > 0.1
          &&
          (!free_dofs || (   (*free_dofs)[edge_to_vertices[edge][0]]
                          && (*free_dofs)[edge_to_vertices[edge][1]])))
      {
        if (vertex_to_edge[edge_to_vertices[edge][0]] == unused &&
            vertex_to_edge[edge_to_vertices[edge][1]] == unused)
        {
          edge_collapse[edge] = true;
          vertex_to_edge[edge_to_vertices[edge][0]] = edge;
          vertex_to_edge[edge_to_vertices[edge][1]] = edge;
          changed = true;
        }

        else {
          for (int j = 0; j < 2; ++j) {
            int vertex1 = edge_to_vertices[edge][j];
            int vertex2 = edge_to_vertices[edge][1-j];

            if (vertex_to_edge[vertex1] != unused &&
                vertex_to_edge[vertex1] != ground &&
                vertex_to_edge[vertex2] == unused &&
                edge_collapse_weight[edge] >
                edge_collapse_weight[vertex_to_edge[vertex1]])
            {
              int old_edge = vertex_to_edge[vertex1];
              edge_collapse[edge] = true;
              edge_collapse[old_edge] = false;

              vertex_to_edge[edge_to_vertices[old_edge][0]] = unused;
              vertex_to_edge[edge_to_vertices[old_edge][1]] = unused;
              vertex_to_edge[vertex1] = edge;
              vertex_to_edge[vertex2] = edge;

              changed = true;
              break;
            }
            else if (vertex_to_edge[vertex1] == ground &&
                      vertex_to_edge[vertex2] == unused &&
                      edge_collapse_weight[edge] >
                      vertex_collapse_weight[vertex1])
            {
              vertex_collapse[vertex1] = false;
              vertex_to_edge[vertex1] = edge;
              vertex_to_edge[vertex2] = edge;
              edge_collapse[edge] = true;
              changed = true;
              break;
            }
          }
        }

      }
    }

    for (int vertex = 0; vertex < nr_vertices; ++vertex) {
      if (vertex_to_edge[vertex] == unused
          && (vertex_collapse_weight[vertex] > 0.1 ||
              (free_dofs && !(*free_dofs)[vertex])))
      {
        vertex_to_edge[vertex] = ground;
        vertex_collapse[vertex] = true;
        changed = true;
      }
      else if (vertex_to_edge[vertex] >= 0 &&
                vertex_collapse_weight[vertex] >
                edge_collapse_weight[vertex_to_edge[vertex]])
      {
        edge_collapse[vertex_to_edge[vertex]] = false;
        vertex_to_edge[vertex] = ground;
        vertex_collapse[vertex] = true;
        changed = true;
      }
    }

  }
  while (changed);

}

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

  for (int i = 0; i < nr_edges; ++i) {
    for (int j = 0; j < 2; ++j) {
      vertex_strength[edge_to_vertices[i][j]] += weights_edges[i];
    }
  }

  for (int i = 0; i < nr_vertices; ++i) {
    double current_weight = weights_vertices[i];
    vertex_strength[i] += current_weight;

    // when vertex weight is not 0.0, then also vertex_strength of that vertex
    // can't be 0.0
    if (current_weight != 0.0)
      vertex_collapse_weight[i] = current_weight / vertex_strength[i];
  }


  for (int i = 0; i < nr_edges; ++i) {
    double vstr1 = vertex_strength[edge_to_vertices[i][0]];
    double vstr2 = vertex_strength[edge_to_vertices[i][1]];

    // share of the edge weight to the vertex strength
    // same as: weights_edges[i] / vstr1 + weights_edges[i] / vstr2
    edge_collapse_weight[i] = weights_edges[i] * (vstr1+vstr2) / (vstr1 * vstr2);
  }
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

  for (int vertex = 0; vertex < nv; ++vertex)
    connected[vertex] = vertex;

  for (int edge = 0; edge < nr_edges; ++edge)
  {
    if (edge_collapse[edge])
    {
      int vertex1 = edge_to_vertices[edge][0];
      int vertex2 = edge_to_vertices[edge][1];
      if (vertex2>vertex1) {
        connected[vertex2] = vertex1;
      }
      else {
        connected[vertex1] = vertex2;
      }
    }
  }

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

  *testout << "vertex_coarse before | after fillup:" << endl;
  for (int vertex = 0; vertex < nv; ++vertex)
  {
    *testout << vertex << ":  " << vertex_coarse[vertex] << " | ";
    if (connected[vertex] != vertex) {
      vertex_coarse[vertex] = vertex_coarse[connected[vertex]];
    }
    *testout << vertex_coarse[vertex] << endl;
  }

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
  ngstd::ClosedHashTable<INT<2>, int> edge_coarse_table(nr_edges);
  edge_coarse.SetSize(nr_edges);

  // compute fine edge to coarse edge map (edge_coarse)
  for (int edge = 0; edge < nr_edges; ++edge)
  {
    int vertex1 = vertex_coarse[edge_to_vertices[edge][0]];
    int vertex2 = vertex_coarse[edge_to_vertices[edge][1]];

    // only edges where both coarse vertices are different and don't
    // collapse to ground will be coarse edges
    if (vertex1 != -1 && vertex2 != -1 && vertex1 != vertex2) {
      edge_coarse_table.Set (INT<2>(vertex1, vertex2).Sort(), -1);
    }
  }

  t1.Stop();

  t2a.Start();
  int cnt = 0;
  for (int i = 0; i < edge_coarse_table.Size(); ++i)
    if (edge_coarse_table.UsedPos(i))
      edge_coarse_table.SetData (i, cnt++);
  t2a.Stop();

  coarse_edge_to_vertices.SetSize(cnt);

  // compute coarse edge to vertex mapping to user for next recursive
  // coarsening
  t2.Start();
  ParallelFor(edge_coarse_table.Size(), [&] (int i)
  {
    if (edge_coarse_table.UsedPos(i)) {
      auto both =  edge_coarse_table.GetBoth(i);
      coarse_edge_to_vertices[both.second] = both.first;
    }
  });
  t2.Stop();

  t3.Start();
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
  t3.Stop();
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

  for (int i = 0; i < edge_to_vertices.Size(); i++) {
    if (edge_coarse[i] != -1) {
      weights_edges_coarse[edge_coarse[i]] += weights_edges[i];
    }
  }
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

  for (int fine_vertex = 0; fine_vertex < nr_vertices; ++fine_vertex) {
    int coarse_vertex = vertex_coarse[fine_vertex];
    if (coarse_vertex != -1) {
      weights_vertices_coarse[coarse_vertex] += weights_vertices[fine_vertex];
    }
  }

  for (int fine_edge = 0; fine_edge < nr_edges; ++fine_edge) {
    for (int i = 0; i < 2; ++i) {
      int cvertex1 = vertex_coarse[ edge_to_vertices[fine_edge][i] ];
      int cvertex2 = vertex_coarse[ edge_to_vertices[fine_edge][1-i] ];
      if (cvertex1 == -1 && cvertex2 != -1) {
        *testout << "edge " << fine_edge << " between cvert " << cvertex1
                 << " and " << cvertex2 << endl;
        weights_vertices_coarse[cvertex2] += weights_edges[fine_edge];
      }
    }
  }
}

}  // h1amg
