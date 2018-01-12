#include <ngstd.hpp>
using namespace ngstd;
#include <la.hpp>
using namespace ngla;

#include "h1.hpp"
#include "h1_helpers.hpp"

#include "h1_smoothed_prol.hpp"
#include "dist1collapser.hpp"
#include "edge.hpp"
// #include "sample_sort.hpp"

#include "build_h1.hpp"


#include "concurrentqueue.h"
typedef moodycamel::ConcurrentQueue<size_t> TQueue; 
typedef moodycamel::ProducerToken TPToken; 
typedef moodycamel::ConsumerToken TCToken; 


namespace h1amg
{


  static TQueue queue;

  template <typename TFUNC>
  void RunParallelDependency (FlatTable<int> dag,
                              TFUNC func)
  {
    Array<atomic<int>> cnt_dep(dag.Size());

    for (auto & d : cnt_dep) 
      d.store (0, memory_order_relaxed);

    static Timer t_cntdep("count dep");
    t_cntdep.Start();
    ParallelFor (Range(dag),
                 [&] (int i)
                 {
                   for (int j : dag[i])
                     cnt_dep[j]++;
                 });
    t_cntdep.Stop();    

    atomic<size_t> num_ready(0), num_final(0);
    ParallelForRange (cnt_dep.Size(), [&] (IntRange r)
                      {
                        size_t my_ready = 0, my_final = 0;
                        for (size_t i : r)
                          {
                            if (cnt_dep[i] == 0) my_ready++;
                            if (dag[i].Size() == 0) my_final++;
                          }
                        num_ready += my_ready;
                        num_final += my_final;
                      });

    Array<int> ready(num_ready);
    ready.SetSize0();
    for (int j : Range(cnt_dep))
      if (cnt_dep[j] == 0) ready.Append(j);

    
    if (!task_manager)
      // if (true)
      {
        while (ready.Size())
          {
            int size = ready.Size();
            int nr = ready[size-1];
            ready.SetSize(size-1);
            
            func(nr);
            
            for (int j : dag[nr])
              {
                cnt_dep[j]--;
                if (cnt_dep[j] == 0)
                  ready.Append(j);
              }
          }
        return;
      }

    atomic<int> cnt_final(0);
    SharedLoop2 sl(Range(ready));

    task_manager -> CreateJob 
      ([&] (const TaskInfo & ti)
       {
        TPToken ptoken(queue); 
        TCToken ctoken(queue); 
        
        for (int i : sl)
          queue.enqueue (ptoken, ready[i]);

        while (1)
           {
             if (cnt_final >= num_final) break;

             int nr;
             if(!queue.try_dequeue_from_producer(ptoken, nr)) 
               if(!queue.try_dequeue(ctoken, nr))  
                 continue; 

	     // cnt_dep[nr]--;   // only for mem-sync
	     
             if (dag[nr].Size() == 0)
               cnt_final++;

             func(nr);

             for (int j : dag[nr])
               {
                 if (--cnt_dep[j] == 0)
                   queue.enqueue (ptoken, j);
               }
           }
       });
  }



  

  

using UPtrSMdbl = unique_ptr<SparseMatrixTM<double>>;
using SPtrSMdbl = shared_ptr<SparseMatrixTM<double>>;

UPtrSMdbl CreateProlongation(const Array<int>& vertex_coarse, int ncv, bool complx);

shared_ptr<H1AMG_Mat> BuildH1AMG(
    SPtrSMdbl sysmat, const Array<INT<2>>& edge_to_vertices,
    const Array<double>& weights_edges, const Array<double>& weights_vertices,
    shared_ptr<BitArray> free_dofs, const H1Options& h1_options)
{
  static Timer Tbuild_h1("H1-AMG::BuildH1AMG");
  RegionTimer Rbuild_h1(Tbuild_h1);

  cout << IM(5) << "H1 Sysmat nze: " << sysmat->NZE() << endl;
  cout << IM(5) << "H1 Sysmat nze per row: " << sysmat->NZE() / (double)sysmat->Height() << endl;
  auto ne = edge_to_vertices.Size();
  auto nv = weights_vertices.Size();

  Array<int> vertex_coarse;

  Array<double> vertex_strength;
  Array<double> edge_collapse_weight;
  Array<double> vertex_collapse_weight;

  Array<INT<2> > coarse_edge_to_vertices;
  Array<int> edge_coarse;

  Array<double> weights_edges_coarse;
  Array<double> weights_vertices_coarse;

  ComputeCollapseWeights(
      edge_to_vertices, weights_edges, weights_vertices, vertex_strength, edge_collapse_weight,
      vertex_collapse_weight);

  // IterativeCollapse(
  //     edge_to_vertices, edge_collapse_weight, vertex_collapse_weight, free_dofs, edge_collapse,
  //     vertex_collapse);

  /*
  static Timer Tdist1sorted("Dist1 Sorted Collapsing");
  Tdist1sorted.Start();
  static Timer t1("Dist1 Sorted Collapsing sorting");
  static Timer t2("Dist1 Sorted Collapsing work");

  Dist1Collapser collapser(nv, ne);

  t1.Start();
  Array<int> indices(ne);
  Array<Edge> edges(ne);
  ParallelFor (ne, [&] (size_t edge)
               {
                 indices[edge] = edge;
                 edges[edge] = Edge(edge, edge_to_vertices[edge][0], edge_to_vertices[edge][1]);
               });

  ngstd::SampleSortI(edge_collapse_weight, indices);
  t1.Stop();

  t2.Start();
  int vcnt = 0;
  for (int i = ne-1; i >= 0; --i) {
    auto edge = edges[indices[i]];

    if (vcnt >= nv/2.)
      break;
    if (edge_collapse_weight[edge.id] >= 0.01 && !collapser.AnyVertexCollapsed(edge)) {
      ++vcnt;
      collapser.CollapseEdge(edge);
    }
  }
  t2.Stop();
  Tdist1sorted.Stop();
  */


  Array<int> indices(ne);
  ParallelFor (ne, [&] (size_t edge)
               {
                 indices[edge] = edge;
               });

  ngstd::SampleSortI(edge_collapse_weight, indices);
  Array<int> invindices(ne);
  ParallelFor (ne, [&] (size_t edge)
               {
                 invindices[indices[edge]] = edge;
               });
  
  
  static Timer Tdist1sorted("Dist1 Sorted Collapsing");
  Tdist1sorted.Start();

  Array<bool> vertex_collapse(nv);
  Array<bool> edge_collapse(ne);
  edge_collapse = false;
  vertex_collapse = false;

  Array<int> vertex_collapse_to(nv);
  for (size_t i = 0; i < nv; i++)
    vertex_collapse_to[i] = i;
  
  Dist1Collapser collapser(nv, ne);
  Array<Edge> edges(ne);

  ParallelFor (ne, [&] (size_t edge)
               {
                 edges[edge] = Edge(edge, edge_to_vertices[edge][0], edge_to_vertices[edge][1]);
               });

  TableCreator<int> v2e_creator(nv);
  for ( ; !v2e_creator.Done(); v2e_creator++)
    ParallelFor (ne, [&] (size_t e)
      {
        for (int j = 0; j < 2; j++)
          v2e_creator.Add (edge_to_vertices[e][j], e);
      });
  Table<int> v2e = v2e_creator.MoveTable();

  ParallelFor (v2e.Size(), [&] (size_t vnr)
               {
		 // QuickSortI (edge_collapse_weight, v2e[vnr]);
		 QuickSortI (invindices, v2e[vnr]);
               }, TasksPerThread(5));
  
  // build edge dependency
  TableCreator<int> edge_dag_creator(ne);
  for ( ; !edge_dag_creator.Done(); edge_dag_creator++)  
    ParallelFor (v2e.Size(), [&] (size_t vnr)
      {
        auto vedges = v2e[vnr];
        for (int j = 0; j+1 < vedges.Size(); j++)
	  edge_dag_creator.Add (vedges[j+1], vedges[j]);
      }, TasksPerThread(5));
  Table<int> edge_dag = edge_dag_creator.MoveTable();

  RunParallelDependency (edge_dag,
                         [&] (int edgenr)
                         {
                           auto v0 = edge_to_vertices[edgenr][0];
                           auto v1 = edge_to_vertices[edgenr][1];
                           if (edge_collapse_weight[edgenr] >= 0.01 && !vertex_collapse[v0] && !vertex_collapse[v1])
                             {
                               edge_collapse[edgenr] = true;
                               vertex_collapse[v0] = true;
                               vertex_collapse[v1] = true;
                               auto minv = min(v0,v1);
                               vertex_collapse_to[v0] = minv;
                               vertex_collapse_to[v1] = minv;
                             }
                             /*
                           auto edge = edges[edgenr];
                           if (edge_collapse_weight[edge.id] >= 0.01 && !collapser.AnyVertexCollapsed(edge))
                             collapser.CollapseEdge(edge);
                             */
                         });
  Tdist1sorted.Stop();

  for (int i = 0; i < nv; i++)
    vertex_collapse[i] = (vertex_collapse_to[i] != i);

  
  /*
  static Timer tvcoll("Interm Vcollapse");
  tvcoll.Start();

  for (int vert = 0; vert < nv; ++vert) {
    vertex_collapse[vert] = (collapser.GetCollapsedToVertex(vert) != vert);
  }
  tvcoll.Stop();

  static Timer tecoll("Interm Ecollapse");
  tecoll.Start();
  for (auto edge : edges) {
    edge_collapse[edge.id] = collapser.IsEdgeCollapsed(edge);
  }
  tecoll.Stop();
  */
  
  /*
  *testout << "edge_weights = " << edge_collapse_weight << endl;  
  *testout << "edge_collapse = " << edge_collapse << endl;

  {
    int cnt = 0;
    for (bool b : edge_collapse)
      if (b) cnt++;
    *testout << "marked: " << cnt << endl;
  }
  */
  
  int nr_coarse_vertices = ComputeFineToCoarseVertex(
      edge_to_vertices, nv, edge_collapse, vertex_collapse, vertex_coarse);

  ComputeFineToCoarseEdge(edge_to_vertices, vertex_coarse, edge_coarse, coarse_edge_to_vertices);

  ComputeCoarseWeightsEdges(
      edge_to_vertices, coarse_edge_to_vertices, edge_coarse, weights_edges, weights_edges_coarse);

  ComputeCoarseWeightsVertices(edge_to_vertices, vertex_coarse, nr_coarse_vertices, weights_edges,
      weights_vertices, weights_vertices_coarse);

  static Timer Tblock_table("H1-AMG::BlockJacobiTable");
  Tblock_table.Start();

  Array<int> nentries(nr_coarse_vertices);
  nentries = 0;

  for (auto fvert = 0; fvert < nv; ++ fvert) {
    auto cvert = vertex_coarse[fvert];
    if (cvert != -1) {
      nentries[cvert] += 1;
    }
  }

  auto blocks = make_shared<Table<int>>(nentries);
  Array<int> cnt(nr_coarse_vertices);
  cnt = 0;

  for (auto fvert = 0; fvert < nv; ++ fvert) {
    auto cvert = vertex_coarse[fvert];
    if (cvert != -1) {
      (*blocks)[cvert][cnt[cvert]++] = fvert;
    }
  }
  auto bjacobi = sysmat->CreateBlockJacobiPrecond(blocks, 0, 1, free_dofs);
  Tblock_table.Stop();

  SPtrSMdbl prol;
  int level_diff = h1_options.maxlevel - h1_options.level;
  if (h1_options.smoothed && level_diff % h1_options.special_level == 0) {
    prol = H1SmoothedProl(
        vertex_coarse, nr_coarse_vertices, edge_to_vertices, weights_edges, sysmat->IsComplex());
  } else if (h1_options.semi_smoothed && level_diff % h1_options.special_level == 0) {
    auto triv_prol = CreateProlongation(vertex_coarse, nr_coarse_vertices, sysmat->IsComplex());
    prol = CreateSmoothedProlongation(edge_to_vertices, weights_edges, nv, move(triv_prol));
  }
  else {
    prol = CreateProlongation(vertex_coarse, nr_coarse_vertices, sysmat->IsComplex());
  }

  // build coarse mat
  static Timer Trestrict_sysmat("H1-AMG::RestrictSysmat");
  Trestrict_sysmat.Start();
  auto coarsemat = dynamic_pointer_cast<SparseMatrixTM<double>>(sysmat->Restrict(*prol));
  Trestrict_sysmat.Stop();

  cout << IM(5) << "H1 level " << h1_options.level
       << ", Nr. vertices: " << nv << ", Nr. edges: " << ne << endl
       << "e/v: " << ne/double(nv) << endl
       << "coarse/fine verts: " << nr_coarse_vertices/double(nv)
       << ", coarse/fine edges: " << coarse_edge_to_vertices.Size()/double(ne)<< endl;

  int smoother_its = 1;
  if (h1_options.variable_vcycle) { smoother_its = pow(2, level_diff); }

  auto amg_h1 = make_shared<H1AMG_Mat>(sysmat, bjacobi, std::move(prol), smoother_its);

  if (nr_coarse_vertices <= h1_options.min_verts
      || nr_coarse_vertices >= h1_options.vertex_factor * nv
      || h1_options.level <= 1)
  {
    auto sptr_cmat = shared_ptr<BaseSparseMatrix>(coarsemat);
    sptr_cmat->SetInverseType (SPARSECHOLESKY);
    auto inv = sptr_cmat->InverseMatrix();
    amg_h1->SetRecursive(inv, sptr_cmat);
  }
  else {
    auto new_options = H1Options(h1_options);
    new_options.level = h1_options.level-1;
    // pass NULL, because all non-free dofs should be collapsed to ground by now
    auto recAMG = BuildH1AMG(
        coarsemat, coarse_edge_to_vertices, weights_edges_coarse, weights_vertices_coarse, nullptr,
        new_options);
    amg_h1->SetRecursive(recAMG);
  }

  return amg_h1;
}

UPtrSMdbl CreateProlongation(const Array<int>& vertex_coarse, int ncv, bool complx)
{
  static Timer Tcreate_prol("H1-AMG::CreateProlongation");
  RegionTimer Rcreate_prol(Tcreate_prol);

  auto nv = vertex_coarse.Size();
  Array<int> non_zero_per_row(nv);
  non_zero_per_row = 0;

  for (auto i = 0; i < nv; ++i) {
    if (vertex_coarse[i] != -1) { non_zero_per_row[i] = 1; }
  }

  UPtrSMdbl prol = nullptr;
  if (!complx) {
    prol = UPtrSMdbl(new SparseMatrix<double>(non_zero_per_row, ncv));
  } else {
    prol = UPtrSMdbl(new SparseMatrix<double, Complex, Complex>(non_zero_per_row, ncv));
  }

  for (auto i = 0; i < nv; ++i) {
    if (vertex_coarse[i] != -1) {
      (*prol)(i, vertex_coarse[i]) = 1;
    }
  }

  return move(prol);
}

}  // h1amg
