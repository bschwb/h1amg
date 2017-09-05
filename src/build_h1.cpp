#include <ngstd.hpp>
using namespace ngstd;
#include <la.hpp>
using namespace ngla;

#include "h1.hpp"
#include "h1_helpers.hpp"

#include "h1_smoothed_prol.hpp"
#include "dist1collapser.hpp"
#include "edge.hpp"
#include "sample_sort.hpp"

#include "build_h1.hpp"

namespace h1amg
{

using UPtrSMdbl = unique_ptr<SparseMatrixTM<double>>;

UPtrSMdbl CreateProlongation(const Array<int>& vertex_coarse, int ncv, bool complx);

shared_ptr<H1AMG_Mat> BuildH1AMG(
    const BaseSparseMatrix* sysmat, const Array<INT<2>>& edge_to_vertices,
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

  static Timer Tdist1sorted("Dist1 Sorted Collapsing");
  Tdist1sorted.Start();
  static Timer t1("Dist1 Sorted Collapsing sorting");
  static Timer t2("Dist1 Sorted Collapsing work");

  Dist1Collapser collapser(nv, ne);

  t1.Start();
  Array<int> indices(ne);
  Array<Edge> edges(ne);
  for (size_t edge = 0; edge < ne; ++edge) {
    indices[edge] = edge;
    edges[edge] = Edge(edge, edge_to_vertices[edge][0], edge_to_vertices[edge][1]);
  }

  SampleSortI(edge_collapse_weight, indices);
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

  static Timer tvcoll("Interm Vcollapse");
  tvcoll.Start();
  Array<bool> vertex_collapse(nv);
  for (int vert = 0; vert < nv; ++vert) {
    vertex_collapse[vert] = (collapser.GetCollapsedToVertex(vert) != vert);
  }
  tvcoll.Stop();

  static Timer tecoll("Interm Ecollapse");
  tecoll.Start();
  Array<bool> edge_collapse(ne);
  for (auto edge : edges) {
    edge_collapse[edge.id] = collapser.IsEdgeCollapsed(edge);
  }
  tecoll.Stop();

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

  UPtrSMdbl prol;
  int level_diff = h1_options.maxlevel - h1_options.level;
  if (h1_options.smoothed && level_diff % h1_options.special_level == 0) {
    prol = H1SmoothedProl(
        vertex_coarse, nr_coarse_vertices, edge_to_vertices, weights_edges, sysmat->IsComplex());
  } else if (h1_options.semi_smoothed && level_diff % h1_options.special_level == 0) {
    prol = H1SemiSmoothedProl(
        vertex_coarse, nr_coarse_vertices, edge_to_vertices, weights_edges, sysmat->IsComplex());
  }
  else {
    prol = CreateProlongation(vertex_coarse, nr_coarse_vertices, sysmat->IsComplex());
  }

  // build coarse mat
  static Timer Trestrict_sysmat("H1-AMG::RestrictSysmat");
  Trestrict_sysmat.Start();
  auto coarsemat = sysmat->Restrict(*prol);
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
