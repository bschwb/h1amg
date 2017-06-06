#include <ngstd.hpp>
using namespace ngstd;
#include <la.hpp>
using namespace ngla;

#include "h1.hpp"
#include "h1_helpers.hpp"

#include "h1_smoothed_prol.hpp"

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
  cout << IM(5) << "H1 Sysmat nze: " << sysmat->NZE() << endl;
  cout << IM(5) << "H1 Sysmat nze per row: " << sysmat->NZE() / (double)sysmat->Height() << endl;
  auto nr_edges = edge_to_vertices.Size();
  auto nr_vertices = weights_vertices.Size();

  cout << IM(5) << "H1 level " << h1_options.level
       << ", nr_edges = " << nr_edges
       << ", nr_vertices = " << nr_vertices << endl;

  Array<int> vertex_coarse;

  Array<bool> edge_collapse;
  Array<bool> vertex_collapse;

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

  IterativeCollapse(
      edge_to_vertices, edge_collapse_weight, vertex_collapse_weight, free_dofs, edge_collapse,
      vertex_collapse);

  int nr_coarse_vertices = ComputeFineToCoarseVertex(
      edge_to_vertices, vertex_collapse.Size(), edge_collapse, vertex_collapse, vertex_coarse);

  ComputeFineToCoarseEdge(edge_to_vertices, vertex_coarse, edge_coarse, coarse_edge_to_vertices);

  ComputeCoarseWeightsEdges(
      edge_to_vertices, coarse_edge_to_vertices, edge_coarse, weights_edges, weights_edges_coarse);

  ComputeCoarseWeightsVertices(edge_to_vertices, vertex_coarse, nr_coarse_vertices, weights_edges,
      weights_vertices, weights_vertices_coarse);

  Array<int> nentries(nr_coarse_vertices);
  nentries = 0;

  for (auto fvert = 0; fvert < nr_vertices; ++ fvert) {
    auto cvert = vertex_coarse[fvert];
    if (cvert != -1) {
      nentries[cvert] += 1;
    }
  }

  auto blocks = make_shared<Table<int>>(nentries);
  Array<int> cnt(nr_coarse_vertices);
  cnt = 0;

  for (auto fvert = 0; fvert < nr_vertices; ++ fvert) {
    auto cvert = vertex_coarse[fvert];
    if (cvert != -1) {
      (*blocks)[cvert][cnt[cvert]++] = fvert;
    }
  }
  auto bjacobi = sysmat->CreateBlockJacobiPrecond(blocks, 0, 1, free_dofs);

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
  auto coarsemat = sysmat->Restrict(*prol);
  cout << IM(5) << "H1 Coarsemat nze: " << coarsemat->NZE() << endl;
  cout << IM(5) << "H1 Coarsemat nze per row: "
       << coarsemat->NZE() / (double)coarsemat->Height() << endl;

  int smoother_its = 1;
  if (h1_options.variable_vcycle) { smoother_its = pow(2, level_diff); }

  auto amg_h1 = make_shared<H1AMG_Mat>(sysmat, bjacobi, std::move(prol), smoother_its);

  if (nr_coarse_vertices <= h1_options.min_verts
      || nr_coarse_vertices >= h1_options.vertex_factor * nr_vertices
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
