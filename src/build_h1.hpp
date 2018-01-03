#ifndef H1AMG_BUILDH1_HPP_
#define H1AMG_BUILDH1_HPP_

#include <ngstd.hpp>
#include <la.hpp>

#include "h1.hpp"

namespace h1amg
{

struct H1Options {
  int maxlevel = 10;
  int level = 10;
  bool variable_vcycle = false;
  float vertex_factor = 0.8;
  int min_verts = 1000;
  bool smoothed = false;
  bool semi_smoothed = true;
  int special_level = 3;

  H1Options()
  { }

  H1Options(int a_maxlevel, int a_level, bool a_variable_vcycle)
    : maxlevel(a_maxlevel), level(a_level), variable_vcycle(a_variable_vcycle)
  { }
};

// Return a algebraic multigrid preconditioner for H1 problems
// Defaults:
//   * 10 levels
//   * sparsecholesky inverse at lowest level
//   * vertex_factor = 0.8
//   * min_verts = 20
//   * special prolongation every 3rd level
std::shared_ptr<H1AMG_Mat> BuildH1AMG(
    shared_ptr<ngla::SparseMatrixTM<double>> sysmat, const ngstd::Array<INT<2>>& edge_to_vertices,
    const ngstd::Array<double>& weights_edges, const ngstd::Array<double>& weights_vertices,
    shared_ptr<ngstd::BitArray> free_dofs, const H1Options& h1_options);


}  // h1amg

#endif  // H1AMG_BUILDH1_HPP_
