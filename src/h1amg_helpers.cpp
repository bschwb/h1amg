#include <bla.hpp>

#include "h1amg_helpers.hpp"

namespace h1amg
{

double SplitMatrixConstant(
    ngbla::FlatMatrix<double> source_matrix, ngbla::FlatMatrix<double>& constant_matrix,
    ngbla::FlatMatrix<double>& nullspace_matrix)
{
  const int nr_dofs = source_matrix.Height();

  double sum_matrix = 0;
  for (auto val : source_matrix.AsVector()) {
    sum_matrix += val;
  }
  double c = sum_matrix / (nr_dofs * nr_dofs);

  const double TOL = 1e-10;
  if (abs(c) <= TOL) c = 0;

  constant_matrix = c;
  nullspace_matrix = source_matrix - constant_matrix;

  return c;
}

}  // h1amg
