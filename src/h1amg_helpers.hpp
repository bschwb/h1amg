#ifndef H1AMG_HELPERS_HPP_
#define H1AMG_HELPERS_HPP_

namespace h1amg {

// Splits given matrix in one with constant functions in nullspace
// and a constant part.
// Returns a constant, such that constant_matrix = c * one_matrix
double SplitMatrixConstant(
    ngbla::FlatMatrix<double> source_matrix, ngbla::FlatMatrix<double>& constant_matrix,
    ngbla::FlatMatrix<double>& nullspace_matrix);

}  // h1amg

#endif  // H1AMG_HELPERS_HPP_
