#include <bla.hpp>
using namespace ngbla;

#include "complex_mat.hpp"

namespace h1amg
{

FlatMatrix<double> AddedParts(const FlatMatrix<Complex>& complex_mat, LocalHeap& lh)
{
  auto height = complex_mat.Height();
  auto width = complex_mat.Width();
  FlatMatrix<double> output_mat(height, width, lh);

  for (auto i=0; i < height; ++i) {
    for (auto j=0; j < width; ++j) {
      output_mat(i, j) = complex_mat(i, j).real() + complex_mat(i, j).imag();
    }
  }

  return output_mat;
}

}  // h1amg
