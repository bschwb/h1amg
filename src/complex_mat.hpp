#ifndef H1AMG_COMPLEX_MAT_HPP_
#define H1AMG_COMPLEX_MAT_HPP_

#include <bla.hpp>

namespace h1amg
{

// Construct a FlatMatrix<double> by entrywise adding real and imaginary part
// of a FlatMatrix<Complex>.
ngbla::FlatMatrix<double> AddedParts(
    const ngbla::FlatMatrix<Complex>& complex_mat, ngstd::LocalHeap& lh);

}  // h1amg

#endif  // H1AMG_COMPLEX_MAT_HPP_
