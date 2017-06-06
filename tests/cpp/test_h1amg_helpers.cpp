#include <catch.hpp>

#include <comp.hpp>

#include <h1amg_helpers.hpp>
using namespace h1amg;

TEST_CASE("Split Matrix Constant")
{
  ngbla::Matrix<double> constant_matrix;
  ngbla::Matrix<double> nullspace_matrix;
  ngbla::Matrix<double> source_matrix;

  SECTION("Split Identity Matrix")
  {
    ngbla::Matrix<double> identity = ngbla::Identity(3);
    source_matrix = 3*identity;
    REQUIRE(1 == SplitMatrixConstant(source_matrix, constant_matrix, nullspace_matrix));
  }

  SECTION("Split Const Nullspace Matrix")
  {
    ngbla::Matrix<double> source_matrix =
        {{ 1, 0, -1},
         { 0, 0,  0},
         {-1, 0,  1}};

    REQUIRE(0 == SplitMatrixConstant(source_matrix, constant_matrix, nullspace_matrix));
  }

}
