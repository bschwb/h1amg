#include <catch.hpp>

#include <comp.hpp>
#include <array.hpp>

#include <h1_helpers.hpp>
using namespace h1amg;

// ----------------------------------------------------------------------------
// -- ComputeCollapseWeights

TEST_CASE("Only Edge Weights")
{
  ngstd::Array<INT<2>> edge_to_vertices { INT<2>(0, 1) };
  ngstd::Array<double> weights_edges { 1 };
  ngstd::Array<double> weights_vertices { 0, 0 };

  ngstd::Array<double> vertex_strength;
  ngstd::Array<double> edge_collapse_weight;
  ngstd::Array<double> vertex_collapse_weight;
  ngstd::Array<double> edge_collapse_weight_result { 2 };
  ngstd::Array<double> vertex_collapse_weight_result { 0, 0 };

  ComputeCollapseWeights(
      edge_to_vertices, weights_edges, weights_vertices, vertex_strength, edge_collapse_weight,
      vertex_collapse_weight);

  REQUIRE(edge_collapse_weight_result == edge_collapse_weight);
  REQUIRE(vertex_collapse_weight_result == vertex_collapse_weight);
}

TEST_CASE("One Edge With Vertex Weights")
{
  ngstd::Array<INT<2>> edge_to_vertices { INT<2>(0, 1) };
  ngstd::Array<double> weights_edges { 3 };
  ngstd::Array<double> weights_vertices { 1, 1 };

  ngstd::Array<double> vertex_strength;
  ngstd::Array<double> edge_collapse_weight;
  ngstd::Array<double> vertex_collapse_weight;
  ngstd::Array<double> edge_collapse_weight_result { 1.5  };
  ngstd::Array<double> vertex_collapse_weight_result { 0.25, 0.25 };

  ComputeCollapseWeights(edge_to_vertices, weights_edges, weights_vertices,
                         vertex_strength, edge_collapse_weight,
                         vertex_collapse_weight);

  REQUIRE(edge_collapse_weight_result == edge_collapse_weight);
  REQUIRE(vertex_collapse_weight_result == vertex_collapse_weight);
}
