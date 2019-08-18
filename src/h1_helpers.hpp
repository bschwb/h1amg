#ifndef H1AMG_MAT_HELPERS_HPP_
#define H1AMG_MAT_HELPERS_HPP_

#include <ngstd.hpp>

namespace h1amg
{
  using ngstd::INT;

// Figure out best edges to collapse or vertices to collapse with ground
// with iterative improvement.
// Given is a edge to vertices mapping, collapse weights for the edges and
// vertices.
void IterativeCollapse(
    const ngstd::Array<INT<2>>& edge_to_vertices, const ngstd::Array<double>& edge_collapse_weight,
    const ngstd::Array<double>& vertex_collapse_weight, std::shared_ptr<ngstd::BitArray> free_dofs,
    ngstd::Array<bool>& edge_collapse, ngstd::Array<bool>& vertex_collapse);

// Compute the balanced collapse weights for edges and vertices from the
// given weights.
// All collapse weights added up should amount to the number of vertices.
void ComputeCollapseWeights(
    const ngstd::Array<INT<2>>& edge_to_vertices, const ngstd::Array<double>& weights_edges,
    const ngstd::Array<double>& weights_vertices, ngstd::Array<double>& vertex_strength,
    ngstd::Array<double>& edge_collapse_weight, ngstd::Array<double>& vertex_collapse_weight);

// Compute how fine vertices get mapped to coarse vertices, given which edges
// and vertices to collapse.
// -1 in vertex_coarse means vertex maps to ground.
// Returns number of coarse vertices.
int ComputeFineToCoarseVertex(
    const ngstd::Array<INT<2>>& edge_to_vertices, int nverts,
    const ngstd::Array<bool>& edge_collapse, const ngstd::Array<bool>& vertex_collapse,
    ngstd::Array<int>& vertex_coarse);

// Compute how fine edges map to coarse edges, using the alreade computed
// vertex mapping.
void ComputeFineToCoarseEdge(
    const ngstd::Array<INT<2>>& edge_to_vertices, const ngstd::Array<int>& vertex_coarse,
    ngstd::Array<int>& edge_coarse, ngstd::Array<INT<2>>& coarse_edge_to_vertices);


// Computes coarse edge weights to use for the next iteration.
void ComputeCoarseWeightsEdges(
  const ngstd::Array<INT<2>>& edge_to_vertices,
  const ngstd::Array<INT<2>>& coarse_edge_to_vertices, const ngstd::Array<int>& edge_coarse,
  const ngstd::Array<double>& weights_edges, ngstd::Array<double>& weights_edges_coarse);

// Compute coarse vertex weights to use for the next iteration.
void ComputeCoarseWeightsVertices(
  const ngstd::Array<INT<2>>& edge_to_vertices, const ngstd::Array<int>& vertex_coarse,
  const int nr_coarse_vertices, const ngstd::Array<double>& weights_edges,
  const ngstd::Array<double>& weights_vertices, ngstd::Array<double>& weights_vertices_coarse);

}  // h1amg

#endif  // H1AMG_MAT_HELPERS_HPP_
