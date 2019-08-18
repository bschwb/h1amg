#ifndef MYAMG_DIST1COLLAPSER_HPP_
#define MYAMG_DIST1COLLAPSER_HPP_

#include <vector>
#include <tuple>

#include <ngstd.hpp>

#include "edge.hpp"

namespace h1amg
{

class Dist1Collapser
{
public:
  // Construct a Dist1Collapser
  Dist1Collapser(size_t t_nr_vertices, size_t t_nr_edges)
    : m_vertex_collapsed_to(t_nr_vertices),
      m_vertex_collapsed_edge(t_nr_vertices), // , std::nullopt),
      m_edge_collapsed(t_nr_edges) // , false)
  {
    static ngstd::Timer Tcollaps_init("Collapse Init");
    ngstd::RegionTimer Rcollaps_init(Tcollaps_init);

    m_edge_collapsed = false;

    for (size_t vertex_i=0; vertex_i < t_nr_vertices; ++vertex_i) {
      m_vertex_collapsed_to[vertex_i] = vertex_i;
    }
  }

  // Return true if any vertex collapses to an other vertex.
  inline bool AnyVertexCollapsed(Edge edge) const
  {
    return (m_vertex_collapsed_edge[edge.v1] != std::nullopt)
        || (m_vertex_collapsed_edge[edge.v2] != std::nullopt);
  }

  // Return index of vertex which the specified vertex collapses too.
  // If the vertex doesn't collapse or collapses to itself, than this is the
  // same number as the input.
  inline size_t GetCollapsedToVertex(size_t t_vertex_number) const
  { return m_vertex_collapsed_to[t_vertex_number]; }

  // Return edge index of the edge which collapses with one of the vertices.
  // If both vertices collapse to the same edge return that edge
  // (then the method was also called with that edge)
  // If both vertices collapse to different edge or none at all than return std::nullopt.
  std::optional<Edge> GetCollapsedEdge(Edge edge) const;

  // Store information of collapsing.
  void CollapseEdge(Edge edge);

  // Revert collapse informations about edge.
  void UncollapseEdge(Edge edge);

  bool IsEdgeCollapsed(Edge edge)
  {
    return m_edge_collapsed[edge.id];
  }

private:
  Array<size_t> m_vertex_collapsed_to;
  Array<std::optional<Edge>> m_vertex_collapsed_edge;

  Array<bool> m_edge_collapsed;
};

}  // myamg

#endif  // MYAMG_DIST1COLLAPSER_HPP_
