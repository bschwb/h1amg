#include <cassert>
#include <stdexcept>
#include <algorithm>
using namespace std;

#include <ngstd.hpp>
using namespace ngstd;

#include "edge.hpp"

#include "dist1collapser.hpp"

namespace h1amg
{

void Dist1Collapser::CollapseEdge(Edge edge)
{
  if (m_edge_collapsed[edge.id] == true) {
    throw std::logic_error("Can't collapse Edge " + to_string(edge.id) + "! Edge is already collapsed.");
  }
  if (m_vertex_collapsed_edge[edge.v1] != nullopt) {
    throw std::logic_error(
      "Can't collapse Edge " + to_string(edge.id) + "! Source Vertex " + to_string(edge.v1) +
      " is already collapsed.");
  }
  if (m_vertex_collapsed_edge[edge.v2] != nullopt) {
    throw std::logic_error(
      "Can't collapse Edge " + to_string(edge.id) + "! Target Vertex " + to_string(edge.v2) +
      " is already collapsed.");
  }

  auto min_vert = min(edge.v1, edge.v2);
  m_vertex_collapsed_to[edge.v1] = min_vert;
  m_vertex_collapsed_to[edge.v2] = min_vert;

  m_vertex_collapsed_edge[edge.v1] = edge;
  m_vertex_collapsed_edge[edge.v2] = edge;

  m_edge_collapsed[edge.id] = true;
}


void Dist1Collapser::UncollapseEdge(Edge edge)
{
  m_vertex_collapsed_to[edge.v1] = edge.v1;
  m_vertex_collapsed_to[edge.v2] = edge.v2;

  m_vertex_collapsed_edge[edge.v1] = nullopt;
  m_vertex_collapsed_edge[edge.v2] = nullopt;

  m_edge_collapsed[edge.id] = false;
}


optional<Edge> Dist1Collapser::GetCollapsedEdge(Edge edge) const
{
  optional<Edge> v1_edge = m_vertex_collapsed_edge[edge.v1];
  optional<Edge> v2_edge = m_vertex_collapsed_edge[edge.v2];

  if (v1_edge == v2_edge) {
    if (v1_edge) { assert(*v1_edge == edge); }
    return v1_edge;
  }
  else if (v1_edge && v2_edge) { return nullopt; }
  else if (v1_edge) { return v1_edge; }
  else { return v2_edge; }

  assert(false && "This can't happen");
}

}  // myamg
