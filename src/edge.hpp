#ifndef MYAMG_EDGE_HPP_
#define MYAMG_EDGE_HPP_

#include <ostream>

namespace h1amg
{
using Vertex = size_t;

struct Edge {
  Edge()
  { }

  explicit Edge(std::size_t a_id, Vertex a_v1, Vertex a_v2)
    : id(a_id), v1(a_v1), v2(a_v2)
  { }

  std::size_t id;
  Vertex v1;
  Vertex v2;
};


inline bool operator<(const Edge& lhs, const Edge& rhs)
{ return lhs.id < rhs.id; }

inline bool operator==(const Edge& lhs, const Edge& rhs)
{ return (lhs.id == rhs.id) && (lhs.v1 == rhs.v1) && (lhs.v2 == rhs.v2); }

inline std::ostream& operator<<(std::ostream& str, const Edge& edge)
{
  str << "Edge " << edge.id << ": (" << edge.v1 << ", " << edge.v2 << ")";
  return str;
}



}  // myamg

#endif  // MYAMG_EDGE_HPP_
