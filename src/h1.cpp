#include <cassert>
#include <memory>
using namespace std;

#include <la.hpp>
using namespace ngla;
#include <ngstd.hpp>
using namespace ngstd;

#include "h1_helpers.hpp"
#include "h1.hpp"


namespace h1amg
{

H1AMG_Mat::H1AMG_Mat(
    const BaseSparseMatrix* a_system, SPtrBJacobi a_bjacobi,
    UPtrSMdbl a_prolongation, const int a_smoother_its)
  : m_system(a_system), m_bjacobi(a_bjacobi),
    m_prolongation(std::move(a_prolongation)),
    m_trans_prolongation(UPtrSMdbl(TransposeMatrix(*m_prolongation))),
    m_smoother_its(a_smoother_its)
{
  if(!m_system) {
    throw Exception("Nullptr to m_system in H1AMG_Mat constructor.");
  }
  if(!m_bjacobi) {
    throw Exception("Nullptr to m_bjacobi in H1AMG_Mat constructor.");
  }
  if(!m_prolongation) {
    throw Exception("Nullptr to m_prolongation in H1AMG_Mat constructor.");
  }

  if(m_system->Width() != m_system->Height()) {
    throw Exception("System matrix is not a square matrix in H1AMG_Mat constructor.");
  }
  if(m_bjacobi->Width() != m_bjacobi->Height()) {
    throw Exception("Jacobi smoother is not a square matrix in H1AMG_Mat constructor.");
  }
  if(m_system->Height() != m_bjacobi->Width()) {
    throw Exception("Dimensions of system matrix and jacobi smoother don't match \
                     in H1AMG_Mat constructor.");
  }

  if(m_prolongation->Height() != m_system->Width()) {
    throw Exception("Dimensions of system matrix and prolongation don't match \
                     in H1AMG_Mat constructor.");
  }
  if( (m_prolongation->Height()>0) && (m_prolongation->Height() <= m_prolongation->Width()) ) {
    throw Exception("Prolongation is not narrowing in H1AMG_Mat constructor.");
  }
}

void H1AMG_Mat::SetRecursive(
    shared_ptr<BaseMatrix> a_recursive, shared_ptr<BaseSparseMatrix> a_cmat)
{

  if(!a_recursive) {
    throw Exception("Nullptr to t_recursive in H1AMG_Mat::SetRecursive method.");
  }
  if(a_recursive->Width() != a_recursive->Height()) {
    throw Exception("Recursive matrix is not a square matrix in H1AMG_Mat::SetRecursive.");
  }
  if(a_recursive->Height() != m_prolongation->Width()) {
    cout << "prol width || rec mat height:  " << m_prolongation->Width() << " || "
         << a_recursive->Height() << endl;
    throw Exception("Dimensions of prolongation and recursive matrix don't match in \
                     H1AMG_Mat::SetRecursive.");
  }

  m_recursive = a_recursive;
}

void H1AMG_Mat::Mult(const BaseVector& b, BaseVector& x) const
{
  assert(m_system);
  assert(m_bjacobi);
  assert(m_prolongation);
  assert(m_recursive);

  static Timer Th1_mult("H1-AMG::Mult");
  RegionTimer Rh1_mult(Th1_mult);

  x = 0;
  m_bjacobi->GSSmooth(x, b, m_smoother_its);

  auto residuum = m_system->CreateVector();
  residuum = b - (*m_system) * x;

  auto coarse_residuum = m_recursive->CreateVector();
  coarse_residuum = (*m_trans_prolongation) * residuum;

  auto coarse_x = m_recursive->CreateVector();
  m_recursive->Mult(coarse_residuum, coarse_x);

  x += (*m_prolongation) * coarse_x;

  m_bjacobi->GSSmoothBack (x, b, m_smoother_its);
}

}  // h1amg
