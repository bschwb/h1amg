#ifndef H1AMG_H1_HPP_
#define H1AMG_H1_HPP_

#include <memory>

#include <la.hpp>
#include <ngstd.hpp>


namespace h1amg
{

using SPtrBJacobi = std::shared_ptr<ngla::BaseBlockJacobiPrecond>;
using UPtrSMdbl = std::unique_ptr<ngla::SparseMatrixTM<double>>;

// Provides a matrix which acts as an AMG preconditioner for H1 problems of
// order 1.
class H1AMG_Mat : public ngla::BaseMatrix
{
public:
  // Preconditions: * a_system and a_jacobi_smoother act as square matrices and
  //                  have the same dimensions.
  //                * height of a_prolongation is the same as width of a_system
  //                  and width of a_prolongation is smaller
  H1AMG_Mat(
      const ngla::BaseSparseMatrix* a_system, SPtrBJacobi a_bjacobi,
      UPtrSMdbl a_prolongation, const int a_smoother_its);

  virtual ~H1AMG_Mat() override
  { }

  // Preconditions: * a_recursive acts as a square matrix and its dimensions
  //                  are the same as width of t_prolongation
  void SetRecursive(
      std::shared_ptr<ngla::BaseMatrix> a_recursive,
      std::shared_ptr<ngla::BaseSparseMatrix> a_cmat = nullptr);

  virtual bool IsComplex() const override
  { return m_system->IsComplex(); }

  virtual int VHeight() const override
  { return m_system->Height(); }

  virtual int VWidth() const override
  { return m_system->Width(); }

  // Return y as the result of the preconditioner applied to x.
  virtual void Mult(const ngla::BaseVector& b, ngla::BaseVector& x) const override;

  virtual ngla::AutoVector CreateVector() const override
  { return m_system->CreateVector(); }

private:
  const ngla::BaseSparseMatrix* m_system = nullptr;
  std::shared_ptr<ngla::BaseSparseMatrix> m_cmat = nullptr;

  SPtrBJacobi m_bjacobi = nullptr;

  UPtrSMdbl m_prolongation = nullptr;
  UPtrSMdbl m_trans_prolongation = nullptr;
  std::shared_ptr<ngla::BaseMatrix> m_recursive = nullptr;

  const int m_smoother_its;
};

}  // h1amg

#endif  // H1AMG_H1_HPP_
