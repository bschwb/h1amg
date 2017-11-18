#ifndef H1AMG_HPP_
#define H1AMG_HPP_

#include <memory>
#include <mutex>

#include <comp.hpp>
#include <la.hpp>

#include "h1.hpp"

namespace h1amg
{

class H1AMG : public ngcomp::Preconditioner
{
private:
  shared_ptr<ngcomp::BilinearForm> bfa;

  shared_ptr<ngla::BaseMatrix> amg_matrix;

  size_t m_ndof;
  ngstd::HashTable<INT<2>,double> dof_pair_weights;
  ngstd::ParallelHashTable<INT<2>,double> par_dof_pair_weights;
  ngstd::Array<double> weights_vertices;

  ngstd::Array<std::mutex> m_hashlocks;

  shared_ptr<ngla::BitArray> freedofs;

  H1Options m_h1_options;

public:
  H1AMG(
      const ngcomp::PDE& a_pde, const ngstd::Flags& a_flags,
      const string a_name = "H1AMG Preconditioner");

  H1AMG(
      shared_ptr<ngcomp::BilinearForm> a_bfa, const Flags& a_flags,
      const string a_name = "H1AMG Preconditioner");

  virtual ~H1AMG() override
  { }

  virtual void Update() override
  { }

  virtual const ngla::BaseMatrix& GetAMatrix() const override
  { return bfa->GetMatrix(); }

  virtual const ngla::BaseMatrix& GetMatrix() const override
  { return *amg_matrix; }

  virtual void InitLevel(shared_ptr<BitArray> afreedofs = nullptr) override;

  virtual void FinalizeLevel(const BaseMatrix* mat) override;

  void AddElementMatrixCommon(
      ngstd::FlatArray<int> dnums, const ngbla::FlatMatrix<double>& elmat, ngstd::LocalHeap& lh);

  virtual void AddElementMatrix(
      ngstd::FlatArray<int> dnums, const ngbla::FlatMatrix<double>& elmat, ngcomp::ElementId ei,
      ngstd::LocalHeap& lh) override;

  virtual void AddElementMatrix(
      ngstd::FlatArray<int> dnums, const ngbla::FlatMatrix<Complex>& elmat, ngcomp::ElementId ei,
      ngstd::LocalHeap& lh) override;

  virtual const char* ClassName() const override
  { return "H1AMG Preconditioner"; }
};

}  // h1amg

#endif  // H1AMG_HPP_
