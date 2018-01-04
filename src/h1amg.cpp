#include <ngstd.hpp>
using namespace ngstd;
#include <bla.hpp>
using namespace ngbla;
#include <comp.hpp>
using namespace ngcomp;

#include "complex_mat.hpp"

#include "h1.hpp"
#include "build_h1.hpp"
#include "h1amg_helpers.hpp"

#include "h1amg.hpp"


namespace h1amg
{

  template <typename THASH, typename T, typename FUNC>
  void ParallelIterate (const ngstd::ParallelHashTable<THASH,T> & hashtable, FUNC func)
  {
    Array<size_t> base(hashtable.NumBuckets());
    size_t sum = 0;
    for (size_t i = 0; i < hashtable.NumBuckets(); i++)
      {
        base[i] = sum;
        sum += hashtable.Used(i); 
      }
    ParallelFor(hashtable.NumBuckets(),
                [&] (size_t nr)
                {
                  size_t cnt = base[nr];
                  hashtable.Iterate(nr,
                                    [&cnt, func] (THASH key, T val)
                                    {
                                      func(cnt, key, val);
                                      cnt++;
                                    });
                });
  }



H1AMG::H1AMG(const PDE& a_pde, const Flags& a_flags, const string a_name)
  : H1AMG(a_pde.GetBilinearForm(a_flags.GetStringFlag("bilinearform", "")), a_flags, a_name)
{ }

H1AMG::H1AMG(shared_ptr<BilinearForm> a_bfa, const Flags& a_flags, const string a_name)
  : Preconditioner(a_bfa, a_flags, a_name), bfa(a_bfa), m_ndof(bfa->GetFESpace()->GetNDof()),
  dof_pair_weights(2*m_ndof), m_hashlocks(2*m_ndof)
{
  if (bfa) {
    while (bfa->GetLowOrderBilinearForm()) {
      bfa = bfa->GetLowOrderBilinearForm();
    }
  }
  m_h1_options.maxlevel = int(flags.GetNumFlag("levels", 10));
  m_h1_options.level = int(flags.GetNumFlag("levels", 10));
  m_h1_options.variable_vcycle = flags.GetDefineFlag("variable_vcycle");
  m_h1_options.smoothed = flags.GetDefineFlag("smoothed");
  m_h1_options.semi_smoothed = !flags.GetDefineFlag("not_semi_smoothed");
  m_h1_options.special_level = int(flags.GetNumFlag("special_level", 3));
}

void H1AMG::InitLevel(shared_ptr<BitArray> afreedofs)
{
  static Timer Tinit_level("H1-AMG::InitLevel");
  RegionTimer Rinit_level(Tinit_level);
  *testout << "initlevel amg" << endl;
  freedofs = afreedofs;
  int nr_dofs = freedofs->Size(); //bfa->GetFESpace()->GetNDof();
  weights_vertices = Array<double>(nr_dofs);
  weights_vertices = 0;
}

void H1AMG::FinalizeLevel(const BaseMatrix* mat)
{
  static Timer Tfinlevel("H1-AMG::FinalizeLevel");
  RegionTimer Rfinlevel(Tfinlevel);

  *testout << "finalize lvl  amg" << endl;

  static Timer Tcnt_edges("H1-AMG::FinalizeLevel::CountEdges");
  Tcnt_edges.Start();
  /*
  int cnt = 0;
  for (auto key_val: dof_pair_weights)
  { ++cnt; }
  cout << "cnt " << cnt << " =?= " << par_dof_pair_weights.Used() << endl;
  */
  size_t cnt = par_dof_pair_weights.Used();
  Tcnt_edges.Stop();

  static Timer Tcreate_e2v("H1-AMG::FinalizeLevel::CreateE2VMapping");
  Tcreate_e2v.Start();
  Array<double> weights (cnt);
  Array<INT<2> > edge_to_vertices (cnt);

  /*
  int i = 0;
  for (auto key_val: dof_pair_weights) {
    weights[i] = key_val.second;
    edge_to_vertices[i] = key_val.first;
    ++i;
  }
  */

  ParallelIterate (par_dof_pair_weights,
                   [&weights,&edge_to_vertices] (size_t i, INT<2> key, double weight)
                   {
                     weights[i] = weight;
                     edge_to_vertices[i] = key;
                   });
  /*
  int i = 0;
  par_dof_pair_weights.Iterate
    ( [&weights,&edge_to_vertices,&i] (INT<2> key, double weight)
      {
        weights[i] = weight;
        edge_to_vertices[i] = key;
        ++i;
      } );
  */  

  
  Tcreate_e2v.Stop();

  amg_matrix = BuildH1AMG(
      dynamic_pointer_cast<SparseMatrixTM<double>> (const_cast<BaseMatrix*>(mat)->shared_from_this()),
      edge_to_vertices, weights, weights_vertices,
      freedofs, m_h1_options);

  cout << "matrices done" << endl;

  if (timing) { Timing(); }
  if (test) { Test(); }
}

void H1AMG::AddElementMatrixCommon(
    FlatArray<int> dnums, const FlatMatrix<double>& elmat, LocalHeap& lh)
{
  static Timer addelmat("amg - addelmat",2);
  static Timer addelmat1("amg - addelmat1",2);
  static Timer addelmat2("amg - addelmat2",2);
  static Timer addelmat3("amg - addelmat3",2);

  int tid = TaskManager::GetThreadId();
  ThreadRegionTimer reg(addelmat,tid);
  NgProfiler::StartThreadTimer (addelmat1, tid);
  
  auto ndof = elmat.Height();
  FlatMatrix<double> constant_elmat(ndof, lh);
  FlatMatrix<double> nullspace_elmat(ndof, lh);

  auto vertex_weight = SplitMatrixConstant(elmat, constant_elmat, nullspace_elmat);

  if (vertex_weight != 0.0) {
    vertex_weight /= ndof;
    for (auto i=0; i < ndof; ++i)
      weights_vertices[ dnums[i] ] += vertex_weight;
  }

  FlatMatrix<double> approx_elmat(ndof, lh);
  approx_elmat = 0;

  NgProfiler::StopThreadTimer (addelmat1, tid);

  FlatMatrix<double> schur_complement(2, lh);
  BitArray used(ndof, lh);
  for (auto i=0; i < ndof; ++i) {
    for (auto j=0; j < ndof; ++j) {
      auto first_dof = dnums[i];
      auto second_dof = dnums[j];
      if (first_dof < second_dof) {
        schur_complement = 0;
        used.Clear();

        // the schur-complement is calculated in respect to the two current
        // dofs
        used.Set(i);
        used.Set(j);

        CalcSchurComplement(nullspace_elmat, schur_complement, used, lh);

        auto schur_entry = schur_complement(0);

        if (!std::isnan(schur_entry)) {
          approx_elmat(i, i) += schur_entry;
          approx_elmat(j, j) += schur_entry;
          approx_elmat(i, j) -= schur_entry;
          approx_elmat(j, i) -= schur_entry;

          INT<2> i2(first_dof, second_dof);
          auto hash = HashValue(i2, dof_pair_weights.Size());
          /*
          {
            ThreadRegionTimer reg(addelmat2,tid);            
            lock_guard<mutex> guard(m_hashlocks[hash]);
            dof_pair_weights[i2] += schur_entry;
          }
          */
          {
            ThreadRegionTimer reg(addelmat3,tid);
            par_dof_pair_weights.Do (i2, [schur_entry] (auto & v) { v += schur_entry; });
          }
        }
        else {
          INT<2> i2(first_dof, second_dof);
          dof_pair_weights[i2] += 0;
        }
      }
    }
  }

}

void H1AMG::AddElementMatrix(
    FlatArray<int> dnums, const FlatMatrix<double>& elmat, ElementId ei, LocalHeap& lh)
{
  HeapReset hr(lh);
  AddElementMatrixCommon(dnums, elmat, lh);
}


void H1AMG::AddElementMatrix(
    FlatArray<int> dnums, const FlatMatrix<Complex>& elmat, ElementId ei, LocalHeap& lh)
{
  HeapReset hr(lh);
  auto combined_elmat = AddedParts(elmat, lh);
  AddElementMatrixCommon(dnums, combined_elmat, lh);
}

RegisterPreconditioner<H1AMG> initmyamg("h1amg");

}  // h1amg
