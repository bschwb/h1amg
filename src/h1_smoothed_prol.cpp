#include <cassert>
#include <memory>
using namespace std;

#include <la.hpp>
using namespace ngla;
#include <ngstd.hpp>
using namespace ngstd;

#include "h1_smoothed_prol.hpp"

namespace h1amg
{

Table<int> Coarse2FineVertexTable(const Array<int>& vertex_coarse, int ncv);
SparseMatrix<double> EdgeConnectivityMatrix(const Array<INT<2>>& e2v, int nv);

// Build matrix graph before computing elmats
UPtrSMdbl H1SmoothedProl(
    const Array<int>& vertex_coarse, int ncv, const Array<INT<2>>& e2v, const Array<double>& ew,
    bool complx)
{
  int nverts = vertex_coarse.Size();
  int nedges = e2v.Size();

  auto c2f = Coarse2FineVertexTable(vertex_coarse, ncv);
  auto econ = EdgeConnectivityMatrix(e2v, nverts);
  const SparseMatrix<double> & cecon(econ);

  DynamicTable<int> mat_graph(nverts);

  Array<Array<int>> row_dofs(ncv);
  Array<Array<int>> col_dofs(ncv);

  Array<int> c_2_lc(ncv); //coarse to loc coarse
  c_2_lc = -1;
  Matrix<double> m(2);

  Array<int> ext_d(150);
  Array<int> ext_dc(150);
  Array<Array<int> > exd(2);
  exd[0].SetSize(150);
  exd[1].SetSize(150);
  Array<int> all_coarse(50);
  Array<int> ext_e(150);
  Array<Array<int> > v_patch_es(2);
  v_patch_es[0].SetSize(150);
  v_patch_es[1].SetSize(150);
  Array<int> int_e(50);

  Array<double> mat_2_mem(10000);
  Array<double> mat_3_mem(10000);

  Array<INT<2> > carry_over;

  //Build matrix graph
  for (auto patch : Range(ncv)) {

    if (c2f[patch][0] != -1 && c2f[patch][1] != -1) {
      ext_dc.SetSize(0);
      for (auto k : Range(2)) {
        auto rid = econ.GetRowIndices(c2f[patch][k]);
        for(auto d2 : rid) {
          if(c2f[patch][1-k] != d2) {
            if ((!ext_dc.Contains(vertex_coarse[d2])) && (vertex_coarse[d2] != -1)) {
              ext_dc.Append(vertex_coarse[d2]);
            }
          }
        }
      }
      if (ext_dc.Size()==0) {
        for (auto l:Range(2)) {
          if (c2f[patch][1-l]!=-1) {
            mat_graph.Add(c2f[patch][1-l], patch);
            carry_over.Append(INT<2>(c2f[patch][1-l], patch));
          }
        }
      }
      else {
        for (auto d:ext_dc) {
          mat_graph.Add(c2f[patch][0], d);
          mat_graph.Add(c2f[patch][1], d);
        }
        mat_graph.Add(c2f[patch][0], patch);
        mat_graph.Add(c2f[patch][1], patch);

        col_dofs[patch].SetSize(ext_dc.Size()+1);

        for (auto j : Range(ext_dc.Size())) {
          col_dofs[patch][j] = ext_dc[j];
        }

        col_dofs[patch][ext_dc.Size()] = patch;

        row_dofs[patch].SetSize(2);
        row_dofs[patch][0] = c2f[patch][0];
        row_dofs[patch][1] = c2f[patch][1];
      }
    }
    else {
      for (auto l : Range(2)) {
        if (c2f[patch][1-l] != -1) {
          mat_graph.Add(c2f[patch][1-l],patch);
          carry_over.Append(INT<2>(c2f[patch][1-l],patch));
        }
      }
    }
  }


  Array<int> nrd(nverts);
  for (auto k : Range(nverts)) {
    QuickSort(mat_graph[k]);
    nrd[k] = mat_graph[k].Size();
  }

  UPtrSMdbl prol = nullptr;
  if (!complx) {
    prol = UPtrSMdbl(new SparseMatrix<double>(nrd, ncv));
  } else {
    prol = UPtrSMdbl(new SparseMatrix<double, Complex, Complex>(nrd, ncv));
  }

  prol->AsVector() = 0.0;

  for (auto r : Range(nverts)) {
    auto cs = prol->GetRowIndices(r);
    for (auto k : Range(cs.Size())) {
      cs[k] = mat_graph[r][k];
    }
  }

  for (auto t : carry_over) {
    (*prol)(t[0],t[1]) = 1.0;
  }

  //build and add elmats
  for (auto patch : Range(ncv)) {
    if (c2f[patch][0] != -1 && c2f[patch][1] != -1) {
      auto d_patch = c2f[patch];
      int np = d_patch.Size();

      ext_d.SetSize(0);
      ext_dc.SetSize(0);
      exd[0].SetSize(0);
      exd[1].SetSize(0);
      for (auto k : Range(2)) {
        auto d = d_patch[k];
        auto rid = econ.GetRowIndices(d);

        for (auto d2 : rid) {
          if ((!d_patch.Contains(d2)) && (vertex_coarse[d2] != -1))
          {
            ext_d.Append(d2);
            auto vcd2 = vertex_coarse[d2];

            exd[k].Append(vcd2);

            if (!ext_dc.Contains(vcd2)) {
              ext_dc.Append(vcd2);
            }
          }
        }
      }
      if (ext_dc.Size()!=0) {
        int nexf = ext_d.Size();
        int nexc = ext_dc.Size();

        all_coarse.SetSize(0);
        for (auto d : ext_dc) {
          all_coarse.Append(d);
        }
        all_coarse.Append(patch);

        //no need to sort but also not expensive
        //if you put this back in also put back the other one!
        //QuickSort(all_coarse);
        for (auto k:Range(all_coarse.Size())) {
          c_2_lc[all_coarse[k]] = k;
        }


        int_e.SetSize(0);
        for (auto d : d_patch) {
          for(auto d2 : d_patch) {
            if (d != d2) {
              if (!int_e.Contains(cecon(d,d2))) {
                int_e.Append(cecon(d,d2));
              }
            }
          }
        }

        ext_e.SetSize(0);
        v_patch_es[0].SetSize(0);
        v_patch_es[1].SetSize(0);
        for (auto k : Range(2))
        {
          auto d = d_patch[k];
          auto ods  = econ.GetRowIndices(d);
          auto enrs = econ.GetRowValues(d);
          for (auto j:Range(ods.Size()) ) {
            if ((!int_e.Contains(enrs[j])) && (vertex_coarse[ods[j]] != -1)) {
              v_patch_es[k].Append(enrs[j]);
              if (!ext_e.Contains(enrs[j])) {
                ext_e.Append(enrs[j]);
              }
            }
          }
        }
        int nfic = ext_e.Size();

        m = ew[int_e[0]];
        for (auto k : Range(d_patch.Size())) {
          for (auto j : Range(v_patch_es[k].Size())) {
            m(1-k,1-k) += 2*ew[((v_patch_es[k]))[j]];
          }
        }
        m *= 1.0/(m(0,0)*m(1,1)-m(1,0)*m(0,1));

        // np x nfic  nfic x nexc+1
        FlatMatrix<double> m2(d_patch.Size(), nexc+1, &(mat_2_mem[0]));
        m2 = 0.0;

        for (auto k : Range(v_patch_es.Size())) {
          for(auto j : Range(v_patch_es[k].Size())) {
            m2(k,c_2_lc[((exd[k]))[j]]) += ew[((v_patch_es[k]))[j]];
            m2(k,c_2_lc[vertex_coarse[d_patch[0]]]) += ew[((v_patch_es[k]))[j]];
          }
        }
        FlatMatrix<double> m3(d_patch.Size(), nexc+1, &(mat_3_mem[0]));
        m3 = m*m2;

        prol->AddElementMatrix(row_dofs[patch], col_dofs[patch], m3);
      }
    }
  }

  return std::move(prol);
}


Table<int> Coarse2FineVertexTable(const Array<int>& vertex_coarse, int ncv)
{
  static Timer Tsemi_c2fvt("H1SmoothedProl - Coarse2FineVertexTable");
  RegionTimer Rsemi_c2fvt(Tsemi_c2fvt);

  auto nv = vertex_coarse.Size();

  // two coarse vertex per fine vertex
  Array<int> twos(ncv);
  twos = 2;
  Table<int> c2f(twos);

  // Set all entries of table to -1 initial entry
  // by iterating over rows (FlatArrays can be assigned a scalar value and every entry gets that
  // value)
  for (auto r : c2f) {
    r = -1;
  }

  // invert mapping
  static Timer Tsemi_invmap("H1SmoothedProl - c2fvt - Invert Mapping");
  Tsemi_invmap.Start();
  for (const auto k : Range(nv)) {
    if (vertex_coarse[k] != -1) {
      c2f[vertex_coarse[k]][(c2f[vertex_coarse[k]][0]==-1)?0:1] = k;
    }
  }
  Tsemi_invmap.Stop();

  // sort entries
  static Timer Tsemi_sort("H1SmoothedProl - c2fvt - Sort Entries");
  Tsemi_sort.Start();
  for (auto r : c2f) {
    if (r[0]>r[1]) {
      auto d = r[0];
      r[0] = r[1];
      r[1] = d;
    }
  }
  Tsemi_sort.Stop();

  return c2f;
}

// its a square matrix
// econ(i, j) = #e for which e=(vert_i, vert_j)
// so rows and cols correspond to vertices and the entry is the edge number
SparseMatrix<double> EdgeConnectivityMatrix(const Array<INT<2>>& e2v, int nv)
{
  static Timer Tsemi_edge_con("H1SmoothedProl - EdgeConnectivityMatrix");
  RegionTimer Rsemi_edge_con(Tsemi_edge_con);
  auto ne = e2v.Size();

  // count nr of connected edges of vertex
  static Timer Tsemi_cnt_con("H1SmoothedProl - edgeconmat- Cnt Connected edges");
  Tsemi_cnt_con.Start();

  Array<int> econ_s(nv);
  econ_s = 0;

  for (auto k : Range(ne)) {
    econ_s[e2v[k][0]]++;
    econ_s[e2v[k][1]]++;
  }
  Tsemi_cnt_con.Stop();

  // build the table for the connections
  static Timer Tsemi_con_table("H1SmoothedProl - edgeconmat - Connection table");
  Tsemi_con_table.Start();
  Table<int> tab_econ(econ_s);
  econ_s = 0; // used again for counting!!!!
  for (auto k : Range(ne)) {
    tab_econ[e2v[k][0]][econ_s[e2v[k][0]]++] = k;
    tab_econ[e2v[k][1]][econ_s[e2v[k][1]]++] = k;
  }
  Tsemi_con_table.Stop();

  // build the edge connection matrix
  static Timer Tsemi_con_matrix("H1SmoothedProl - edgeconmat - create matrix");
  Tsemi_con_matrix.Start();
  SparseMatrix<double> econ(econ_s, nv);
  econ.AsVector() = -1;
  for (auto k:Range(nv)) {
    auto ind = econ.GetRowIndices(k);
    auto val = econ.GetRowValues(k);
    Array<int> ods(ind.Size());
    int cnt = 0;

    for (auto j:Range(econ_s[k])) {
      auto enr = tab_econ[k][j];
      auto d = (e2v[enr][0]!=k)?e2v[enr][0]:e2v[enr][1];
      ods[cnt++] = d;
    }

    Array<int> indices(ind.Size());
    for (auto k : Range(ind.Size())) { indices[k]=k; }
    QuickSortI(ods,indices);

    for (auto j : Range(econ_s[k])) {
      ind[j] = ods[indices[j]];
      val[j] = tab_econ[k][indices[j]];
    }
  }
  Tsemi_con_matrix.Stop();

  // Sanity checks
  static Timer Tsemi_con_checks("H1SmoothedProl - edgeconmat - Checks");
  Tsemi_con_checks.Start();
  for (auto k : Range(econ.Height())) {
    auto vs = econ.GetRowValues(k);
    for (auto d : econ.GetRowIndices(k)) {
      if (d<0) {
        cout << endl << endl << "STOOOOOOP2" << endl << endl;
      }
    }
    for (auto j : Range(vs.Size())) {
      if (vs[j]==-1) {
        cout << endl << endl << "STOOOOOOOP" << endl << endl;
      }
    }
  }
  Tsemi_con_checks.Stop();

  return econ;
}

UPtrSMdbl BuildOffDiagSubstitutionMatrix(
    const Array<INT<2>>& e2v, const Array<double>& eweights, int nv);

UPtrSMdbl CreateSmoothedProlongation(
    const Array<INT<2>>& e2v, const Array<double>& eweights, int nv,
    const UPtrSMdbl triv_prol)
{
  static Timer Tsmoothed("H1 Create Smoothed Prolongation");
  RegionTimer Rsmoothed(Tsmoothed);

  UPtrSMdbl subst = BuildOffDiagSubstitutionMatrix(e2v, eweights, nv);
  cout << "Smoothed prol start" << endl;

  static Timer Tsmoothed_jacobi("H1 Create Smoothed Prol - Jacobi Mat");
  Tsmoothed_jacobi.Start();
  double avg = 0.5;
  ParallelFor (nv, [&] (int vert) {
    auto row_indices = subst->GetRowIndices(vert);
    auto row_vals = subst->GetRowValues(vert);
    double sum = 0.;
    for (auto v : row_vals) {
      sum += v;
    }
    for (auto i : row_indices) {
      // minus sign is in sum because we summed over negatives
      (*subst)(vert, i) *= avg/sum;
    }
    (*subst)(vert, vert) += 1. - avg;
  });
  Tsmoothed_jacobi.Stop();

  assert(subst->Width() == triv_prol->Height());
  UPtrSMdbl smoothed_prol =
    UPtrSMdbl(static_cast<SparseMatrixTM<double>*>(MatMult(*subst, *triv_prol)));

  return std::move(smoothed_prol);
}

UPtrSMdbl BuildOffDiagSubstitutionMatrix(
  const Array<INT<2>>& e2v, const Array<double>& eweights, int nv)
{
  static Timer Tsubst("H1 Build Subst Matrix");
  RegionTimer Rsubst(Tsubst);

  auto ne = e2v.Size();

  static Timer Tsubst_nze("H1 Build Subst Matrix - count nze");
  Tsubst_nze.Start();
  Array<int> nze(nv);
  ParallelFor(nv, [&] (int vert) {
    nze[vert] = 1;
  });

  ParallelFor(ne, [&] (int edge) {
    AsAtomic(nze[e2v[edge][0]])++;
    AsAtomic(nze[e2v[edge][1]])++;
  });

  // TODO: Check for 2d or 3d problem to decide maximum number of nze per row.
  // 1 entry for the corresponding coarse vertex plus 2d: 3, 3d: 4
  int max_nze_per_row = 5;  // use this for 3d. 4 + 1
  ParallelFor(nv, [&] (int vert) {
      nze[vert] = min(nze[vert], max_nze_per_row);
  });
  Tsubst_nze.Stop();

  static Timer Tsubst_mem("H1 Build Subst Matrix - memory");
  Tsubst_mem.Start();
  auto subst = std::make_unique<ngla::SparseMatrix<double>>(nze, nv);
  Tsubst_mem.Stop();

  static Timer Tsubst_table("H1 Build Subst Matrix - build table");
  Tsubst_table.Start();
  TableCreator<pair<int, double>> creator(nv);
  for (; !creator.Done(); creator++) {
    ParallelFor(nv, [&] (auto vert) {
        creator.Add(vert, make_pair<>(vert, 0));
    });
    ParallelFor(ne, [&] (auto edge) {
        auto ew = eweights[edge];
        auto v1 = e2v[edge][0];
        auto v2 = e2v[edge][1];

        creator.Add(v1, make_pair<>(v2, -ew));
        creator.Add(v2, make_pair<>(v1, -ew));
    });
  }
  auto table = creator.MoveTable();
  Tsubst_table.Stop();

  // row indices for sparse matrix need to be sorted
  static Timer Tsubst_sort("H1 Build Subst Matrix - sort rows in table");
  Tsubst_sort.Start();
  auto less_second = [] (pair<int, double> lhs, pair<int, double> rhs) {
    return lhs.second < rhs.second;
  };
  ParallelFor(table.Size(), [&] (auto row) {
    QuickSort(table[row], less_second);
  });
  Tsubst_sort.Stop();

  static Timer Tsubst_write("H1 Build Subst Matrix - write into matrix");
  Tsubst_write.Start();
  auto less_first = [] (pair<int, double> lhs, pair<int, double> rhs) {
    return lhs.first < rhs.first;
  };
  ParallelFor(table.Size(), [&] (auto row_i) {
    auto row = table[row_i];
    auto row_ind = subst->GetRowIndices(row_i);
    auto row_vals = subst->GetRowValues(row_i);
    ArrayMem<pair<int, double>, 5> new_row(row_ind.Size());
    new_row[0].first = row_i;
    new_row[0].second = 0;
    for (int i = 1; i < new_row.Size(); ++i) {
      new_row[i].first = row[i-1].first;
      new_row[i].second = row[i-1].second;
    }
    QuickSort(new_row, less_first);
    for (int i = 0; i < row_ind.Size(); ++i) {
      row_ind[i] = new_row[i].first;
      row_vals[i] = new_row[i].second;
    }
  });
  Tsubst_write.Stop();

  return std::move(subst);
}


}  // h1amg
