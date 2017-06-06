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
  for (const auto k : Range(nv)) {
    if (vertex_coarse[k] != -1) {
      c2f[vertex_coarse[k]][(c2f[vertex_coarse[k]][0]==-1)?0:1] = k;
    }
  }

  // sort entries
  for (auto r : c2f) {
    if (r[0]>r[1]) {
      auto d = r[0];
      r[0] = r[1];
      r[1] = d;
    }
  }

  return c2f;
}

// its a square matrix
// econ(i, j) = #e for which e=(vert_i, vert_j)
// so rows and cols correspond to vertices and the entry is the edge number
SparseMatrix<double> EdgeConnectivityMatrix(const Array<INT<2>>& e2v, int nv)
{
  auto ne = e2v.Size();

  // count nr of connected edges of vertex
  Array<int> econ_s(nv);
  econ_s = 0;

  for (auto k : Range(ne)) {
    econ_s[e2v[k][0]]++;
    econ_s[e2v[k][1]]++;
  }

  // build the table for the connections
  Table<int> tab_econ(econ_s);
  econ_s = 0; // used again for counting!!!!
  for (auto k : Range(ne)) {
    tab_econ[e2v[k][0]][econ_s[e2v[k][0]]++] = k;
    tab_econ[e2v[k][1]][econ_s[e2v[k][1]]++] = k;
  }

  // build the edge connection matrix
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

  // Sanity checks
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

  return econ;
}


UPtrSMdbl H1SemiSmoothedProl(
    const Array<int>& vertex_coarse, int ncv, const Array<INT<2>>& e2v, const Array<double>& ew,
    bool complx)
{
  static Timer tt ("fictitious prolongation (single) - total");

  int nverts = vertex_coarse.Size();
  int nedges = e2v.Size();

  Array<int> twos(ncv);
  twos = 2;
  Table<int> c2f(twos);
  for (auto r:c2f) {
    r = -1;
  }
  for (auto k:Range(nverts)) {
    if (vertex_coarse[k]!=-1) {
      c2f[vertex_coarse[k]][(c2f[vertex_coarse[k]][0]==-1)?0:1] = k;
    }
  }

  auto econ = EdgeConnectivityMatrix(e2v, nverts);
  const SparseMatrix<double> & cecon(econ);

  TableCreator<int> create_graph(nverts);
  while (!create_graph.Done()) {
    for (auto dof_f : Range(nverts)) {
      if (vertex_coarse[dof_f] != -1) {
        if (c2f[vertex_coarse[dof_f]][1] == -1) {
          create_graph.Add(dof_f, vertex_coarse[dof_f]);
        }
        else {
          auto row = econ.GetRowIndices(dof_f);
          Array<int> other_c;
          for (auto d:row) {
            if (vertex_coarse[d]!=-1 && !other_c.Contains(vertex_coarse[d])) {
              other_c.Append(vertex_coarse[d]);
            }
          }
          QuickSort(other_c);
          for (auto d:other_c) {
            create_graph.Add(dof_f, d);
          }
        }
      }
    }
    create_graph++;
  }

  auto graph = create_graph.MoveTable();
  Array<int> nze_perow(nverts);
  for (auto k:Range(nverts)) {
    nze_perow[k] = graph[k].Size();
  }

  UPtrSMdbl prol = nullptr;
  if (!complx) {
    prol = UPtrSMdbl(new SparseMatrix<double>(nze_perow, ncv));
  } else {
    prol = UPtrSMdbl(new SparseMatrix<double, Complex, Complex>(nze_perow, ncv));
  }

  for (auto rnr : Range(nverts)) {
    auto d = prol->GetRowIndices(rnr);
    auto gr = graph[rnr];
    for (auto k : Range(gr.Size())) {
      d[k] = gr[k];
    }
  }

  Array<int> c2lc(ncv);
  c2lc = -1;
  for (auto dof_c : Range(ncv)) {
    if (c2f[dof_c][1]==-1) {
      prol->GetRowValues(c2f[dof_c][0])[0] = 1.0;
    }
    else {
      for (auto l : Range(2)) {
        auto dof_f = c2f[dof_c][l];
        auto other_dof = c2f[dof_c][1-l];
        double wt_other_dof = ew[(int)cecon(dof_f, other_dof)];
        Array<int> other_ds;
        Array<double> wts;
        double s = wt_other_dof;
        for (auto j : cecon.GetRowIndices(dof_f)) {
          if (j != other_dof && vertex_coarse[j] != -1) {
            other_ds.Append(j);
            auto ed_wt = ew[(int)cecon(dof_f, j)];
            wts.Append(ed_wt);
            s += ed_wt;
          }
        }
        for (auto k : Range(graph[dof_f].Size())) {
          c2lc[graph[dof_f][k]] = k;
        }

        auto sinv = 1.0/s;
        auto vals = prol->GetRowValues(dof_f);
        vals = 0.0;
        int indl = c2lc[dof_c];
        vals[indl] += wt_other_dof;

        for (auto k : Range(other_ds.Size())) {
          vals[indl] += 0.5*wts[k];
          vals[c2lc[vertex_coarse[other_ds[k]]]] += 0.5*wts[k];
        }
        for (auto k : Range(vals.Size())) {
          vals[k] *= sinv;
        }
      }
    }
  }

  return std::move(prol);
}


}  // h1amg
