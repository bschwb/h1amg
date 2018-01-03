#ifndef H1AMG_H1_SMOOTHED_PROL_HPP_
#define H1AMG_H1_SMOOTHED_PROL_HPP_

#include <la.hpp>
#include <ngstd.hpp>

namespace h1amg
{

using UPtrSMdbl = std::unique_ptr<ngla::SparseMatrixTM<double>>;
using SPtrSMdbl = std::shared_ptr<ngla::SparseMatrixTM<double>>;

SPtrSMdbl H1SmoothedProl(
    const ngstd::Array<int>& vertex_coarse, int ncv, const ngstd::Array<ngstd::INT<2>>& e2v,
    const ngstd::Array<double>& ew, bool complx=false);

SPtrSMdbl CreateSmoothedProlongation(
    const ngstd::Array<ngstd::INT<2>>& e2v, const ngstd::Array<double>& eweights, int nv,
    const UPtrSMdbl triv_prol);

}  // h1amg

#endif  // H1AMG_H1_SMOOTHED_PROL_HPP_
