/*****************************************************************************
*
*     Program: EPTlib
*     Author: Alessandro Arduino <a.arduino@inrim.it>
*
*  MIT License
*
*  Copyright (c) 2020-2021  Alessandro Arduino
*  Istituto Nazionale di Ricerca Metrologica (INRiM)
*  Strada delle cacce 91, 10135 Torino
*  ITALY
*
*  Permission is hereby granted, free of charge, to any person obtaining a copy
*  of this software and associated documentation files (the "Software"), to deal
*  in the Software without restriction, including without limitation the rights
*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
*  copies of the Software, and to permit persons to whom the Software is
*  furnished to do so, subject to the following conditions:
*
*  The above copyright notice and this permission notice shall be included in all
*  copies or substantial portions of the Software.
*
*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
*  SOFTWARE.
*
*****************************************************************************/

#include "eptlib/finite_difference.h"

#include <complex>
#include <cstring>
#include <iostream>

#include "eptlib/linalg/linalg_util.h"
#include "eptlib/linalg/linalg_qr.h"
#include "eptlib/linalg/linalg_householder.h"

using namespace eptlib;

// FDSavitzkyGolayFilter constructor
FDSavitzkyGolayFilter::
FDSavitzkyGolayFilter(const Shape &shape) :
    shape_(shape),
    m_vox_(std::accumulate(shape.GetSize().begin(),shape.GetSize().end(),1,std::multiplies<int>())),
    n_unk_(0) {
    double toll = 1e-10;
    std::array<int,NDIM> mm = shape_.GetSize();
    // check to have odd shape size
    assert(m_vox_%2);
    // get the central voxel address
    std::array<int,NDIM> ii0;
    for (int d = 0; d<NDIM; ++d) {
        ii0[d] = mm[d]/2;
    }
    // initialise the design matrix
    int n_row = shape_.GetVolume();
    constexpr int n_col_start = 1 + NDIM + (NDIM*(NDIM+1))/2;
    int n_col = n_col_start;
    linalg::MatrixReal F(n_col,std::vector<double>(n_row));
    // initialise the derivative indexes
    std::array<DifferentialOperator,n_col_start> der_idx_inv;
    der_idx_inv[0] = DifferentialOperator::Field;
    der_idx_inv[1] = DifferentialOperator::GradientXX;
    der_idx_inv[2] = DifferentialOperator::GradientYY;
    der_idx_inv[3] = DifferentialOperator::GradientZZ;
    der_idx_inv[4] = DifferentialOperator::GradientX;
    der_idx_inv[5] = DifferentialOperator::GradientY;
    der_idx_inv[6] = DifferentialOperator::GradientZ;
    for (int c = 7; c<n_col_start; ++c) {
        der_idx_inv[c] = DifferentialOperator::END;
    }
    // fill the design matrix
    std::array<int,NDIM> ii;
    std::array<double,NDIM> di;
    int r = 0;
    for (ii[2] = 0; ii[2]<mm[2]; ++ii[2]) {
        for (ii[1] = 0; ii[1]<mm[1]; ++ii[1]) {
            for (ii[0] = 0; ii[0]<mm[0]; ++ii[0]) {
                if (shape_[ii]) {
                    for (int d = 0; d<NDIM; ++d) {
                        di[d] = ii[d]-ii0[d];
                    }
                    // design matrix for even quantities
                    F[0][r] = 1.0;
                    for (int d = 0; d<NDIM; ++d) {
                        F[1+d][r] = di[d]*di[d]/2.0;
                    }
                    // design matrix for odd quantities
                    int c = 0;
                    for (int d = 0; d<NDIM; ++d) {
                        F[1+NDIM + d][r] = di[d];
                        for (int d2 = d+1; d2<NDIM; ++d2) {
                            F[1+2*NDIM + c][r] = di[d]*di[d2];
                            ++c;
                        }
                    }
                    ++r;
                }
            }
        }
    }
    // check that all lines are filled
    std::vector<double> v(n_col);
    for (int c = 0; c<n_col; ++c) {
        v[c] = linalg::Norm2(F[c].data(),n_row);
    }
    double ref = linalg::MaxAbs(v.data(),n_col);
    for (int c = 0; c<n_col; ++c) {
        if (v[c]<toll*ref) {
            --n_col;
            F[c].swap(F[n_col]);
            std::swap(v[c],v[n_col]);
            std::swap(der_idx_inv[c],der_idx_inv[n_col]);
            --c;
        }
    }
    F.resize(n_col);
    n_unk_ = n_col;
    // identify derivative indexes in the design matrix
    for (int d = 0; d<der_idx_.size(); ++d) {
        der_idx_[d] = -1;
    }
    for (int c = 0; c<n_col; ++c) {
        if (der_idx_inv[c]!=DifferentialOperator::END) {
            der_idx_[static_cast<int>(der_idx_inv[c])] = c;
        }
    }
    // compute the QR factorisation
    linalg::HouseholderQR(&qr_,F,n_row,n_col);
    // invert the matrix
    linalg::MatrixReal A(n_row,std::vector<double>(n_col));
    for (int k = 0; k<n_row; ++k) {
        double *b = new double[n_row];
        std::memset(b,0,n_row*sizeof(double));
        b[k] = 1.0;
        linalg::QRSolve(A[k].data(),qr_,b,n_row,n_col);
        delete[] b;
    }
    // Compute the kernels
    for (int d = 0; d<NDIM; ++d) {
        int grad = static_cast<int>(DifferentialOperator::GradientX)+d;
        int lapl = static_cast<int>(DifferentialOperator::GradientXX)+d;
        if (der_idx_[grad]!=-1) {
            grad_kernel_[d].resize(n_row,0.0);
            for (int idx = 0; idx<n_row; ++idx) {
                grad_kernel_[d][idx] = A[idx][der_idx_[grad]];
            }
        }
        if (der_idx_[lapl]!=-1) {
            lapl_kernel_[d].resize(n_row,0.0);
            for (int idx = 0; idx<n_row; ++idx) {
                lapl_kernel_[d][idx] = A[idx][der_idx_[lapl]];
            }
        }
    }
    int field = static_cast<int>(DifferentialOperator::Field);
    if (der_idx_[field]) {
        field_kernel_.resize(n_row,0.0);
        for (int idx = 0; idx<n_row; ++idx) {
            field_kernel_[idx] = A[idx][der_idx_[field]];
        }
    }
    return;
}

// FDSavitzkyGolayFilter apply and compute the fitting quality
template <typename NumType>
EPTlibError FDSavitzkyGolayFilter::
Apply(const DifferentialOperator diff_op, NumType *dst, double *var, const NumType *src, const std::array<int,NDIM> &nn, const std::array<double,NDIM> &dd) const {
    if (diff_op!=DifferentialOperator::Laplacian) {
        if (der_idx_[static_cast<int>(diff_op)]==-1) {
            return EPTlibError::MissingData;
        }
    } else {
        bool can_compute_lapl = false;
        for (int d = 0; d<NDIM; ++d) {
            int lapl = static_cast<int>(DifferentialOperator::GradientXX)+d;
            if (der_idx_[lapl]!=-1) {
                can_compute_lapl = true;
            }
        }
        if (!can_compute_lapl) {
            return EPTlibError::MissingData;
        }
    }
    bool get_var = var!=nullptr;
    std::array<int,NDIM> ii;
    std::array<int,NDIM> rr;
    std::array<int,NDIM> inc;
    for (int d = 0; d<NDIM; ++d) {
        rr[d] = shape_.GetSize()[d]/2;
    }
    inc[0] = 1;
    inc[1] = nn[0]-shape_.GetSize()[0];
    inc[2] = nn[0]*(nn[1]-shape_.GetSize()[1]);
    // get the variance coefficient
    double var_coeff = 0.0;
    if (get_var) {
        std::vector<double> u(n_unk_,0.0);
        if (diff_op==DifferentialOperator::Field) {
            // zero order
            u[der_idx_[static_cast<int>(diff_op)]] = 1.0;
            var_coeff = EvaluateVariance(u);
        } else if (diff_op==DifferentialOperator::GradientX||
            diff_op==DifferentialOperator::GradientY||
            diff_op==DifferentialOperator::GradientZ) {
            // first order
            int d = static_cast<int>(diff_op) - static_cast<int>(DifferentialOperator::GradientX);
            u[der_idx_[static_cast<int>(diff_op)]] = 1.0/dd[d];
            var_coeff = EvaluateVariance(u);
        } else if (diff_op==DifferentialOperator::GradientXX||
            diff_op==DifferentialOperator::GradientYY||
            diff_op==DifferentialOperator::GradientZZ) {
            // second order
            int d = static_cast<int>(diff_op) - static_cast<int>(DifferentialOperator::GradientXX);
            u[der_idx_[static_cast<int>(diff_op)]] = 1.0/dd[d]/dd[d];
            var_coeff = EvaluateVariance(u);
        } else if (diff_op==DifferentialOperator::Laplacian) {
            // laplacian
            for (int d = 0; d<NDIM; ++d) {
                int lapl = static_cast<int>(DifferentialOperator::GradientXX)+d;
                if (der_idx_[lapl]!=-1) {
                    u[der_idx_[lapl]] = 1.0/dd[d]/dd[d];
                }
            }
            var_coeff = EvaluateVariance(u);
        }
    }
    // loop over field voxels
    for (ii[2] = rr[2]; ii[2]<nn[2]-rr[2]; ++ii[2]) {
        for (ii[1] = rr[1]; ii[1]<nn[1]-rr[1]; ++ii[1]) {
            for (ii[0] = rr[0]; ii[0]<nn[0]-rr[0]; ++ii[0]) {
                std::array<int,NDIM> ii_l;
                std::vector<NumType> field_crop(shape_.GetVolume(),0.0);
                // inner loop over kernel voxels
                int idx_f = 0;
                int idx_l = 0;
                int idx_g = ii[0]-rr[0] + nn[0]*(ii[1]-rr[1] + nn[1]*(ii[2]-rr[2]));
                for (ii_l[2] = -rr[2]; ii_l[2]<=rr[2]; ++ii_l[2]) {
                    for (ii_l[1] = -rr[1]; ii_l[1]<=rr[1]; ++ii_l[1]) {
                        for (ii_l[0] = -rr[0]; ii_l[0]<=rr[0]; ++ii_l[0]) {
                            if (shape_[idx_l]) {
                                // set the field_crop to the field value
                                field_crop[idx_f] = src[idx_g];
                                ++idx_f;
                            }
                            ++idx_l;
                            idx_g += inc[0];
                        }
                        idx_g += inc[1];
                    }
                    idx_g += inc[2];
                }
                // compute the derivative...
                int idx = ii[0] + nn[0]*(ii[1] + nn[1]*ii[2]);
                if (get_var) {
                    // ...getting the quality index
                    std::vector<NumType> a(n_unk_);
                    var[idx] = var_coeff*linalg::QRSolve(a.data(),qr_,field_crop.data(),field_crop.size(),n_unk_);
                    if (diff_op==DifferentialOperator::Field) {
                        // zero order
                        dst[idx] = a[der_idx_[static_cast<int>(diff_op)]];
                    } else if (diff_op==DifferentialOperator::GradientX||
                        diff_op==DifferentialOperator::GradientY||
                        diff_op==DifferentialOperator::GradientZ) {
                        // first order
                        int d = static_cast<int>(diff_op) - static_cast<int>(DifferentialOperator::GradientX);
                        dst[idx] = a[der_idx_[static_cast<int>(diff_op)]]/dd[d];
                    } else if (diff_op==DifferentialOperator::GradientXX||
                        diff_op==DifferentialOperator::GradientYY||
                        diff_op==DifferentialOperator::GradientZZ) {
                        // second order
                        int d = static_cast<int>(diff_op) - static_cast<int>(DifferentialOperator::GradientXX);
                        dst[idx] = a[der_idx_[static_cast<int>(diff_op)]]/dd[d]/dd[d];
                    } else if (diff_op==DifferentialOperator::Laplacian) {
                        // laplacian
                        dst[idx] = 0.0;
                        for (int d = 0; d<NDIM; ++d) {
                            int lapl = static_cast<int>(DifferentialOperator::GradientXX)+d;
                            if (der_idx_[lapl]!=-1) {
                                dst[idx] += a[der_idx_[lapl]]/dd[d]/dd[d];
                            }
                        }
                    }
                } else {
                    // ...without the quality index
                    if (diff_op==DifferentialOperator::Field) {
                        // zero order
                        dst[idx] = ZeroOrder(field_crop);
                    }
                    if (diff_op==DifferentialOperator::GradientX||
                        diff_op==DifferentialOperator::GradientY||
                        diff_op==DifferentialOperator::GradientZ) {
                        // first order
                        int d = static_cast<int>(diff_op) - static_cast<int>(DifferentialOperator::GradientX);
                        dst[idx] = FirstOrder(d,field_crop,dd);
                    } else if (diff_op==DifferentialOperator::GradientXX||
                        diff_op==DifferentialOperator::GradientYY||
                        diff_op==DifferentialOperator::GradientZZ) {
                        // second order
                        int d = static_cast<int>(diff_op) - static_cast<int>(DifferentialOperator::GradientXX);
                        dst[idx] = SecondOrder(d,field_crop,dd);
                    } else if (diff_op==DifferentialOperator::Laplacian) {
                        // laplacian
                        dst[idx] = Laplacian(field_crop,dd);
                    }
                }
            }
        }
    }
    return EPTlibError::Success;
}
// FDSavitzkyGolayFilter apply
template <typename NumType>
EPTlibError FDSavitzkyGolayFilter::
Apply(const DifferentialOperator diff_op, NumType *dst, const NumType *src, const std::array<int,NDIM> &nn, const std::array<double,NDIM> &dd) const {
    return Apply(diff_op,dst,nullptr,src,nn,dd);
}

// FDSavitzkyGolayFilter apply
EPTlibError FDSavitzkyGolayFilter::
ApplyWrappedPhase(const DifferentialOperator diff_op, double *dst, const double *src, const std::array<int,NDIM> &nn, const std::array<double,NDIM> &dd) const {
    if (diff_op!=DifferentialOperator::Laplacian) {
        if (der_idx_[static_cast<int>(diff_op)]==-1) {
            return EPTlibError::MissingData;
        }
    } else {
        bool can_compute_lapl = false;
        for (int d = 0; d<NDIM; ++d) {
            int lapl = static_cast<int>(DifferentialOperator::GradientXX)+d;
            if (der_idx_[lapl]!=-1) {
                can_compute_lapl = true;
            }
        }
        if (!can_compute_lapl) {
            return EPTlibError::MissingData;
        }
    }
    std::array<int,NDIM> ii;
    std::array<int,NDIM> rr;
    std::array<int,NDIM> inc;
    for (int d = 0; d<NDIM; ++d) {
        rr[d] = shape_.GetSize()[d]/2;
    }
    inc[0] = 1;
    inc[1] = nn[0]-shape_.GetSize()[0];
    inc[2] = nn[0]*(nn[1]-shape_.GetSize()[1]);
    // loop over field voxels
    for (ii[2] = rr[2]; ii[2]<nn[2]-rr[2]; ++ii[2]) {
        for (ii[1] = rr[1]; ii[1]<nn[1]-rr[1]; ++ii[1]) {
            for (ii[0] = rr[0]; ii[0]<nn[0]-rr[0]; ++ii[0]) {
                std::array<int,NDIM> ii_l;
                std::vector<std::complex<double> > field_crop(shape_.GetVolume(),0.0);
                // inner loop over kernel voxels
                int idx_f = 0;
                int idx_l = 0;
                int idx_g = ii[0]-rr[0] + nn[0]*(ii[1]-rr[1] + nn[1]*(ii[2]-rr[2]));
                for (ii_l[2] = -rr[2]; ii_l[2]<=rr[2]; ++ii_l[2]) {
                    for (ii_l[1] = -rr[1]; ii_l[1]<=rr[1]; ++ii_l[1]) {
                        for (ii_l[0] = -rr[0]; ii_l[0]<=rr[0]; ++ii_l[0]) {
                            if (shape_[idx_l]) {
                                // set the field_crop to the field value
                                field_crop[idx_f] = std::exp(std::complex<double>(0.0,src[idx_g]));
                                ++idx_f;
                            }
                            ++idx_l;
                            idx_g += inc[0];
                        }
                        idx_g += inc[1];
                    }
                    idx_g += inc[2];
                }
                // compute the derivative...
                int idx = ii[0] + nn[0]*(ii[1] + nn[1]*ii[2]);
                // ...without the quality index
                if (diff_op==DifferentialOperator::Field) {
                    // zero order
                    dst[idx] = std::log(ZeroOrder(field_crop)).imag();
                }
                if (diff_op==DifferentialOperator::GradientX||
                    diff_op==DifferentialOperator::GradientY||
                    diff_op==DifferentialOperator::GradientZ) {
                    // first order
                    int d = static_cast<int>(diff_op) - static_cast<int>(DifferentialOperator::GradientX);
                    std::complex<double> tmp = ZeroOrder(field_crop);
                    std::complex<double> der = FirstOrder(d,field_crop,dd);
                    dst[idx] = (der/tmp).imag();
                } else if (diff_op==DifferentialOperator::GradientXX||
                    diff_op==DifferentialOperator::GradientYY||
                    diff_op==DifferentialOperator::GradientZZ) {
                    // second order
                    int d = static_cast<int>(diff_op) - static_cast<int>(DifferentialOperator::GradientXX);
                    std::complex<double> tmp = ZeroOrder(field_crop);
                    std::complex<double> der = FirstOrder(d,field_crop,dd);
                    std::complex<double> dder = SecondOrder(d,field_crop,dd);
                    dst[idx] = (der*der/tmp/tmp + dder/tmp).imag();
                } else if (diff_op==DifferentialOperator::Laplacian) {
                    // laplacian
                    dst[idx] = 0.0;
                    for (int d = 0; d<NDIM; ++d) {
                        std::complex<double> tmp = ZeroOrder(field_crop);
                        std::complex<double> der = FirstOrder(d,field_crop,dd);
                        std::complex<double> dder = SecondOrder(d,field_crop,dd);
                        dst[idx] += (der*der/tmp/tmp + dder/tmp).imag();
                    }
                }
            }
        }
    }
    return EPTlibError::Success;
}

// FDSavitzkyGolayFilter apply kernel zero order derivative
template <typename NumType>
NumType FDSavitzkyGolayFilter::
ZeroOrder(const std::vector<NumType> &field_crop) const {
    NumType dst = std::inner_product(field_kernel_.begin(),field_kernel_.end(),field_crop.begin(),static_cast<NumType>(0.0));
    return dst;
}
// FDSavitzkyGolayFilter apply kernel first order derivative
template <typename NumType>
NumType FDSavitzkyGolayFilter::
FirstOrder(const int d, const std::vector<NumType> &field_crop, const std::array<double,NDIM> &dd) const {
    NumType dst = std::inner_product(grad_kernel_[d].begin(),grad_kernel_[d].end(),field_crop.begin(),static_cast<NumType>(0.0))/dd[d];
    return dst;
}
// FDSavitzkyGolayFilter apply kernel second order derivative
template <typename NumType>
NumType FDSavitzkyGolayFilter::
SecondOrder(const int d, const std::vector<NumType> &field_crop, const std::array<double,NDIM> &dd) const {
    NumType dst = std::inner_product(lapl_kernel_[d].begin(),lapl_kernel_[d].end(),field_crop.begin(),static_cast<NumType>(0.0))/dd[d]/dd[d];
    return dst;
}
// FDSavitzkyGolayFilter apply kernel laplacian
template <typename NumType>
NumType FDSavitzkyGolayFilter::
Laplacian(const std::vector<NumType> &field_crop, const std::array<double,NDIM> &dd) const {
    NumType dst = 0.0;
    for (int d = 0; d<NDIM; ++d) {
        dst += std::inner_product(lapl_kernel_[d].begin(),lapl_kernel_[d].end(),field_crop.begin(),static_cast<NumType>(0.0))/dd[d]/dd[d];
    }
    return dst;
}

// FDSavitzkyGolayFilter getters
const Shape& FDSavitzkyGolayFilter::
GetShape() const {
    return shape_;
}

// FDSavitzkyGolayFilter evaluate the uncertainty
double FDSavitzkyGolayFilter::
EvaluateVariance(const std::vector<double> &u) const {
    assert(u.size()==n_unk_);
    std::vector<double> r(n_unk_,0.0);
    // solve the equation with the transpose of R
    r[0] = u[0]/qr_[0][0];
    for (int i = 1; i<n_unk_; ++i) {
        double s = 0;
        for (int j = 0; j<i; ++j) {
            s += r[j]*qr_[i][j];
        }
        r[i] = (u[i]-s)/qr_[i][i];
    }
    // compute the variance
    double sigma = linalg::Dot(r.data(),r.data(),n_unk_);
    return sigma;
}

// FDSavitzkyGolayFilter specialisations
template EPTlibError FDSavitzkyGolayFilter::Apply<double>(const DifferentialOperator diff_op,double *dst,double *var, const double *src, const std::array<int,NDIM> &nn, const std::array<double,NDIM> &dd) const;
template EPTlibError FDSavitzkyGolayFilter::Apply<std::complex<double> >(const DifferentialOperator diff_op,std::complex<double> *dst,double *var, const std::complex<double> *src, const std::array<int,NDIM> &nn, const std::array<double,NDIM> &dd) const;

template EPTlibError FDSavitzkyGolayFilter::Apply<double>(const DifferentialOperator diff_op,double *dst, const double *src, const std::array<int,NDIM> &nn, const std::array<double,NDIM> &dd) const;
template EPTlibError FDSavitzkyGolayFilter::Apply<std::complex<double> >(const DifferentialOperator diff_op,std::complex<double> *dst, const std::complex<double> *src, const std::array<int,NDIM> &nn, const std::array<double,NDIM> &dd) const;

template double FDSavitzkyGolayFilter::FirstOrder<double>(const int d,const std::vector<double> &field_crop,const std::array<double,NDIM> &dd) const;
template std::complex<double> FDSavitzkyGolayFilter::FirstOrder<std::complex<double> >(const int d,const std::vector<std::complex<double> > &field_crop,const std::array<double,NDIM> &dd) const;

template double FDSavitzkyGolayFilter::SecondOrder<double>(const int d,const std::vector<double> &field_crop,const std::array<double,NDIM> &dd) const;
template std::complex<double> FDSavitzkyGolayFilter::SecondOrder<std::complex<double> >(const int d,const std::vector<std::complex<double> > &field_crop,const std::array<double,NDIM> &dd) const;

template double FDSavitzkyGolayFilter::Laplacian<double>(const std::vector<double> &field_crop,const std::array<double,NDIM> &dd) const;
template std::complex<double> FDSavitzkyGolayFilter::Laplacian<std::complex<double> >(const std::vector<std::complex<double> > &field_crop,const std::array<double,NDIM> &dd) const;
