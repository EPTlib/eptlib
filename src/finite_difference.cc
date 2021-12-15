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
    int idx0 = m_vox_/2;
    std::array<int,NDIM> ii0;
    for (int d = 0; d<NDIM; ++d) {
        ii0[d] = mm[d]/2;
    }
    // initialise the design matrix
    int n_row = shape_.GetVolume();
    int n_col = 1 + NDIM + (NDIM*(NDIM+1))/2;
    linalg::MatrixReal F(n_col,std::vector<real_t>(n_row));
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
                        F[1+d][r] = di[d]*di[d];
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
    std::vector<real_t> v(n_col);
    for (int c = 0; c<n_col; ++c) {
        v[c] = linalg::Norm2(F[c].data(),n_row);
    }
    real_t ref = linalg::MaxAbs(v.data(),n_col);
    for (int d = 0; d<1+2*NDIM; ++d) {
        if (v[d]<toll*ref) {
            throw std::runtime_error("Impossible to set-up the Savitzky-Golay filter with the provided shape: no points along direction "+std::to_string(d)+".");
        }
    }
    for (int c = 1+2*NDIM; c<n_col; ++c) {
        if (v[c]<toll*ref) {
            --n_col;
            F[c].swap(F[n_col]);
            --c;
        }
    }
    F.resize(n_col);
    n_unk_ = n_col;
    // compute the QR factorisation
    linalg::HouseholderQR(&qr_,F,n_row,n_col);
    // invert the matrix
    linalg::MatrixReal A(n_row,std::vector<real_t>(n_col));
    for (int k = 0; k<n_row; ++k) {
        real_t *b = new real_t[n_row];
        std::memset(b,0,n_row*sizeof(real_t));
        b[k] = 1.0;
        linalg::QRSolve(A[k].data(),qr_,b,n_row,n_col);
        delete[] b;
    }
    // Compute the kernel for the laplacian
    for (int d = 0; d<NDIM; ++d) {
        lapl_kernel_[d].resize(m_vox_,0.0);
        r = 0;
        for (int idx = 0; idx<m_vox_; ++idx) {
            if (shape_[idx]) {
                lapl_kernel_[d][idx] = 2.0*A[r][1+d];
                ++r;
            }
        }
    }
    // Compute the kernel for the gradient
    for (int d = 0; d<NDIM; ++d) {
        grad_kernel_[d].resize(m_vox_,0.0);
        r = 0;
        for (int idx = 0; idx<m_vox_; ++idx) {
            if (shape_[idx]) {
                grad_kernel_[d][idx] = A[r][1+NDIM+d];
                ++r;
            }
        }
    }
    return;
}

// FDSavitzkyGolayFilter apply and compute the fitting quality
template <typename NumType>
EPTlibError FDSavitzkyGolayFilter::
ApplyChi2(const DifferentialOperator diff_op, NumType *dst, real_t *chi2, const NumType *src, const std::array<int,NDIM> &nn, const std::array<double,NDIM> &dd) const {
    const int n_vox = std::accumulate(nn.begin(),nn.end(),1,std::multiplies<int>());
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
                // compute the derivative
                std::vector<NumType> a(n_unk_);
                int idx = ii[0] + nn[0]*(ii[1] + nn[1]*ii[2]);
                chi2[idx] = linalg::QRSolve(a.data(),qr_,field_crop.data(),field_crop.size(),n_unk_);
                if (diff_op==DifferentialOperator::GradientX||
                    diff_op==DifferentialOperator::GradientY||
                    diff_op==DifferentialOperator::GradientZ) {
                    int d = static_cast<int>(diff_op) -
                        static_cast<int>(DifferentialOperator::GradientX);
                    dst[idx] = a[1+NDIM+d]/dd[d];
                } else if (diff_op==DifferentialOperator::GradientXX||
                    diff_op==DifferentialOperator::GradientYY||
                    diff_op==DifferentialOperator::GradientZZ) {
                    int d = static_cast<int>(diff_op) -
                        static_cast<int>(DifferentialOperator::GradientXX);
                    dst[idx] = 2.0*a[1+d]/dd[d]/dd[d];
                } else if (diff_op==DifferentialOperator::Laplacian) {
                    dst[idx] = 0.0;
                    for (int d = 0; d<NDIM; ++d) {
                        dst[idx] += 2.0*a[1+d]/dd[d]/dd[d];
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
    const int n_vox = std::accumulate(nn.begin(),nn.end(),1,std::multiplies<int>());
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
                std::vector<NumType> field_crop(m_vox_,0.0);
                // inner loop over kernel voxels
                int idx_l = 0;
                int idx_g = ii[0]-rr[0] + nn[0]*(ii[1]-rr[1] + nn[1]*(ii[2]-rr[2]));
                for (ii_l[2] = -rr[2]; ii_l[2]<=rr[2]; ++ii_l[2]) {
                    for (ii_l[1] = -rr[1]; ii_l[1]<=rr[1]; ++ii_l[1]) {
                        for (ii_l[0] = -rr[0]; ii_l[0]<=rr[0]; ++ii_l[0]) {
                            if (shape_[idx_l]) {
                                // set the field_crop to the field value
                                field_crop[idx_l] = src[idx_g];
                            }
                            ++idx_l;
                            idx_g += inc[0];
                        }
                        idx_g += inc[1];
                    }
                    idx_g += inc[2];
                }
                // compute the derivative
                int idx = ii[0] + nn[0]*(ii[1] + nn[1]*ii[2]);
                if (diff_op==DifferentialOperator::GradientX||
                    diff_op==DifferentialOperator::GradientY||
                    diff_op==DifferentialOperator::GradientZ) {
                    int d = static_cast<int>(diff_op) -
                        static_cast<int>(DifferentialOperator::GradientX);
                    dst[idx] = FirstOrder(d,field_crop,dd);
                } else if (diff_op==DifferentialOperator::GradientXX||
                    diff_op==DifferentialOperator::GradientYY||
                    diff_op==DifferentialOperator::GradientZZ) {
                    int d = static_cast<int>(diff_op) -
                        static_cast<int>(DifferentialOperator::GradientXX);
                    dst[idx] = SecondOrder(d,field_crop,dd);
                } else if (diff_op==DifferentialOperator::Laplacian) {
                    dst[idx] = Laplacian(field_crop,dd);
                }
            }
        }
    }
    return EPTlibError::Success;
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
const std::array<std::vector<double>,NDIM>& FDSavitzkyGolayFilter::
GetLaplKernel() const {
    return lapl_kernel_;
}
const std::array<std::vector<double>,NDIM>& FDSavitzkyGolayFilter::
GetGradKernel() const {
    return grad_kernel_;
}

// FDSavitzkyGolayFilter specialisations
template EPTlibError FDSavitzkyGolayFilter::ApplyChi2<double>(const DifferentialOperator diff_op,double *dst,double *chi2, const double *src, const std::array<int,NDIM> &nn, const std::array<double,NDIM> &dd) const;
template EPTlibError FDSavitzkyGolayFilter::ApplyChi2<std::complex<double> >(const DifferentialOperator diff_op,std::complex<double> *dst,double *chi2, const std::complex<double> *src, const std::array<int,NDIM> &nn, const std::array<double,NDIM> &dd) const;

template EPTlibError FDSavitzkyGolayFilter::Apply<double>(const DifferentialOperator diff_op,double *dst, const double *src, const std::array<int,NDIM> &nn, const std::array<double,NDIM> &dd) const;
template EPTlibError FDSavitzkyGolayFilter::Apply<std::complex<double> >(const DifferentialOperator diff_op,std::complex<double> *dst, const std::complex<double> *src, const std::array<int,NDIM> &nn, const std::array<double,NDIM> &dd) const;

template double FDSavitzkyGolayFilter::FirstOrder<double>(const int d,const std::vector<double> &field_crop,const std::array<double,NDIM> &dd) const;
template std::complex<double> FDSavitzkyGolayFilter::FirstOrder<std::complex<double> >(const int d,const std::vector<std::complex<double> > &field_crop,const std::array<double,NDIM> &dd) const;

template double FDSavitzkyGolayFilter::SecondOrder<double>(const int d,const std::vector<double> &field_crop,const std::array<double,NDIM> &dd) const;
template std::complex<double> FDSavitzkyGolayFilter::SecondOrder<std::complex<double> >(const int d,const std::vector<std::complex<double> > &field_crop,const std::array<double,NDIM> &dd) const;

template double FDSavitzkyGolayFilter::Laplacian<double>(const std::vector<double> &field_crop,const std::array<double,NDIM> &dd) const;
template std::complex<double> FDSavitzkyGolayFilter::Laplacian<std::complex<double> >(const std::vector<std::complex<double> > &field_crop,const std::array<double,NDIM> &dd) const;

// Apply the FD filter to a wrapped phase input field.
EPTlibError eptlib::WrappedPhaseDerivative(const DifferentialOperator diff_op,
    double *dst, const double *src, const std::array<int,NDIM> &nn,
    const std::array<double,NDIM> &dd, const FDSavitzkyGolayFilter &fd_filter) {
    int n_vox = std::accumulate(nn.begin(),nn.end(),1,std::multiplies<int>());
    // initialize the exponential map
    std::vector<std::complex<double> > emap(n_vox);
    for (int idx = 0; idx<n_vox; ++idx) {
        emap[idx] = std::exp(std::complex<double>(0,src[idx]));
    }
    // compute the derivative...
    if (diff_op==DifferentialOperator::GradientX||
        diff_op==DifferentialOperator::GradientY||
        diff_op==DifferentialOperator::GradientZ) {
        // ...first order
        std::vector<std::complex<double> > d_emap(n_vox);
        fd_filter.Apply(diff_op,d_emap.data(),emap.data(),nn,dd);
        for (int idx = 0; idx<n_vox; ++idx) {
            dst[idx] = std::imag(d_emap[idx]/emap[idx]);
        }
    } else if (diff_op==DifferentialOperator::GradientXX||
        diff_op==DifferentialOperator::GradientYY||
        diff_op==DifferentialOperator::GradientZZ) {
        // ...second order
        DifferentialOperator diff_op1 = static_cast<DifferentialOperator>(static_cast<int>(diff_op)-static_cast<int>(DifferentialOperator::GradientXX));
        std::vector<std::complex<double> > d_emap(n_vox);
        std::vector<std::complex<double> > dd_emap(n_vox);
        fd_filter.Apply(diff_op1,d_emap.data(),emap.data(),nn,dd);
        fd_filter.Apply(diff_op,dd_emap.data(),emap.data(),nn,dd);
        for (int idx = 0; idx<n_vox; ++idx) {
            dst[idx] = std::imag(d_emap[idx]*d_emap[idx]/emap[idx]/emap[idx] + dd_emap[idx]/emap[idx]);
        }
    } else if (diff_op==DifferentialOperator::Laplacian) {
        // ...laplacian
        std::vector<std::complex<double> > dx_emap(n_vox);
        std::vector<std::complex<double> > dy_emap(n_vox);
        std::vector<std::complex<double> > dz_emap(n_vox);
        std::vector<std::complex<double> > dd_emap(n_vox);
        fd_filter.Apply(DifferentialOperator::GradientX,dx_emap.data(),emap.data(),nn,dd);
        fd_filter.Apply(DifferentialOperator::GradientY,dy_emap.data(),emap.data(),nn,dd);
        fd_filter.Apply(DifferentialOperator::GradientZ,dz_emap.data(),emap.data(),nn,dd);
        fd_filter.Apply(DifferentialOperator::Laplacian,dd_emap.data(),emap.data(),nn,dd);
        for (int idx = 0; idx<n_vox; ++idx) {
            dst[idx] = std::imag((dx_emap[idx]*dx_emap[idx]+dy_emap[idx]*dy_emap[idx]+dz_emap[idx]*dz_emap[idx])/emap[idx]/emap[idx] + dd_emap[idx]/emap[idx]);
        }
    }
    return EPTlibError::Success;
}

// Apply the FD filter to a wrapped phase input field and compute the fitting quality.
EPTlibError eptlib::WrappedPhaseDerivativeChi2(const DifferentialOperator diff_op,
    double *dst, real_t *chi2, const double *src, const std::array<int,NDIM> &nn,
    const std::array<double,NDIM> &dd, const FDSavitzkyGolayFilter &fd_filter) {
    int n_vox = std::accumulate(nn.begin(),nn.end(),1,std::multiplies<int>());
    // initialize the exponential map
    std::vector<std::complex<double> > emap(n_vox);
    for (int idx = 0; idx<n_vox; ++idx) {
        emap[idx] = std::exp(std::complex<double>(0,src[idx]));
    }
    // compute the derivative...
    if (diff_op==DifferentialOperator::GradientX||
        diff_op==DifferentialOperator::GradientY||
        diff_op==DifferentialOperator::GradientZ) {
        // ...first order
        std::vector<std::complex<double> > d_emap(n_vox);
        fd_filter.ApplyChi2(diff_op,d_emap.data(),chi2,emap.data(),nn,dd);
        for (int idx = 0; idx<n_vox; ++idx) {
            dst[idx] = std::imag(d_emap[idx]/emap[idx]);
        }
    } else if (diff_op==DifferentialOperator::GradientXX||
        diff_op==DifferentialOperator::GradientYY||
        diff_op==DifferentialOperator::GradientZZ) {
        // ...second order
        DifferentialOperator diff_op1 = static_cast<DifferentialOperator>(static_cast<int>(diff_op)-static_cast<int>(DifferentialOperator::GradientXX));
        std::vector<std::complex<double> > d_emap(n_vox);
        std::vector<std::complex<double> > dd_emap(n_vox);
        fd_filter.Apply(diff_op1,d_emap.data(),emap.data(),nn,dd);
        fd_filter.ApplyChi2(diff_op,dd_emap.data(),chi2,emap.data(),nn,dd);
        for (int idx = 0; idx<n_vox; ++idx) {
            dst[idx] = std::imag(d_emap[idx]*d_emap[idx]/emap[idx]/emap[idx] + dd_emap[idx]/emap[idx]);
        }
    } else if (diff_op==DifferentialOperator::Laplacian) {
        // ...laplacian
        std::vector<std::complex<double> > dx_emap(n_vox);
        std::vector<std::complex<double> > dy_emap(n_vox);
        std::vector<std::complex<double> > dz_emap(n_vox);
        std::vector<std::complex<double> > dd_emap(n_vox);
        fd_filter.Apply(DifferentialOperator::GradientX,dx_emap.data(),emap.data(),nn,dd);
        fd_filter.Apply(DifferentialOperator::GradientY,dy_emap.data(),emap.data(),nn,dd);
        fd_filter.Apply(DifferentialOperator::GradientZ,dz_emap.data(),emap.data(),nn,dd);
        fd_filter.ApplyChi2(DifferentialOperator::Laplacian,dd_emap.data(),chi2,emap.data(),nn,dd);
        for (int idx = 0; idx<n_vox; ++idx) {
            dst[idx] = std::imag((dx_emap[idx]*dx_emap[idx]+dy_emap[idx]*dy_emap[idx]+dz_emap[idx]*dz_emap[idx])/emap[idx]/emap[idx] + dd_emap[idx]/emap[idx]);
        }
    }
    return EPTlibError::Success;
}
