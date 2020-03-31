/*****************************************************************************
*
*     Program: EPTlib
*     Author: Alessandro Arduino <a.arduino@inrim.it>
*
*  MIT License
*
*  Copyright (c) 2020  Alessandro Arduino
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

#include <Eigen/Dense>

using namespace eptlib;

// FDLaplacianKernel constructor
FDLaplacianKernel::
FDLaplacianKernel(const Shape &shape) :
    shape_(shape),
    m_vox_(std::accumulate(shape.GetSize().begin(),shape.GetSize().end(),1,std::multiplies<int>())) {
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
    int n_col = 1 + NDIM;
    if (!shape_.IsSymmetric()) {
        n_col += (NDIM*(NDIM+1))/2;
    }
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> F(n_row,n_col);
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
                    F(r,0) = 1.0;
                    for (int d = 0; d<NDIM; ++d) {
                        F(r,1+d) = di[d]*di[d];
                    }
                    if (!shape_.IsSymmetric()) {
                        int c = 0;
                        for (int d = 0; d<NDIM; ++d) {
                            F(r,1+NDIM+d) = di[d];
                            for (int d2 = d+1; d2<NDIM; ++d2) {
                                F(r,1+2*NDIM+c) = di[d]*di[d2];
                                ++c;
                            }
                        }
                    }
                    ++r;
                }
            }
        }
    }
    // solve the normal equations
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> A = (F.transpose()*F).inverse()*F.transpose();
    // Compute the kernel
    for (int d = 0; d<NDIM; ++d) {
        kernel_[d].resize(m_vox_);
        r = 0;
        for (int idx = 0; idx<m_vox_; ++idx) {
            if (shape_[idx]) {
                kernel_[d][idx] = 2.0*A(1+d,r);
                ++r;
            }
        }
    }
    return;
}

// FDLaplacianKernel apply
template <typename NumType>
EPTlibError_t FDLaplacianKernel::
ComputeLaplacian(NumType *dst, const NumType *src, const std::array<int,NDIM> &nn,
    const std::array<double,NDIM> &dd) {
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
                // compute the laplacian
                int idx = ii[0] + nn[0]*(ii[1] + nn[1]*ii[2]);
                dst[idx] = 0.0;
                for (int d = 0; d<NDIM; ++d) {
                    dst[idx] += std::inner_product(kernel_[d].begin(),kernel_[d].end(),field_crop.begin(),static_cast<NumType>(0.0))/dd[d]/dd[d];
                }
            }
        }
    }
    return EPTlibError::Success;
}

// FDLaplacianKernel specialisations
template EPTlibError_t FDLaplacianKernel::ComputeLaplacian<double>(double *dst, const double *src, const std::array<int,NDIM> &nn, const std::array<double,NDIM> &dd);
template EPTlibError_t FDLaplacianKernel::ComputeLaplacian<std::complex<double>>(std::complex<double> *dst, const std::complex<double> *src, const std::array<int,NDIM> &nn, const std::array<double,NDIM> &dd);
