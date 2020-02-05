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

#include <Eigen/Dense>

using namespace eptlib;

// FDLaplacianKernel constructor
FDLaplacianKernel::
FDLaplacianKernel(const Shape &shape) :
    n_dim_(shape.GetNDim()), n_vox_(shape.GetSize()), shape_(shape) {
    // check to have odd shape size
    assert(shape_.GetSize()%2);
    // copy the content of nn and dd
    nn_.resize(n_dim_);
    for (int d = 0; d<n_dim_; ++d) {
        nn_[d] = shape_.GetSize(d);
    }
    // get the central voxel address
    int idx0 = shape_.GetSize()/2;
    std::vector<int> ii0(n_dim_);
    IdxToMultiIdx(ii0,idx0,nn_);
    // initialise the design matrix
    int n_row = shape_.GetVolume();
    int n_col = 1 + n_dim_;
    if (!shape_.IsSymmetric()) {
        n_col += (n_dim_*(n_dim_+1))/2;
    }
    Eigen::Matrix<real_t,Eigen::Dynamic,Eigen::Dynamic> F(n_row,n_col);
    // fill the design matrix
    std::vector<int> ii(n_dim_);
    std::vector<real_t> di(n_dim_);
    int r = 0;
    for (int idx = 0; idx<n_vox_; ++idx) {
        if (shape_[idx]) {
            IdxToMultiIdx(ii,idx,nn_);
            std::transform(ii.begin(),ii.end(),ii0.begin(),di.begin(),
                [](const int &a,const int &b)->real_t{return static_cast<real_t>(a-b);});
            F(r,0) = 1.0;
            for (int d = 0; d<n_dim_; ++d) {
                F(r,1+d) = di[d]*di[d];
            }
            if (!shape_.IsSymmetric()) {
                int c = 0;
                for (int d = 0; d<n_dim_; ++d) {
                    F(r,1+n_dim_+d) = di[d];
                    for (int d2 = d+1; d2<n_dim_; ++d2) {
                        F(r,1+2*n_dim_+c) = di[d]*di[d2];
                        ++c;
                    }
                }
            }
            ++r;
        }
    }
    // solve the normal equations
    Eigen::Matrix<real_t,Eigen::Dynamic,Eigen::Dynamic> A = (F.transpose()*F).inverse()*F.transpose();
    // Compute the kernel
    kernel_.resize(n_dim_);
    for (int d = 0; d<n_dim_; ++d) {
        kernel_[d].resize(n_vox_);
        r = 0;
        for (int idx = 0; idx<n_vox_; ++idx) {
            if (shape_[idx]) {
                kernel_[d][idx] = 2.0*A(1+d,r);
                ++r;
            }
        }
    }
    return;
}
