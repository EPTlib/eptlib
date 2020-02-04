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

#ifndef EPTLIB_FINITE_DIFFERENCES_H_
#define EPTLIB_FINITE_DIFFERENCES_H_

#include <cassert>

#include <Eigen/Dense>

#include "eptlib/config.h"
#include "eptlib/shape.h"
#include "eptlib/util.h"

namespace eptlib {

/**
 * Class for the application of finite differences in computing the Laplacian.
 * 
 * It implements the Savitzky-Golay filter.
 */
class FDLaplacianKernel {
    public:
        /**
         * Constructor.
         * 
         * @param shape mask over which apply the finite difference scheme.
         */
        FDLaplacianKernel(const Shape &shape);

        /**
         * Apply the FD Laplacian filter to an input field.
         * 
         * @tparam NumType numeric typename.
         * @tparam T,U iterator typenames.
         * 
         * @param[out] dst pointer to the output destination.
         * @param[in] src pointer to the input source.
         * @param[in] u_nn number of voxels in each direction of the field.
         * @param[in] u_dd size of voxels in each direction.
         * 
         * @return a Success or Unknown error.
         */
        template <typename NumType, typename T, typename U>
        EPTlibError_t ApplyFilter(NumType *dst, const NumType *src, const T &u_nn, const U &u_dd);
    private:
        /// Number of spatial dimensions.
        int n_dim_;
        /// Total number of voxels.
        int n_vox_;
        /// Number of voxels in each direction.
        std::vector<int> nn_;
        /// Shape of the kernel for Laplacian approximation.
        Shape shape_;
        /// Kernel for Laplacian approximation.
        std::vector<std::vector<real_t>> kernel_;
};

// ---------------------------------------------------------------------------
// -------------------------  Implementation detail  -------------------------
// ---------------------------------------------------------------------------

template <typename NumType, typename T, typename U>
EPTlibError_t FDLaplacianKernel::
ApplyFilter(NumType *dst, const NumType *src, const T &u_nn, const U &u_dd) {
    // check input coherence
    assert(u_nn.size()==n_dim_);
    // variables declaration
    const int u_n_vox = std::accumulate(u_nn.begin(),u_nn.end(),1,std::multiplies<int>());
    std::vector<int> ii(n_dim_);
    std::vector<int> ii_l(n_dim_);
    std::vector<int> ii_g(n_dim_);
    std::vector<NumType> field_crop(n_vox_);
    bool kernel_in_domain;
    // loop over field voxels
    for (int idx = 0; idx<u_n_vox; ++idx) {
        IdxToMultiIdx(ii,idx,u_nn);
        kernel_in_domain = true;
        // inner loop over kernel voxels
        for (int idx_k = 0; idx_k<n_vox_; ++idx_k) {
            if (shape_[idx_k]) {
                // check if the element fall within the field domain
                IdxToMultiIdx(ii_l,idx_k,nn_);
                for (int d = 0; d<n_dim_; ++d) {
                    ii_l[d] -= nn_[d]/2;
                    ii_g[d] = ii[d]+ii_l[d];
                    if (ii_g[d]<0 || ii_g[d]>=u_nn[d]) {
                        // outside field domain
                        kernel_in_domain = false;
                        break;
                    }
                }
                if (kernel_in_domain) {
                    // set the `field_crop' to the field value
                    int idx_g = MultiIdxToIdx(ii_g,u_nn);
                    field_crop[idx_k] = src[idx_g];
                } else {
                    // outside field domain
                    break;
                }
            } else {
                // set to zero the `field_crop'
                field_crop[idx_k] = 0.0;
            }
        }
        // set to zero the `dst'
        dst[idx] = 0.0;
        // if the previous loop ended with `kernel_in_domain==true'...
        if (kernel_in_domain) {
            // ...then apply the filter
            for (int d = 0; d<n_dim_; ++d) {
                dst[idx] += std::inner_product(kernel_[d].begin(),kernel_[d].end(),field_crop.begin(),
                    static_cast<NumType>(0.0))/u_dd[d]/u_dd[d];
            }
        }
    }
    return EPTlibError::Success;
}

}  // namespace eptlib

#endif  // EPTLIB_FINITE_DIFFERENCES_H_
