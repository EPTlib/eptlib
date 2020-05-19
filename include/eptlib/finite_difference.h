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

#include <array>
#include <vector>

#include "eptlib/shape.h"
#include "eptlib/util.h"

namespace eptlib {

/**
 * Differential operators that can be approximated with Savitzky-Golay.
 */
typedef enum DifferentialOperator {
    /// First order derivative along X
    GradientX = 0,
    /// First order derivative along Y
    GradientY,
    /// First order derivative along z
    GradientZ,
    /// Second order derivative along XX
    GradientXX,
    /// Second order derivative along YY
    GradientYY,
    /// Second order derivative along ZZ
    GradientZZ,
    /// Laplacian
    Laplacian,
} DifferentialOperator_t;

/**
 * Class for the application of finite differences.
 * 
 * It implements the Savitzky-Golay filter.
 */
class FDSavitzkyGolayFilter {
    public:
        /**
         * Constructor.
         * 
         * @param shape mask over which apply the finite difference scheme.
         */
        FDSavitzkyGolayFilter(const Shape &shape);
        /**
         * Apply the FD filter to an input field.
         * 
         * @tparam NumType numeric typename.
         * 
         * @param[in] diff_op differential operator type.
         * @param[out] dst pointer to the output destination.
         * @param[in] src pointer to the input source.
         * @param[in] nn number of voxels in each direction.
         * @param[in] dd size of voxels in each direction.
         * 
         * @return a Success or Unknown error.
         */
        template <typename NumType>
        EPTlibError_t Apply(const DifferentialOperator_t diff_op, NumType *dst,
            const NumType *src, const std::array<int,NDIM> &nn, const std::array<double,NDIM> &dd);
    private:
        /// Shape of the kernel for Laplacian approximation.
        Shape shape_;
        /// Total number of voxels.
        int m_vox_;
        /// Kernel for Laplacian approximation.
        std::array<std::vector<double>,NDIM> lapl_kernel_;
        /// Kernel for gradient approximation.
        std::array<std::vector<double>,NDIM> grad_kernel_;
        /// Apply the kernel for first order derivatives.
        template <typename NumType>
        NumType FirstOrder(const int d, const std::vector<NumType> &field_crop, const std::array<double,NDIM> &dd);
        /// Apply the kernel for second order derivatives.
        template <typename NumType>
        NumType SecondOrder(const int d, const std::vector<NumType> &field_crop, const std::array<double,NDIM> &dd);
        /// Apply the kernel for Laplacian computation.
        template <typename NumType>
        NumType Laplacian(const std::vector<NumType> &field_crop, const std::array<double,NDIM> &dd);
};

}  // namespace eptlib

#endif  // EPTLIB_FINITE_DIFFERENCES_H_
