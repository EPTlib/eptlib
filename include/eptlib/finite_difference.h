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

#ifndef EPTLIB_FINITE_DIFFERENCES_H_
#define EPTLIB_FINITE_DIFFERENCES_H_

#include <array>
#include <vector>

#include "eptlib/shape.h"
#include "eptlib/util.h"
#include "eptlib/linalg/linalg_util.h"

namespace eptlib {

/**
 * Differential operators that can be approximated with Savitzky-Golay.
 */
enum class DifferentialOperator {
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
};

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
         * @brief Apply the FD filter to an input field and compute the fitting quality.
         * 
         * @tparam NumType numeric typename.
         * 
         * @param[in] diff_op differential operator type.
         * @param[out] dst pointer to the output destination.
         * @param[out] chi2 pointer to the output fitting quality.
         * @param[in] src pointer to the input source.
         * @param[in] nn number of voxels in each direction.
         * @param[in] dd size of voxels in each direction.
         * 
         * @return a Success or Unknown error.
         */
        template <typename NumType>
        EPTlibError ApplyChi2(const DifferentialOperator diff_op,NumType *dst,real_t *chi2,
            const NumType *src,const std::array<int,NDIM> &nn,const std::array<double,NDIM> &dd) const;
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
        EPTlibError Apply(const DifferentialOperator diff_op, NumType *dst,
            const NumType *src, const std::array<int,NDIM> &nn, const std::array<double,NDIM> &dd) const;
        /**
         * Apply the kernel for first order derivatives.
         * 
         * @tparam NumType numeric typename.
         * 
         * @param d derivative direction.
         * @param field_crop field values within the computation kernel.
         * @param dd size of voxels in each direction.
         * 
         * @return the approximated derivative.
         */
        template <typename NumType>
        NumType FirstOrder(const int d, const std::vector<NumType> &field_crop, const std::array<double,NDIM> &dd) const;
        /**
         * Apply the kernel for second order derivatives.
         * 
         * @tparam NumType numeric typename.
         * 
         * @param d derivative direction.
         * @param field_crop field values within the computation kernel.
         * @param dd size of voxels in each direction.
         * 
         * @return the approximated derivative.
         */
        template <typename NumType>
        NumType SecondOrder(const int d, const std::vector<NumType> &field_crop, const std::array<double,NDIM> &dd) const;
        /**
         * Apply the kernel for Laplacian computation.
         * 
         * @tparam NumType numeric typename.
         * 
         * @param field_crop field values within the computation kernel.
         * @param dd size of voxels in each direction.
         * 
         * @return the approximated Laplacian.
         */
        template <typename NumType>
        NumType Laplacian(const std::vector<NumType> &field_crop, const std::array<double,NDIM> &dd) const;
        /**
         * Return a const reference to the kernel shape.
         * 
         * @return a const reference to the kernel shape.
         */
        const Shape& GetShape() const;
        /**
         * Return a const reference to the kernel for Laplacian.
         * 
         * @return a const reference to the kernel for Laplacian.
         */
        const std::array<std::vector<double>,NDIM>& GetLaplKernel() const;
        /**
         * Return a const reference to the kernel for gradient.
         * 
         * @return a const reference to the kernel for gradient.
         */
        const std::array<std::vector<double>,NDIM>& GetGradKernel() const;
    private:
        /// Shape of the kernel for Laplacian approximation.
        Shape shape_;
        /// Total number of voxels.
        int m_vox_;
        /// Number of fitting unknowns.
        int n_unk_;
        /// QR decomposition of the design matrix.
        linalg::MatrixReal qr_;
        /// Kernel for Laplacian approximation.
        std::array<std::vector<double>,NDIM> lapl_kernel_;
        /// Kernel for gradient approximation.
        std::array<std::vector<double>,NDIM> grad_kernel_;
};

/**
 * Apply the FD filter to a wrapped phase input field.
 * 
 * @param[in] diff_op differential operator type.
 * @param[out] dst pointer to the output destination.
 * @param[in] src pointer to the input source.
 * @param[in] nn number of voxels in each direction.
 * @param[in] dd size of voxels in each direction.
 * @param[in] fd_filter filter that computes the derivatives.
 * 
 * @return a Success or Unknown error.
 */
EPTlibError WrappedPhaseDerivative(const DifferentialOperator diff_op,
    double *dst, const double *src, const std::array<int,NDIM> &nn,
    const std::array<double,NDIM> &dd, const FDSavitzkyGolayFilter &fd_filter);

/**
 * Apply the FD filter to a wrapped phase input field and compute the fitting quality.
 * 
 * @param[in] diff_op differential operator type.
 * @param[out] dst pointer to the output destination.
 * @param[out] chi2 pointer to the output fitting quality.
 * @param[in] src pointer to the input source.
 * @param[in] nn number of voxels in each direction.
 * @param[in] dd size of voxels in each direction.
 * @param[in] fd_filter filter that computes the derivatives.
 * 
 * @return a Success or Unknown error.
 */
EPTlibError WrappedPhaseDerivativeChi2(const DifferentialOperator diff_op,
    double *dst, real_t *chi2, const double *src, const std::array<int,NDIM> &nn,
    const std::array<double,NDIM> &dd, const FDSavitzkyGolayFilter &fd_filter);

}  // namespace eptlib

#endif  // EPTLIB_FINITE_DIFFERENCES_H_
