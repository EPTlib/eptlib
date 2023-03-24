/*****************************************************************************
*
*     Program: EPTlib
*     Author: Alessandro Arduino <a.arduino@inrim.it>
*
*  MIT License
*
*  Copyright (c) 2023  Alessandro Arduino
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

#ifndef EPTLIB_POLYNOMIAL_FITTING_H_
#define EPTLIB_POLYNOMIAL_FITTING_H_

#include <tuple>
#include <vector>

#include "eptlib/shape.h"

#include "eptlib/linalg/linalg_householder.h"
#include "eptlib/linalg/linalg_qr.h"
#include "eptlib/linalg/linalg_util.h"

namespace eptlib {

namespace polynomial {

    /**
     * @brief Get the number of three-dimensional monomials of given maximum degree.
     * 
     * @param degree Maximum degree of the monomial.
     * 
     * @return the number of monomials. 
     */
    size_t GetNumberOfMonomials(const size_t degree);

    /**
     * @brief Evaluate all the monomials of given degree.
     * 
     * @param x x-coordinate of where the monomials are evaluted.
     * @param y y-coordinate of where the monomials are evaluted.
     * @param z z-coordinate of where the monomials are evaluted.
     * @param n degree of the monomials.
     * 
     * @return vector with the evaluated monomials.
     * 
     * The monomials are all those written as x^i * y^j * z^k, with
     * i+j+k = n. The elements are ordered as follow,
     * 1. (i=n  , j=0, k=0);
     * 2. (i=n-1, j=1, k=0);
     * 3. (i=n-2, j=2, k=0);
     * ...
     */
    std::vector<double> EvaluateMonomials(const double x, const double y, const double z, const size_t n);

    /**
     * @brief Fill the design matrix for polynomial fitting with a basis of monomials.
     * 
     * @param d0 resolution in meter along direction x.
     * @param d1 resolution in meter along direction y.
     * @param d2 resolution in meter along direction z.
     * @param window mask over which the basis is evaluated.
     * @param degree degree of the fitting polynomial.
     * 
     * @return the design matrix.
     * 
     * The columns of the matrix are ordered as follow,
     * (1, x, y, z, x^2, x*y, y^2, x*z, y*z, z^2, x^3, x^2*y, ...).
     */
    eptlib::linalg::MatrixReal DesignMatrixWithMonomialsBasis(const double d0, const double d1, const double d2,
        const eptlib::Shape &window, const size_t degree);

    /**
     * @brief Permute the columns of a matrix to have all the null columns at the end.
     * 
     * @param A matrix to be permuted.
     * @param n_row number of rows of the matrix.
     * @param n_col number of columns of the matrix.
     * 
     * @return a std::tuple with:
     *     1) a vector of indices defining the permutation;
     *     2) the number of non-null columns in the matrix.
     */
    std::tuple<std::vector<size_t>, size_t> PermuteColumns(eptlib::linalg::MatrixReal *A, const size_t n_row, const size_t n_col);
    
    /**
     * @brief 
     * permutation matrix P
     *   F * x = b
     *   (F*P) * (P*x) = b
     *   F_ = F*P, matrix with permuted columns
     *   x_ = P*x, vector with permuted rows
     *   F_ * x_ = b
     * From column c, F_ has null columns
     *   F' = F_(:, 1:c), shrinked matrix
     * From now on, the actual algorithm
     *   F' * x' = b <- I solve this full rank problem
     *   x_ = [x'; zeros(length(x)-c,1)] <-
     *   x = P*x_ <- original solution vector
     * 
     * The last three steps are the one actually computed.
     * I need to know P (which I use to get F_, so I just
     * have to store it) and c.
     * 
     * It is more convenient to store P as a vector of
     * indices p such that x = x_(p). So, p(1) is the index
     * where the first column of the original problem is moved
     * in the permuted problem, p(2) is the index for the second
     * column, and so on.
     */
    template <typename Scalar>
    std::tuple<std::vector<Scalar>, double> PolynomialFitting(const double d0, const double d1, const double d2,
        const eptlib::Shape &window, const size_t degree, const std::vector<Scalar> &values
//        , const std::vector<double> &weight
        )
        {
        eptlib::linalg::MatrixReal F = DesignMatrixWithMonomialsBasis(d0,d1,d2, window, degree);
        size_t n_col = F.size();
        size_t n_row = F[0].size();
        // permute and shrink the design matrix
        auto [p, n_col_shrinked] = PermuteColumns(&F, n_row,n_col);
        // compute the QR decomposition of F
        eptlib::linalg::MatrixReal QR;
        eptlib::linalg::HouseholderQR(&QR, F ,n_row ,n_col_shrinked);
        // solve the linear system
        std::vector<Scalar> x(n_col, 0.0);
        double chi2n = eptlib::linalg::QRSolve(x.data(), QR, values.data(), n_row, n_col_shrinked);
        // permute the solution
        std::vector<Scalar> y(n_col);
        for (size_t idx = 0; idx<n_col; ++idx) {
            y[p[idx]] = x[idx];
        }
        return {y, chi2n};
    }

}  // namespace polynomial

}  // namespace eptlib

#endif  // EPTLIB_POLYNOMIAL_FITTING_H_
