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

#ifndef EPTLIB_LINALG_REGRESSION_H_
#define EPTLIB_LINALG_REGRESSION_H_

#include <tuple>

#include "eptlib/linalg/qr.h"
#include "eptlib/linalg/matrix.h"
#include "eptlib/linalg/vector.h"

namespace eptlib {

namespace linalg {

    /**
     * @brief Solve a linear regression problem using the QR decomposition.
     * 
     * @tparam Scalar numerical type of the forcing term entries.
     * 
     * @param A design matrix of the regression problem.
     * @param b forcing term of the regression problem.
     * 
     * @return a std::tuple with:
     *     1) a vector solving the linear regression;
     *     2) the normalized chi-squared statistic of the linear regression.
     *     If `b` is a complex-valued vector, the statistic of the real and
     *     imaginary parts of the residual are stored as the real and the
     *     imaginary parts of a complex number.
     */
    template <typename Scalar>
    std::tuple<std::vector<Scalar>, Scalar> LinearRegression(const eptlib::linalg::Matrix<double> &A, const std::vector<Scalar> &b) {
        auto [QR, p] = eptlib::linalg::QRDecomposition(A);
        auto [x, chi] = eptlib::linalg::QRSolve(QR, b);
        eptlib::linalg::Permute(x.begin(), x.end(), p);
        size_t m = A.GetNRow();
        size_t r = eptlib::linalg::QRGetRank(QR);
        Scalar chi2n;
        if constexpr (std::is_same_v<Scalar, std::complex<double> >) {
            double chi_r = std::real(chi);
            double chi_i = std::imag(chi);
            double chi2n_r = chi_r > 0.0 ? chi_r * chi_r / (m - r) : 0.0;
            double chi2n_i = chi_i > 0.0 ? chi_i * chi_i / (m - r) : 0.0;
            chi2n = Scalar(chi2n_r, chi2n_i);
        } else {
            chi2n = chi > 0.0 ? chi * chi / (m - r) : 0.0;
        }
        return {x, chi2n};
    }

}  // namespace linalg

}  // namespace eptlib

#endif  // EPTLIB_LINALG_REGRESSION_H_
