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

#include "eptlib/filter/savitzky_golay.h"

#include <complex>
#include <tuple>

#include "eptlib/linalg/linalg_householder.h"
#include "eptlib/linalg/linalg_qr.h"
#include "eptlib/linalg/linalg_util.h"

#include "eptlib/polynomial/fitting.h"

namespace {

    // Solve a linear equation with R^T as matrix of coefficients.
    void RTSolve(double *a, const eptlib::linalg::MatrixReal &QR, const double *b, const size_t n) {
        a[0] = b[0] / QR[0][0];
        for (int i = 1; i<n; ++i) {
            double s = 0;
            for (int j = 0; j<i; ++j) {
                s += a[j] * QR[i][j];
            }
            a[i] = (b[i] - s) / QR[i][i];
        }
        return;
    }

}  //

// SavitzkyGolay constructor
eptlib::filter::SavitzkyGolay::
SavitzkyGolay(const double d0, const double d1, const double d2,
    const eptlib::Shape &window, const size_t degree) :
    window_(window), degree_(degree) {
    dd_[0] = d0;
    dd_[1] = d1;
    dd_[2] = d2;
    eptlib::linalg::MatrixReal F = eptlib::polynomial::DesignMatrixWithMonomialsBasis(d0,d1,d2, window, degree);
    size_t n_col = F.size();
    size_t n_row = F[0].size();
    // initialise the kernel coefficients
    zero_order_derivative_     .resize(window_.GetVolume());
    first_order_derivative_ [0].resize(window_.GetVolume());
    first_order_derivative_ [1].resize(window_.GetVolume());
    first_order_derivative_ [2].resize(window_.GetVolume());
    second_order_derivative_[0].resize(window_.GetVolume());
    second_order_derivative_[1].resize(window_.GetVolume());
    second_order_derivative_[2].resize(window_.GetVolume());
    // initialize the residuals
    residuals_.resize(n_row, std::vector<double>(n_row, 0.0));
    // solve the linear systems (F * a_k = e_k, for k = 1,...,n_row)
    std::vector<double> a(n_col);
    std::vector<double> b(n_row);
    for (int k = 0; k<n_row; ++k) {
        b.assign(n_row, 0.0);
        b[k] = 1.0;
        std::vector<double> a;
        std::tie(a, std::ignore) = eptlib::polynomial::PolynomialFitting(d0,d1,d2, window_, degree_, b);
        // assign the kernel coefficients
        zero_order_derivative_     [k] =     a[0];
        first_order_derivative_ [0][k] =     a[1];
        first_order_derivative_ [1][k] =     a[2];
        first_order_derivative_ [2][k] =     a[3];
        second_order_derivative_[0][k] = 2.0*a[4];
        second_order_derivative_[1][k] = 2.0*a[6];
        second_order_derivative_[2][k] = 2.0*a[9];
        // compute the residual (r_k = F * a_k - e_k)
        for (int col = 0; col<n_col; ++col) {
            for (int row = 0; row<n_row; ++row) {
                residuals_[k][row] += F[col][row] * a[col];
            }
        }
        residuals_[k][k] -= 1.0;
    }
    // compute the QR decomposition of F
    auto [p, n_col_shrinked] = eptlib::polynomial::PermuteColumns(&F, n_row,n_col);
    eptlib::linalg::MatrixReal QR;
    eptlib::linalg::HouseholderQR(&QR, F ,n_row ,n_col);
    // compute the variance coefficients of zero and first order derivatives
    for (int der = 0; der<4; ++der) {
        a.assign(n_col, 0.0);
        a[der] = 1.0;
        RTSolve(a.data(), QR, a.data(), n_col);
        variance_coefficients_[der] = eptlib::linalg::Norm2(a.data(), n_col);
    }
    // compute the variance coefficients of second order derivatives
    for (int der = 4; der<7; ++der) {
        a.assign(n_col, 0.0);
        a[der] = 2.0;
        RTSolve(a.data(), QR, a.data(), n_col);
        variance_coefficients_[der] = eptlib::linalg::Norm2(a.data(), n_col);
    }
    // compute the variance coefficients of laplacian
    a.assign(n_col, 0.0);
    a[4] = 2.0;
    a[5] = 2.0;
    a[6] = 2.0;
    RTSolve(a.data(), QR, a.data(), n_col);
    variance_coefficients_[7] = eptlib::linalg::Norm2(a.data(), n_col);
    // rescale the variance coefficients with n_row-n_col
    for (int idx = 0; idx<variance_coefficients_.size(); ++idx) {
        variance_coefficients_[idx] /= n_row - n_col;
    }
    return;
}

// SavitzkyGolay virtual destructor
eptlib::filter::SavitzkyGolay::
~SavitzkyGolay() {
    return;
}

namespace {

    using namespace std::literals::complex_literals;

    // Compute the exponential e^{i * phi} for each element phi of a real-valued vector.
    std::vector<std::complex<double> > ExpIVector(const std::vector<double> &phi) {
        std::vector<std::complex<double> > exp_iphi(phi.size());
        std::transform(phi.begin(), phi.end(), exp_iphi.begin(),
            [](const double &phi) -> std::complex<double> {
                return std::exp(1.0i * phi);
            }
        );
        return exp_iphi;
    }

}  //

// SavitzkyGolay get the correct filter function for a wrapped phase map
std::function<double(const std::vector<double>&)> eptlib::filter::SavitzkyGolay::
GetFilterWrappedPhase(const eptlib::DifferentialOperator differential_operator) const {
    switch (differential_operator) {
    case DifferentialOperator::Field:
        return [&](const std::vector<double> &crop) -> double {
            return std::arg(this->ZeroOrderDerivative(ExpIVector(crop)));
        };
    case DifferentialOperator::GradientX:
        return [&](const std::vector<double> &crop) -> double {
            auto exp_icrop = ExpIVector(crop);
            std::complex<double> e = this->ZeroOrderDerivative(exp_icrop);
            std::complex<double> de_dx = this->FirstOrderDerivative(0, exp_icrop);
            return (de_dx/e).imag();
        };
    case DifferentialOperator::GradientY:
        return [&](const std::vector<double> &crop) -> double {
            auto exp_icrop = ExpIVector(crop);
            std::complex<double> e = this->ZeroOrderDerivative(exp_icrop);
            std::complex<double> de_dy = this->FirstOrderDerivative(1, exp_icrop);
            return (de_dy/e).imag();
        };
    case DifferentialOperator::GradientZ:
        return [&](const std::vector<double> &crop) -> double {
            auto exp_icrop = ExpIVector(crop);
            std::complex<double> e = this->ZeroOrderDerivative(exp_icrop);
            std::complex<double> de_dz = this->FirstOrderDerivative(2, exp_icrop);
            return (de_dz/e).imag();
        };
    case DifferentialOperator::GradientXX:
        return [&](const std::vector<double> &crop) -> double {
            auto exp_icrop = ExpIVector(crop);
            std::complex<double> e = this->ZeroOrderDerivative(exp_icrop);
            std::complex<double> de_dx = this->FirstOrderDerivative(0, exp_icrop);
            std::complex<double> d2e_dx2 = this->SecondOrderDerivative(0, exp_icrop);
            return (d2e_dx2/e - de_dx*de_dx/e/e).imag();
        };
    case DifferentialOperator::GradientYY:
        return [&](const std::vector<double> &crop) -> double {
            auto exp_icrop = ExpIVector(crop);
            std::complex<double> e = this->ZeroOrderDerivative(exp_icrop);
            std::complex<double> de_dy = this->FirstOrderDerivative(1, exp_icrop);
            std::complex<double> d2e_dy2 = this->SecondOrderDerivative(1, exp_icrop);
            return (d2e_dy2/e - de_dy*de_dy/e/e).imag();
        };
    case DifferentialOperator::GradientZZ:
        return [&](const std::vector<double> &crop) -> double {
            auto exp_icrop = ExpIVector(crop);
            std::complex<double> e = this->ZeroOrderDerivative(exp_icrop);
            std::complex<double> de_dz = this->FirstOrderDerivative(2, exp_icrop);
            std::complex<double> d2e_dz2 = this->SecondOrderDerivative(2, exp_icrop);
            return (d2e_dz2/e - de_dz*de_dz/e/e).imag();
        };
    case DifferentialOperator::Laplacian:
        return [&](const std::vector<double> &crop) -> double {
            auto exp_icrop = ExpIVector(crop);
            std::complex<double> e = this->ZeroOrderDerivative(exp_icrop);
            std::complex<double> de_dx = this->FirstOrderDerivative(0, exp_icrop);
            std::complex<double> de_dy = this->FirstOrderDerivative(1, exp_icrop);
            std::complex<double> de_dz = this->FirstOrderDerivative(2, exp_icrop);
            std::complex<double> lapl_e = this->Laplacian(exp_icrop);
            return (lapl_e/e - (de_dx*de_dx + de_dy*de_dy + de_dz*de_dz)/e/e).imag();
        };
    default:
        return [ ](const std::vector<double> &crop) -> double {
            return 0.0;
        };
    };
}

// SavitzkyGolay apply filter with wrapped phase
eptlib::EPTlibError eptlib::filter::SavitzkyGolay::
ApplyWrappedPhase(const eptlib::DifferentialOperator differential_operator, eptlib::Image<double> *dst, const eptlib::Image<double> &src) const {
    auto filter = GetFilterWrappedPhase(differential_operator);
    return MovingWindow(dst, src, window_, filter);
}
