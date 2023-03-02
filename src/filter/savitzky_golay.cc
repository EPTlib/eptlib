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

#include "eptlib/linalg/linalg_householder.h"
#include "eptlib/linalg/linalg_qr.h"
#include "eptlib/linalg/linalg_util.h"

namespace {

    // Compute the maximum number of monomials in a three-dimensional polynomial
    // of given maximum degree.
    size_t GetNumberOfMonomials(const size_t degree) {
        size_t n_monomials = 0;
        for (size_t deg = 0; deg<=degree; ++deg) {
            n_monomials += (deg+2)*(deg+1)/2;
        }
        return n_monomials;
    }

    // Evaluate in x all the monomials of given degree, namely, all the
    // products x^i * y^j * z^k, for any triple (i, j, k) such that
    // i + j + k = degree.
    std::vector<double> EvaluateMonomials(const double x, const double y, const double z, const size_t degree) {
        std::vector<double> monomials;
        for (size_t k = 0; k<=degree; ++k) {
            for (size_t j = 0; j<=degree-k; ++j) {
                size_t i = degree-j-k;
                double monomial = 1.0;
                for (size_t counter = 0; counter<i; ++counter) {
                    monomial *= x;
                }
                for (size_t counter = 0; counter<j; ++counter) {
                    monomial *= y;
                }
                for (size_t counter = 0; counter<k; ++counter) {
                    monomial *= z;
                }
                monomials.push_back(monomial);
            }
        }
        return monomials;
    }

    // Fill the design matrix for three-dimensional polynomial fitting with the
    // monomial basis (1, x, y, z, x^2, x*y, y^2, x*z, y*z, z^2, x^3, x^2*y, x*y^2, y^3, x^2*z, ...).
    eptlib::linalg::MatrixReal DesignMatrixWithMonomialsBasis(const double d0,
        const double d1, const double d2, const eptlib::Shape &window,
        const size_t degree) {
        // initialise the design matrix
        size_t n_row = window.GetVolume();
        size_t n_col = ::GetNumberOfMonomials(degree);
        eptlib::linalg::MatrixReal F(n_col, std::vector<double>(n_row));
        // fill the design matrix (1, x, y, z, x^2, x*y, y^2, x*z, y*z, z^2, x^3, x^2*y, x*y^2, y^3, x^2*z, ...)
        double x0 = window.GetSize(0)/2*d0;
        double y0 = window.GetSize(1)/2*d1;
        double z0 = window.GetSize(2)/2*d2;
        size_t row = 0;
        for (size_t i2 = 0; i2<window.GetSize(2); ++i2) {
            for (size_t i1 = 0; i1<window.GetSize(1); ++i1) {
                for (size_t i0 = 0; i0<window.GetSize(0); ++i0) {
                    if (window(i0, i1, i2)) {
                        double dx = i0*d0-x0;
                        double dy = i1*d1-y0;
                        double dz = i2*d2-z0;
                        size_t col = 0;
                        for (size_t deg = 0; deg<=degree; ++deg) {
                            std::vector<double> monomials = ::EvaluateMonomials(dx,dy,dz, deg);
                            for (double monomial: monomials) {
                                F[col][row] = monomial;
                                ++col;
                            }
                        }
                        ++row;
                    }
                }
            }
        }
        return F;
    }

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
    window_(window) {
    // define the design matrix F
    eptlib::linalg::MatrixReal F = ::DesignMatrixWithMonomialsBasis(d0,d1,d2, window, degree);
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
    // swap columns to have the essential monomials in the first positions
    F[5].swap(F[6]); // xy <-> yy
    F[6].swap(F[9]); // xy <-> zz
    // remove null columns
    for (size_t col = 0; col<n_col; ++col) {
        if (eptlib::linalg::Norm2(F[col].data(), n_row) == 0.0) {
            --n_col;
            F[col].swap(F[n_col]);
            --col;
        }
    }
    F.resize(n_col);
    // compute the QR decomposition of F
    eptlib::linalg::MatrixReal QR;
    eptlib::linalg::HouseholderQR(&QR, F ,n_row ,n_col);
    // initialize the residuals
    residuals_.resize(n_row, std::vector<double>(n_row, 0.0));
    // solve the linear systems (F * a_k = e_k, for k = 1,...,n_row)
    std::vector<double> a(n_col);
    std::vector<double> b(n_row, 0.0);
    for (int k = 0; k<n_row; ++k) {
        b.assign(n_row, 0.0);
        b[k] = 1.0;
        eptlib::linalg::QRSolve(a.data(), QR, b.data(), n_row, n_col);
        // assign the kernel coefficients
        zero_order_derivative_     [k] =     a[0];
        first_order_derivative_ [0][k] =     a[1];
        first_order_derivative_ [1][k] =     a[2];
        first_order_derivative_ [2][k] =     a[3];
        second_order_derivative_[0][k] = 2.0*a[4];
        second_order_derivative_[1][k] = 2.0*a[5];
        second_order_derivative_[2][k] = 2.0*a[6];
        // compute the residual (r_k = F * a_k - e_k)
        for (int col = 0; col<n_col; ++col) {
            for (int row = 0; row<n_row; ++row) {
                residuals_[k][row] += F[col][row] * a[col];
            }
        }
        residuals_[k][k] -= 1.0;
    }
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
                std::exp(1.0i * phi);
            }
        );
        return exp_iphi;
    }

}  //

// SavitzkyGolay get the correct filter function for a wrapped phase map
std::function<double(const std::vector<double>&)> eptlib::filter::SavitzkyGolay::
GetFilterWrappedPhase(const eptlib::filter::DifferentialOperator differential_operator) const {
    switch (differential_operator) {
    case DifferentialOperator::Field:
        return [&](const std::vector<double> &crop) -> double {
            return std::log(this->ZeroOrderDerivative(ExpIVector(crop))).imag();
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

eptlib::EPTlibError eptlib::filter::SavitzkyGolay::
ApplyWrappedPhase(const eptlib::filter::DifferentialOperator differential_operator, eptlib::Image<double> *dst, const eptlib::Image<double> &src) const {
    auto filter = GetFilterWrappedPhase(differential_operator);
    return MovingWindow(dst, src, window_, filter);
}
