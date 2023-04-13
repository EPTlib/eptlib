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

#include <algorithm>
#include <complex>
#include <tuple>

#include "eptlib/polynomial/fitting.h"

namespace {

    // Solve a linear equation with R^T as matrix of coefficients.
    template <typename Scalar>
    std::vector<Scalar> RTSolve(const eptlib::linalg::Matrix<double> &QR, std::vector<Scalar> b) {
        const size_t rank = eptlib::linalg::QRGetRank(QR);
        std::vector<Scalar> x(rank);
        x[0] = b[0] / QR(0,0);
        for (size_t i = 1; i < rank; ++i) {
            Scalar s = 0.0;
            for (size_t j = 0; j < i; ++j) {
                s += QR(j, i) * x[j];
            }
            x[i] = (b[i] - s) / QR(i, i);
        }
        return x;
    }

}  //

// SavitzkyGolay constructor
eptlib::filter::SavitzkyGolay::
SavitzkyGolay(const double d0, const double d1, const double d2,
    const eptlib::Shape &window, const size_t degree) :
    window_(window),
    degree_(degree),
    residuals_(0,0) {
    dd_[0] = d0;
    dd_[1] = d1;
    dd_[2] = d2;
    zero_order_derivative_     .resize(window_.GetVolume());
    first_order_derivative_ [0].resize(window_.GetVolume());
    first_order_derivative_ [1].resize(window_.GetVolume());
    first_order_derivative_ [2].resize(window_.GetVolume());
    second_order_derivative_[0].resize(window_.GetVolume());
    second_order_derivative_[1].resize(window_.GetVolume());
    second_order_derivative_[2].resize(window_.GetVolume());
    // get coordinates of window voxels
    std::vector<double> x(window_.GetVolume());
    std::vector<double> y(window_.GetVolume());
    std::vector<double> z(window_.GetVolume());
    const double x0 = window_.GetSize(0)/2*dd_[0];
    const double y0 = window_.GetSize(1)/2*dd_[1];
    const double z0 = window_.GetSize(2)/2*dd_[2];
    size_t idx = 0;
    for (size_t k = 0; k < window_.GetSize(2); ++k) {
        for (size_t j = 0; j < window_.GetSize(1); ++j) {
            for (size_t i = 0; i < window_.GetSize(0); ++i) {
                if (window_(i,j,k)) {
                    x[idx] = i*dd_[0] - x0;
                    y[idx] = j*dd_[1] - y0;
                    z[idx] = k*dd_[2] - z0;
                    ++idx;
                }
            }
        }
    }
    // compute the design matrix
    const eptlib::linalg::Matrix<double> F = eptlib::polynomial::DesignMatrixWithMonomialsBasis(x, y, z, degree_);
    const size_t n_row = F.GetNRow();
    const size_t n_col = F.GetNCol();
    // factorize the design matrix
    eptlib::linalg::Matrix<double> QR;
    std::vector<size_t> p;
    std::tie(QR, p) = eptlib::linalg::QRDecomposition(F);
    const size_t rank = eptlib::linalg::QRGetRank(QR);
    // compute the kernel coefficients
    residuals_ = eptlib::linalg::Matrix<double>(n_row - rank, n_row);
    std::vector<double> b(window_.GetVolume(), 0.0);
    for (size_t k = 0; k < b.size(); ++k) {
        std::fill(b.begin(), b.end(), 0.0);
        b[k] = 1.0;
        // solve the linear system
        for (size_t col = 0; col < rank; ++col) {
            eptlib::linalg::HouseholderLeft(b.begin()+col, b.end(), QR.begin(col)+col+1, QR(n_row+1,col));
        }
        std::vector<double> c1(rank);
        std::copy(b.begin(), b.begin()+rank, c1.begin());
        std::vector<double> a = SolveTriU(QR, c1);
        a.resize(n_col, 0.0);
        eptlib::linalg::Permute(a.begin(), a.end(), p);
        // store the residual
        std::copy(b.begin()+rank, b.end(), residuals_.begin(k));
        // store the coefficients
        zero_order_derivative_     [k] =     a[0];
        first_order_derivative_ [0][k] =     a[1];
        first_order_derivative_ [1][k] =     a[2];
        first_order_derivative_ [2][k] =     a[3];
        second_order_derivative_[0][k] = 2.0*a[4];
        second_order_derivative_[1][k] = 2.0*a[6];
        second_order_derivative_[2][k] = 2.0*a[9];
    }
    // compute the variance coefficients of zero and first order derivatives
    std::vector<size_t> derivative_indices = {0,1,2,3,4,6,9};
    std::vector<double> w(n_col);
    std::vector<double> w_p(n_col);
    for (size_t der = 0; der < 4; ++der) {
        std::fill(w.begin(), w.end(), 0.0);
        w[derivative_indices[der]] = 1.0;
        for (size_t col = 0; col < n_col; ++col) {
            w_p[p[col]] = w[col];
        }
        std::vector<double> r = RTSolve(QR, w_p);
        variance_coefficients_[der] = eptlib::linalg::Norm2(r.begin(), r.end());
    }
    // compute the variance coefficients of second order derivatives
    for (size_t der = 4; der < 7; ++der) {
        std::fill(w.begin(), w.end(), 0.0);
        w[derivative_indices[der]] = 2.0;
        for (size_t col = 0; col < n_col; ++col) {
            w_p[p[col]] = w[col];
        }
        std::vector<double> r = RTSolve(QR, w_p);
        variance_coefficients_[der] = eptlib::linalg::Norm2(r.begin(), r.end());
    }
    // compute the variance coefficients of laplacian
    {
        std::fill(w.begin(), w.end(), 0.0);
        w[derivative_indices[4]] = 2.0;
        w[derivative_indices[5]] = 2.0;
        w[derivative_indices[6]] = 2.0;
        for (size_t col = 0; col < n_col; ++col) {
            w_p[p[col]] = w[col];
        }
        std::vector<double> r = RTSolve(QR, w_p);
        variance_coefficients_[7] = eptlib::linalg::Norm2(r.begin(), r.end());
    }
    // reduce the variance coefficients
    for (int idx = 0; idx<variance_coefficients_.size(); ++idx) {
        variance_coefficients_[idx] *= variance_coefficients_[idx] / (n_row - rank);
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
