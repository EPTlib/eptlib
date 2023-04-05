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

#include "gtest/gtest.h"

#include "eptlib/linalg/qr.h"

#include <cmath>
#include <complex>
#include <vector>

TEST(LinalgQRGTest, HouseholderReflector) {
    const size_t n = 10;
    std::vector<double> x(n+1);
    std::iota(x.begin(), x.end()-1, 1.0);

    eptlib::linalg::HouseholderReflector(x.begin(), x.end());

    ASSERT_NEAR(x[0], std::sqrt(385.0)+1.0, 1e-12);
    for (size_t i = 1; i<n; ++i) {
        ASSERT_NEAR(x[i], i+1.0, 1e-12);
    }
    ASSERT_NEAR(x[n], 385.0+std::sqrt(385.0), 1e-12);
}

TEST(LinalgQRGTest, HouseholderLeft) {
    const size_t n = 10;
    std::vector<double> u(n+1);
    std::iota(u.begin(), u.end()-1, 1.0);

    std::vector<double> x(n);
    std::copy(u.begin(), u.end()-1, x.begin());

    std::vector<std::complex<double> > x_c(n);
    std::copy(u.begin(), u.end()-1, x_c.begin());

    eptlib::linalg::HouseholderReflector(u.begin(), u.end());
    eptlib::linalg::HouseholderLeft(x.begin(), x.end(), u.begin(), u[n]);
    eptlib::linalg::HouseholderLeft(x_c.begin(), x_c.end(), u.begin(), u[n]);

    ASSERT_NEAR(x[0], -std::sqrt(385.0), 1e-12);
    ASSERT_NEAR(x_c[0].real(), -std::sqrt(385.0), 1e-12);
    ASSERT_NEAR(x_c[0].imag(), 0.0, 1e-12);
    for (size_t i = 1; i<n; ++i) {
        ASSERT_NEAR(x[i], 0.0, 1e-12);
        ASSERT_NEAR(x_c[i].real(), 0.0, 1e-12);
        ASSERT_NEAR(x_c[i].imag(), 0.0, 1e-12);
    }
}

TEST(LinalgQRGTest, QRSolve) {
    const size_t n = 5;

    eptlib::linalg::Matrix<double> A(n, n);
    A(0,0) = 82; A(0,1) = 91; A(0,2) = 13; A(0,3) = 92; A(0,4) = 64;
    A(1,0) = 10; A(1,1) = 28; A(1,2) = 55; A(1,3) = 96; A(1,4) = 97;
    A(2,0) = 16; A(2,1) = 98; A(2,2) = 96; A(2,3) = 49; A(2,4) = 81;
    A(3,0) = 15; A(3,1) = 43; A(3,2) = 92; A(3,3) = 80; A(3,4) = 96;
    A(4,0) = 66; A(4,1) = 4;  A(4,2) = 85; A(4,3) = 94; A(4,4) = 68;

    std::vector<double> u(n);
    std::iota(u.begin(), u.end(), 1.0);

    std::vector<double> b(n);
    for (size_t row = 0; row < n; ++row) {
        b[row] = 0.0;
        for (size_t col = 0; col < n; ++col) {
            b[row] += A(row, col) * u[col];
        }
    }

    std::vector<std::complex<double> > b_c(n);
    std::copy(b.begin(), b.end(), b_c.begin());

    eptlib::linalg::Matrix<double> QR = eptlib::linalg::QRDecomposition(A);
    auto [x, chi] = eptlib::linalg::QRSolve(QR, b);
    auto [x_c, chi_c] = eptlib::linalg::QRSolve(QR, b_c);

    eptlib::linalg::Matrix<double> QR_cp(0,0);
    std::vector<size_t> p_cp(0);
    std::tie(QR_cp, p_cp) = eptlib::linalg::QRDecompositionWithColumnPivoting(A);
    auto [x_cp, chi_cp] = eptlib::linalg::QRSolve(QR_cp, b);
    eptlib::linalg::Permute(x_cp.begin(), x_cp.end(), p_cp);

    ASSERT_DOUBLE_EQ(chi, 0.0);
    ASSERT_DOUBLE_EQ(chi_c, 0.0);
    ASSERT_DOUBLE_EQ(chi_cp, 0.0);
    for (size_t col = 0; col < n; ++col) {
        ASSERT_NEAR(x[col], u[col], 1e-12);
        ASSERT_NEAR(x_c[col].real(), u[col], 1e-12);
        ASSERT_NEAR(x_c[col].imag(), 0.0, 1e-12);
        ASSERT_NEAR(x_cp[col], u[col], 1e-12);
    }
}
