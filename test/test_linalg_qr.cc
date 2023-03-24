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
