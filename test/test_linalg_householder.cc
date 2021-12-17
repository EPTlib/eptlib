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

#include "gtest/gtest.h"

#include "eptlib/linalg/linalg_householder.h"
#include "eptlib/linalg/linalg_qr.h"

#include <cmath>
#include <complex>
#include <vector>

#include "eptlib/util.h"
#include "eptlib/linalg/linalg_util.h"

using namespace eptlib;
using namespace eptlib::linalg;

TEST(LinalgHouseholderGTest,HouseholderReflector) {
    const int n = 10;
    std::vector<double> x(n+1);
    for (int i = 0; i<n; ++i) {
        x[i] = i+1;
    }
    //
    HouseholderReflector(x.data(),n);
    //
    std::vector<double> x_ref(n+1);
    x_ref[0] = 1+std::sqrt(385.0);
    for (int i = 1; i<n; ++i) {
        x_ref[i] = i+1;
    }
    x_ref[n] = x_ref[0]*std::sqrt(385.0);
    //
    for (int i = 0; i<n+1; ++i) {
        ASSERT_NEAR(x[i],x_ref[i],1e-12);
    }
}
TEST(LinalgHouseholderGTest,HouseholderLeft) {
    const int n = 10;
    std::vector<double> u(n+1);
    std::vector<double> x(n);
    std::vector<std::complex<double> > x_c(n);
    for (int i = 0; i<n; ++i) {
        u[i] = i+1;
        x[i] = i+1;
        x_c[i] = i+1;
    }
    //
    HouseholderReflector(u.data(),n);
    HouseholderLeft(x.data(),u.data(),n);
    HouseholderLeft(x_c.data(),u.data(),n);
    //
    std::vector<double> x_ref(n);
    x_ref[0] = -std::sqrt(385.0);
    for (int i = 1; i<n; ++i) {
        x_ref[i] = 0.0;
    }
    //
    for (int i = 0; i<n+1; ++i) {
        ASSERT_NEAR(x[i],x_ref[i],1e-12);
        ASSERT_NEAR(x_c[i].real(),x_ref[i],1e-12);
        ASSERT_NEAR(x_c[i].imag(),0.0,1e-12);
    }
}

TEST(LinalgQRGTest,QRSolve) {
    const int n = 5;
    MatrixReal A(n,std::vector<double>(n));
    std::vector<double> x_ref(n);
    for (int j = 0; j<n; ++j) {
        x_ref[j] = j+1;
    }
    A[0][0] = 82; A[0][1] = 91; A[0][2] = 13; A[0][3] = 92; A[0][4] = 64;
    A[1][0] = 10; A[1][1] = 28; A[1][2] = 55; A[1][3] = 96; A[1][4] = 97;
    A[2][0] = 16; A[2][1] = 98; A[2][2] = 96; A[2][3] = 49; A[2][4] = 81;
    A[3][0] = 15; A[3][1] = 43; A[3][2] = 92; A[3][3] = 80; A[3][4] = 96;
    A[4][0] = 66; A[4][1] = 4;  A[4][2] = 85; A[4][3] = 94; A[4][4] = 68;
    //
    std::vector<double> b(n);
    for (int i = 0; i<n; ++i) {
        b[i] = 0.0;
        for (int j = 0; j<n; ++j) {
            b[i] += A[j][i]*x_ref[j];
        }
    }
    std::vector<std::complex<double> > b_c(n);
    for (int i = 0; i<n; ++i) {
        b_c[i] = b[i];
    }
    //
    std::vector<double> x(n);
    std::vector<std::complex<double> > x_c(n);
    MatrixReal qr;
    HouseholderQR(&qr,A,n,n);
    QRSolve(x.data(),qr,b.data(),n,n);
    QRSolve(x_c.data(),qr,b_c.data(),n,n);
    //
    for (int i = 0; i<n; ++i) {
        ASSERT_NEAR(x[i],x_ref[i],1e-12);
        ASSERT_NEAR(x_c[i].real(),x_ref[i],1e-12);
        ASSERT_NEAR(x_c[i].imag(),0.0,1e-12);
    }
}
