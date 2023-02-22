/*****************************************************************************
*
*     Program: EPTlib
*     Author: Alessandro Arduino <a.arduino@inrim.it>
*
*  MIT License
*
*  Copyright (c) 2020-2023  Alessandro Arduino
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

#include "eptlib/linalg/linalg_util.h"

#include <cmath>
#include <complex>
#include <vector>

#include "eptlib/util.h"

TEST(LinalgUtilGTest,Norm2) {
    const int n = 10;
    std::vector<double> x(n);
    for (int i = 0; i<n; ++i) {
        x[i] = i+1;
    }
    //
    double norm = eptlib::linalg::Norm2(x.data(),n);
    //
    double ref = std::sqrt(385.0);
    ASSERT_NEAR(norm,ref,1e-12);
}

TEST(LinalgUtilGTest,Dot) {
    const int n = 10;
    std::vector<double> x(n);
    std::vector<std::complex<double> > y(n);
    for (int i = 0; i<n; ++i) {
        x[i] = i+1;
        y[i] = i+1.0;
    }
    //
    double dot = eptlib::linalg::Dot(x.data(),x.data(),n);
    std::complex<double> dot_c = eptlib::linalg::Dot(y.data(),x.data(),n);
    //
    double ref = 385.0;
    ASSERT_NEAR(dot,ref,1e-12);
    ASSERT_NEAR(dot_c.real(),ref,1e-12);
    ASSERT_NEAR(dot_c.imag(),0.0,1e-12);
}

TEST(LinalgUtilGTest,MaxAbs) {
    const int n = 10;
    std::vector<double> x(n);
    for (int i = 0; i<n; ++i) {
        x[i] = -(i+1);
    }
    //
    double maxabs = eptlib::linalg::MaxAbs(x.data(),n);
    //
    double ref = 10.0;
    ASSERT_NEAR(maxabs,ref,1e-12);
}

TEST(LinalgUtilGTest,SolveDiag) {
    const int n = 10;
    eptlib::linalg::MatrixReal A(n,std::vector<double>(n));
    std::vector<double> x_ref(n);
    for (int j = 0; j<n; ++j) {
        for (int i = 0; i<n; ++i) {
            A[j][i] = j+n*i+1;
        }
        x_ref[j] = j+1;
    }
    //
    std::vector<double> b(n);
    for (int i = 0; i<n; ++i) {
        b[i] = A[i][i]*x_ref[i];
    }
    std::vector<std::complex<double> > b_c(n);
    for (int i = 0; i<n; ++i) {
        b_c[i] = b[i];
    }
    //
    std::vector<double> x(n);
    eptlib::linalg::SolveDiag(x.data(),A,b.data(),n);
    std::vector<std::complex<double> > x_c(n);
    eptlib::linalg::SolveDiag(x_c.data(),A,b_c.data(),n);
    //
    for (int i = 0; i<n; ++i) {
        ASSERT_NEAR(x[i],x_ref[i],1e-12);
        ASSERT_NEAR(x_c[i].real(),x_ref[i],1e-12);
        ASSERT_NEAR(x_c[i].imag(),0.0,1e-12);
    }
}

TEST(LinalgUtilGTest,SolveTriU) {
    const int n = 10;
    eptlib::linalg::MatrixReal A(n,std::vector<double>(n));
    std::vector<double> x_ref(n);
    for (int j = 0; j<n; ++j) {
        for (int i = 0; i<n; ++i) {
            A[j][i] = j+n*i+1;
        }
        x_ref[j] = j+1;
    }
    //
    std::vector<double> b(n);
    for (int i = 0; i<n; ++i) {
        b[i] = 0.0;
        for (int j = i; j<n; ++j) {
            b[i] += A[j][i]*x_ref[j];
        }
    }
    std::vector<std::complex<double> > b_c(n);
    for (int i = 0; i<n; ++i) {
        b_c[i] = b[i];
    }
    //
    std::vector<double> x(n);
    eptlib::linalg::SolveTriU(x.data(),A,b.data(),n);
    std::vector<std::complex<double> > x_c(n);
    eptlib::linalg::SolveTriU(x_c.data(),A,b_c.data(),n);
    //
    for (int i = 0; i<n; ++i) {
        ASSERT_NEAR(x[i],x_ref[i],1e-12);
        ASSERT_NEAR(x_c[i].real(),x_ref[i],1e-12);
        ASSERT_NEAR(x_c[i].imag(),0.0,1e-12);
    }
}

TEST(LinalgUtilGTest,SolveTriL) {
    const int n = 10;
    eptlib::linalg::MatrixReal A(n,std::vector<double>(n));
    std::vector<double> x_ref(n);
    for (int j = 0; j<n; ++j) {
        for (int i = 0; i<n; ++i) {
            A[j][i] = j+n*i+1;
        }
        x_ref[j] = j+1;
    }
    //
    std::vector<double> b(n);
    for (int i = 0; i<n; ++i) {
        b[i] = 0.0;
        for (int j = 0; j<=i; ++j) {
            b[i] += A[j][i]*x_ref[j];
        }
    }
    std::vector<std::complex<double> > b_c(n);
    for (int i = 0; i<n; ++i) {
        b_c[i] = b[i];
    }
    //
    std::vector<double> x(n);
    eptlib::linalg::SolveTriL(x.data(),A,b.data(),n);
    std::vector<std::complex<double> > x_c(n);
    eptlib::linalg::SolveTriL(x_c.data(),A,b_c.data(),n);
    //
    for (int i = 0; i<n; ++i) {
        ASSERT_NEAR(x[i],x_ref[i],1e-12);
        ASSERT_NEAR(x_c[i].real(),x_ref[i],1e-12);
        ASSERT_NEAR(x_c[i].imag(),0.0,1e-12);
    }
}
