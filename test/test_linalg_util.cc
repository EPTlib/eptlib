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

#include "eptlib/linalg/linalg_util.h"

#include <cmath>
#include <vector>

#include "eptlib/util.h"

using namespace eptlib;
using namespace eptlib::linalg;

TEST(LinalgUtilGTest,Norm2) {
    const int n = 10;
    std::vector<real_t> x(n);
    for (int i = 0; i<n; ++i) {
        x[i] = i+1;
    }
    //
    real_t norm = Norm2(x.data(),n);
    //
    real_t ref = std::sqrt(385.0);
    ASSERT_NEAR(norm,ref,1e-12);
}
TEST(LinalgUtilGTest,Dot) {
    const int n = 10;
    std::vector<real_t> x(n);
    for (int i = 0; i<n; ++i) {
        x[i] = i+1;
    }
    //
    real_t dot = Dot(x.data(),x.data(),n);
    //
    real_t ref = 385.0;
    ASSERT_NEAR(dot,ref,1e-12);
}
TEST(LinalgUtilGTest,MaxAbs) {
    const int n = 10;
    std::vector<real_t> x(n);
    for (int i = 0; i<n; ++i) {
        x[i] = -(i+1);
    }
    //
    real_t maxabs = MaxAbs(x.data(),n);
    //
    real_t ref = 10.0;
    ASSERT_NEAR(maxabs,ref,1e-12);
}
TEST(LinalgUtilGTest,SolveDiag) {
    const int n = 10;
    MatrixReal A(n,std::vector<real_t>(n));
    std::vector<real_t> x_ref(n);
    for (int j = 0; j<n; ++j) {
        for (int i = 0; i<n; ++i) {
            A[j][i] = j+n*i+1;
        }
        x_ref[j] = j+1;
    }
    //
    std::vector<real_t> b(n);
    for (int i = 0; i<n; ++i) {
        b[i] = A[i][i]*x_ref[i];
    }
    //
    std::vector<real_t> x(n);
    SolveDiag(x.data(),A,b.data(),n);
    //
    for (int i = 0; i<n; ++i) {
        ASSERT_NEAR(x[i],x_ref[i],1e-12);
    }
}
TEST(LinalgUtilGTest,SolveTriU) {
    const int n = 10;
    MatrixReal A(n,std::vector<real_t>(n));
    std::vector<real_t> x_ref(n);
    for (int j = 0; j<n; ++j) {
        for (int i = 0; i<n; ++i) {
            A[j][i] = j+n*i+1;
        }
        x_ref[j] = j+1;
    }
    //
    std::vector<real_t> b(n);
    for (int i = 0; i<n; ++i) {
        b[i] = 0.0;
        for (int j = i; j<n; ++j) {
            b[i] += A[j][i]*x_ref[j];
        }
    }
    //
    std::vector<real_t> x(n);
    SolveTriU(x.data(),A,b.data(),n);
    //
    for (int i = 0; i<n; ++i) {
        ASSERT_NEAR(x[i],x_ref[i],1e-12);
    }
}
TEST(LinalgUtilGTest,SolveTriL) {
    const int n = 10;
    MatrixReal A(n,std::vector<real_t>(n));
    std::vector<real_t> x_ref(n);
    for (int j = 0; j<n; ++j) {
        for (int i = 0; i<n; ++i) {
            A[j][i] = j+n*i+1;
        }
        x_ref[j] = j+1;
    }
    //
    std::vector<real_t> b(n);
    for (int i = 0; i<n; ++i) {
        b[i] = 0.0;
        for (int j = 0; j<=i; ++j) {
            b[i] += A[j][i]*x_ref[j];
        }
    }
    //
    std::vector<real_t> x(n);
    SolveTriL(x.data(),A,b.data(),n);
    //
    for (int i = 0; i<n; ++i) {
        ASSERT_NEAR(x[i],x_ref[i],1e-12);
    }
}
