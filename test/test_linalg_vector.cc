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

#include "eptlib/linalg/vector.h"

#include <numeric>
#include <vector>

TEST(LinalgVectorGTest,Norm2) {
    const int n = 10;
    std::vector<double> x(n);
    std::iota(x.begin(), x.end(), 1.0);
    double norm = eptlib::linalg::Norm2(x.begin(), x.end());
    ASSERT_NEAR(norm, std::sqrt(385.0), 1e-12);
}

TEST(LinalgVectorGTest,MaxAbs) {
    const int n = 10;
    std::vector<double> x(n);
    for (int i = 0; i<n; ++i) {
        x[i] = -(i+1);
    }
    double maxabs = eptlib::linalg::MaxAbs(x.begin(), x.end());
    ASSERT_NEAR(maxabs, 10.0, 1e-12);
}

TEST(LinalgVectorGTest,Permute) {
    const int n = 10;
    std::vector<double> x(n);
    std::vector<size_t> p(n);
    for (int i = 0; i<n; ++i) {
        x[i] = -(i+1);
        p[i] = n-1-i;
    }
    eptlib::linalg::Permute(x.begin(), x.end(), p);
    for (int i = 0; i<n; ++i) {
        ASSERT_EQ(x[i], -n+i);
    }
}
