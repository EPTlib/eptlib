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

#include "eptlib/util.h"

#include <complex>
#include <vector>

TEST(UtilGTest,Sum) {
    std::vector<int> vec_of_int(100);
    std::iota(vec_of_int.begin(), vec_of_int.end(),1);
    ASSERT_EQ(eptlib::Sum(vec_of_int), 5050);
    std::vector<double> vec_of_double(100);
    std::copy(vec_of_int.begin(), vec_of_int.end(), vec_of_double.begin());
    ASSERT_DOUBLE_EQ(eptlib::Sum(vec_of_double), 5050.0);
}

TEST(UtilGTest,Prod) {
    std::vector<int> vec_of_int(5);
    std::iota(vec_of_int.begin(), vec_of_int.end(), 1);
    ASSERT_EQ(eptlib::Prod(vec_of_int), 120);
    std::vector<double> vec_of_double(5);
    std::copy(vec_of_int.begin(), vec_of_int.end(), vec_of_double.begin());
    ASSERT_DOUBLE_EQ(eptlib::Prod(vec_of_double), 120.0);
}

TEST(UtilGTest,ArithmeticMean) {
    std::vector<double> vec_of_double(100);
    std::iota(vec_of_double.begin(), vec_of_double.end(), 1);
    ASSERT_DOUBLE_EQ(eptlib::ArithmeticMean(vec_of_double), 50.5);
}

TEST(UtilGTest,MaxMin) {
    std::vector<double> vec_of_double{3.0,1.0,-10.0,7.0,-0.5};
    ASSERT_DOUBLE_EQ(eptlib::Max(vec_of_double), 7.0);
    ASSERT_DOUBLE_EQ(eptlib::Min(vec_of_double), -10.0);
    ASSERT_DOUBLE_EQ(eptlib::MaxAbs(vec_of_double), 10.0);
    ASSERT_DOUBLE_EQ(eptlib::MinAbs(vec_of_double), 0.5);
    std::vector<std::complex<double> > vec_of_complex{std::complex(3.0,0.0), std::complex(-10.0,0.0), std::complex(0.0,15.0)};
    ASSERT_DOUBLE_EQ(eptlib::MaxAbs(vec_of_complex), 15.0);
    ASSERT_DOUBLE_EQ(eptlib::MinAbs(vec_of_complex), 3.0);
}

TEST(UtilGTest,IJKToIdx) {
    size_t i = 20;
    size_t j = 5;
    size_t k = 15;
    size_t n0 = 100;
    size_t n1 = 50;
    ASSERT_EQ(eptlib::IJKToIdx(i,j,k,n0,n1), 75520);
}

TEST(UtilGTest,IdxToIJK) {
    size_t i,j,k;
    size_t n0 = 100;
    size_t n1 = 50;
    size_t idx = 75520;
    eptlib::IdxToIJK(i,j,k, idx,n0,n1);
    ASSERT_EQ(i,20);
    ASSERT_EQ(j,5);
    ASSERT_EQ(k,15);
}
