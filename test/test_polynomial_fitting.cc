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

#include "eptlib/polynomial/fitting.h"

TEST(PolynomialFittingGTest,GetNumberOfMonomials) {
    ASSERT_EQ(eptlib::polynomial::GetNumberOfMonomials(0), 1);
    ASSERT_EQ(eptlib::polynomial::GetNumberOfMonomials(1), 4);
    ASSERT_EQ(eptlib::polynomial::GetNumberOfMonomials(2), 10);
    ASSERT_EQ(eptlib::polynomial::GetNumberOfMonomials(3), 20);
}

TEST(PolynomialFittingGTest,EvaluateMonomials) {
    double x = 0.1;
    double y = 0.2;
    double z = 0.3;
    std::vector<double> p0 = eptlib::polynomial::EvaluateMonomials(x,y,z, 0);
    std::vector<double> p1 = eptlib::polynomial::EvaluateMonomials(x,y,z, 1);
    std::vector<double> p2 = eptlib::polynomial::EvaluateMonomials(x,y,z, 2);
    std::vector<double> p3 = eptlib::polynomial::EvaluateMonomials(x,y,z, 3);
    ASSERT_DOUBLE_EQ(p0[0], 1.0);
    ASSERT_DOUBLE_EQ(p1[0], x);
    ASSERT_DOUBLE_EQ(p1[1], y);
    ASSERT_DOUBLE_EQ(p1[2], z);
    ASSERT_DOUBLE_EQ(p2[0], x*x);
    ASSERT_DOUBLE_EQ(p2[1], x*y);
    ASSERT_DOUBLE_EQ(p2[2], y*y);
    ASSERT_DOUBLE_EQ(p2[3], x*z);
    ASSERT_DOUBLE_EQ(p2[4], y*z);
    ASSERT_DOUBLE_EQ(p2[5], z*z);
    ASSERT_DOUBLE_EQ(p3[0], x*x*x);
    ASSERT_DOUBLE_EQ(p3[1], x*x*y);
    ASSERT_DOUBLE_EQ(p3[2], x*y*y);
    ASSERT_DOUBLE_EQ(p3[3], y*y*y);
    ASSERT_DOUBLE_EQ(p3[4], x*x*z);
    ASSERT_DOUBLE_EQ(p3[5], x*y*z);
    ASSERT_DOUBLE_EQ(p3[6], y*y*z);
    ASSERT_DOUBLE_EQ(p3[7], x*z*z);
    ASSERT_DOUBLE_EQ(p3[8], y*z*z);
    ASSERT_DOUBLE_EQ(p3[9], z*z*z);
}

TEST(PolynomialFittingGTest,DesignMatrixWithMonomialBasis) {
    std::vector<double> x = {0.1, 0.2, 0.3};
    std::vector<double> y = {0.1, 0.2, 0.3};
    std::vector<double> z = {0.1, 0.2, 0.3};
    eptlib::linalg::Matrix<double> A = eptlib::polynomial::DesignMatrixWithMonomialBasis(x,y,z, 2);
    ASSERT_DOUBLE_EQ(A(0,0), 1.0);       ASSERT_DOUBLE_EQ(A(1,0), 1.0);       ASSERT_DOUBLE_EQ(A(1,0), 1.0);
    ASSERT_DOUBLE_EQ(A(0,1), x[0]);      ASSERT_DOUBLE_EQ(A(1,1), x[1]);      ASSERT_DOUBLE_EQ(A(2,1), x[2]);
    ASSERT_DOUBLE_EQ(A(0,2), y[0]);      ASSERT_DOUBLE_EQ(A(1,2), y[1]);      ASSERT_DOUBLE_EQ(A(2,2), y[2]);
    ASSERT_DOUBLE_EQ(A(0,3), z[0]);      ASSERT_DOUBLE_EQ(A(1,3), z[1]);      ASSERT_DOUBLE_EQ(A(2,3), z[2]);
    ASSERT_DOUBLE_EQ(A(0,4), x[0]*x[0]); ASSERT_DOUBLE_EQ(A(1,4), x[1]*x[1]); ASSERT_DOUBLE_EQ(A(2,4), x[2]*x[2]);
    ASSERT_DOUBLE_EQ(A(0,5), x[0]*y[0]); ASSERT_DOUBLE_EQ(A(1,5), x[1]*y[1]); ASSERT_DOUBLE_EQ(A(2,5), x[2]*y[2]);
    ASSERT_DOUBLE_EQ(A(0,6), y[0]*y[0]); ASSERT_DOUBLE_EQ(A(1,6), y[1]*y[1]); ASSERT_DOUBLE_EQ(A(2,6), y[2]*y[2]);
    ASSERT_DOUBLE_EQ(A(0,7), x[0]*z[0]); ASSERT_DOUBLE_EQ(A(1,7), x[1]*z[1]); ASSERT_DOUBLE_EQ(A(2,7), x[2]*z[2]);
    ASSERT_DOUBLE_EQ(A(0,8), y[0]*z[0]); ASSERT_DOUBLE_EQ(A(1,8), y[1]*z[1]); ASSERT_DOUBLE_EQ(A(2,8), y[2]*z[2]);
    ASSERT_DOUBLE_EQ(A(0,9), z[0]*z[0]); ASSERT_DOUBLE_EQ(A(1,9), z[1]*z[1]); ASSERT_DOUBLE_EQ(A(2,9), z[2]*z[2]);
}

TEST(PolynomialFittingGTest,Fitting) {
    std::vector<double> x = {0, 0, -1, 0, 1, 0, 0};
    std::vector<double> y = {0, -1, 0, 0, 0, 1, 0};
    std::vector<double> z = {-1, 0, 0, 0, 0, 0, 1};
    std::vector<double> values(x.size());
    for (size_t idx = 0; idx < values.size(); ++idx) {
        values[idx] = 1.0 + 2.0*x[idx] + 3.0*y[idx] + 4.0*z[idx] + 5.0*x[idx]*x[idx] + 6.0*y[idx]*y[idx] + 7.0*z[idx]*z[idx];
    }
    auto [param, chi2n] = eptlib::polynomial::Fitting(x,y,z, 2, values);
    ASSERT_NEAR(param[0], 1.0, 1e-14);
    ASSERT_NEAR(param[1], 2.0, 1e-14);
    ASSERT_NEAR(param[2], 3.0, 1e-14);
    ASSERT_NEAR(param[3], 4.0, 1e-14);
    ASSERT_NEAR(param[4], 5.0, 1e-14);
    ASSERT_DOUBLE_EQ(param[5], 0.0);
    ASSERT_NEAR(param[6], 6.0, 1e-14);
    ASSERT_DOUBLE_EQ(param[7], 0.0);
    ASSERT_DOUBLE_EQ(param[8], 0.0);
    ASSERT_NEAR(param[9], 7.0, 1e-14);
    ASSERT_DOUBLE_EQ(chi2n, 0.0);
}

TEST(PolynomialFittingGTest,Fitting2D) {
    std::vector<double> x = {0, -1, 0, 1, 0};
    std::vector<double> y = {-1, 0, 0, 0, 1};
    std::vector<double> z = {0, 0, 0, 0, 0};
    std::vector<double> values(x.size());
    for (size_t idx = 0; idx < values.size(); ++idx) {
        values[idx] = 1.0 + 2.0*x[idx] + 3.0*y[idx] + 4.0*x[idx]*x[idx] + 5.0*y[idx]*y[idx];
    }
    auto [param, chi2n] = eptlib::polynomial::Fitting(x,y,z, 2, values);
    ASSERT_NEAR(param[0], 1.0, 1e-14);
    ASSERT_NEAR(param[1], 2.0, 1e-14);
    ASSERT_NEAR(param[2], 3.0, 1e-14);
    ASSERT_DOUBLE_EQ(param[3], 0.0);
    ASSERT_NEAR(param[4], 4.0, 1e-14);
    ASSERT_DOUBLE_EQ(param[5], 0.0);
    ASSERT_NEAR(param[6], 5.0, 1e-14);
    ASSERT_DOUBLE_EQ(param[7], 0.0);
    ASSERT_DOUBLE_EQ(param[8], 0.0);
    ASSERT_DOUBLE_EQ(param[9], 0.0);
    ASSERT_DOUBLE_EQ(chi2n, 0.0);
}
