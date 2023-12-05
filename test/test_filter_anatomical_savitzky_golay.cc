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

#include "eptlib/filter/anatomical_savitzky_golay.h"

TEST(FilterAnatomicalSavitzkyGolayGTest,AnatomicalSavitzkyGolayApply) {
    const size_t n0 = 10;
    const size_t n1 = 10;
    const size_t n2 = 10;

    const double d0 = 1.0;
    const double d1 = 1.0;
    const double d2 = 1.0;

    eptlib::Image<double> constant_field (n0,n1,n2);
    eptlib::Image<double> linear_field   (n0,n1,n2);
    eptlib::Image<double> quadratic_field(n0,n1,n2);

    for (int k = 0; k<n2; ++k) {
        for (int j = 0; j<n1; ++j) {
            for (int i = 0; i<n0; ++i) {
                double x = i*d0;
                double y = j*d1;
                double z = k*d2;

                constant_field (i,j,k) = 1.0;
                linear_field   (i,j,k) = x + y + z;
                quadratic_field(i,j,k) = x*x + y*y + z*z;
            }
        }
    }

    eptlib::Image<double> ref_img(n0,n1,n2);
    for (int k = 0; k<n2; ++k) {
        for (int j = 0; j<n1; ++j) {
            for (int i = 0; i<n0; ++i) {
                if (i < n0/2) {
                    ref_img(i,j,k) = 1.0;
                } else {
                    ref_img(i,j,k) = 2.0;
                }
            }
        }
    }

    const eptlib::Shape window = eptlib::shapes::Cross(3,2,2);
    const size_t degree = 2;
    eptlib::filter::AnatomicalSavitzkyGolay asg_filter(d0,d1,d2, window, degree);

    eptlib::Image<double> lapl_constant_field (n0,n1,n2);
    eptlib::Image<double> lapl_linear_field   (n0,n1,n2);
    eptlib::Image<double> lapl_quadratic_field(n0,n1,n2);

    eptlib::EPTlibError error;

    error = asg_filter.Apply(eptlib::DifferentialOperator::Laplacian, &lapl_constant_field, constant_field, ref_img);
    ASSERT_EQ(error, eptlib::EPTlibError::Success);

    error = asg_filter.Apply(eptlib::DifferentialOperator::Laplacian, &lapl_linear_field, linear_field, ref_img);
    ASSERT_EQ(error, eptlib::EPTlibError::Success);

    error = asg_filter.Apply(eptlib::DifferentialOperator::Laplacian, &lapl_quadratic_field, quadratic_field, ref_img);
    ASSERT_EQ(error, eptlib::EPTlibError::Success);

    for (int k = 0; k<n2; ++k) {
        for (int j = 0; j<n1; ++j) {
            for (int i = 0; i<n0; ++i) {
                ASSERT_NEAR(lapl_constant_field (i,j,k), 0.0, 1e-13);
                ASSERT_NEAR(lapl_linear_field   (i,j,k), 0.0, 1e-13);
                ASSERT_NEAR(lapl_quadratic_field(i,j,k), 6.0, 1e-12);
            }
        }
    }

    eptlib::Image<double> variance(n0,n1,n2);

    error = asg_filter.Apply(eptlib::DifferentialOperator::Laplacian, &lapl_constant_field, &variance, constant_field, ref_img);
    ASSERT_EQ(error, eptlib::EPTlibError::Success);

    error = asg_filter.Apply(eptlib::DifferentialOperator::Laplacian, &lapl_linear_field, &variance, linear_field, ref_img);
    ASSERT_EQ(error, eptlib::EPTlibError::Success);

    error = asg_filter.Apply(eptlib::DifferentialOperator::Laplacian, &lapl_quadratic_field, &variance, quadratic_field, ref_img);
    ASSERT_EQ(error, eptlib::EPTlibError::Success);

    for (int k = 0; k<n2; ++k) {
        for (int j = 0; j<n1; ++j) {
            for (int i = 0; i<n0; ++i) {
                ASSERT_NEAR(lapl_constant_field (i,j,k), 0.0, 1e-13);
                ASSERT_NEAR(lapl_linear_field   (i,j,k), 0.0, 1e-13);
                ASSERT_NEAR(lapl_quadratic_field(i,j,k), 6.0, 1e-12);
            }
        }
    }

    return;
}
