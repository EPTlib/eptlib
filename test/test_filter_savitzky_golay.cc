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

#include "eptlib/filter/savitzky_golay.h"

#include <algorithm>
#include <array>
#include <functional>
#include <numeric>
#include <typeinfo>
#include <vector>

#include "eptlib/image.h"
#include "eptlib/util.h"

using namespace eptlib;

TEST(FilterSavitzkyGolayGTest,SavitzkyGolayLaplacian) {
    const size_t n0 = 10;
    const size_t n1 = 10;
    const size_t n2 = 10;

    const double d0 = 1.0;
    const double d1 = 1.0;
    const double d2 = 1.0;

    Image<double> constant_field (n0,n1,n2);
    Image<double> linear_field   (n0,n1,n2);
    Image<double> quadratic_field(n0,n1,n2);

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

    const Shape window = shapes::CuboidR(1,1,1);
    const size_t degree = 2;
    SavitzkyGolay sg_filter(d0,d1,d2, window, degree);

    Image<double> lapl_constant_field (n0,n1,n2);
    Image<double> lapl_linear_field   (n0,n1,n2);
    Image<double> lapl_quadratic_field(n0,n1,n2);

    EPTlibError error;

    error = sg_filter.Apply(&lapl_constant_field, constant_field);
    ASSERT_EQ(error, EPTlibError::Success);

    error = sg_filter.Apply(&lapl_linear_field, linear_field);
    ASSERT_EQ(error, EPTlibError::Success);

    error = sg_filter.Apply(&lapl_quadratic_field, quadratic_field);
    ASSERT_EQ(error, EPTlibError::Success);

    for (int k = 0; k<n2; ++k) {
        for (int j = 0; j<n1; ++j) {
            for (int i = 0; i<n0; ++i) {
                ASSERT_NEAR(lapl_constant_field(i,j,k), 0.0, 1e-12);
                ASSERT_NEAR(lapl_linear_field  (i,j,k), 0.0, 1e-12);
                if (i==0 || i==n0-1 || j==0 || j==n1-1 || k==0 || k==n2-1) {
                    ASSERT_NEAR(lapl_quadratic_field(i,j,k), 0.0, 1e-12);
                } else {
                    ASSERT_NEAR(lapl_quadratic_field(i,j,k), 6.0, 1e-12);
                }
            }
        }
    }

}
