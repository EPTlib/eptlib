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

#include "eptlib/filter/moving_window.h"

#include <vector>

#include "eptlib/image.h"
#include "eptlib/shape.h"
#include "eptlib/util.h"

TEST(FilterMovingWindowGTest,MovingWindow) {
    const size_t n0 = 10;
    const size_t n1 = 10;
    const size_t n2 = 10;
    const size_t r0 = 1;
    const size_t r1 = 2;
    const size_t r2 = 3;
    eptlib::Image<double> img_in (n0, n1, n2);
    eptlib::Image<double> img_out(n0, n1, n2);
    eptlib::Shape window = eptlib::shapes::Ellipsoid(r0, r1, r2);
    for (int idx = 0; idx<img_in.GetNVox(); ++idx) {
        img_in (idx) = 1.0;
        img_out(idx) = 0.0;
    }
    auto filter = [](const std::vector<double> &crop_in) -> double {
        return eptlib::Sum(crop_in);
    };
    eptlib::EPTlibError error = eptlib::filter::MovingWindow(&img_out, img_in, window, filter);
    ASSERT_EQ(error, eptlib::EPTlibError::Success);
    for (size_t i2 = 0; i2<n2; ++i2) {
        for (size_t i1 = 0; i1<n1; ++i1) {
            for (size_t i0 = 0; i0<n0; ++i0) {
                if (i0<r0 || i0>=n0-r0 || i1<r1 || i1>=n1-r1 || i2<r2 || i2>=n2-r2) {
                    ASSERT_TRUE(std::isnan(img_out(i0,i1,i2)));
                } else {
                    ASSERT_EQ(img_out(i0,i1,i2), window.GetVolume());
                }
            }
        }
    }
}
