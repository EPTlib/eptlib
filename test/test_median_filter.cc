/*****************************************************************************
*
*     Program: EPTlib
*     Author: Alessandro Arduino <a.arduino@inrim.it>
*
*  MIT License
*
*  Copyright (c) 2020-2022  Alessandro Arduino
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

#include "eptlib/median_filter.h"

using namespace eptlib;

TEST(MedianFilterGTest,ApplyFilter) {
    const std::array<int,NDIM> nn = {5,5,5};
    const std::array<int,NDIM> mm = {3,3,3};
    //
    std::array<int,NDIM> rr;
    for (int d = 0; d<NDIM; ++d) {
        rr[d] = mm[d]/2;
    }
    Shape shape = shapes::Cuboid(mm);
    MedianFilter median(shape);
    //
    const int n_vox = std::accumulate(nn.begin(),nn.end(),1,std::multiplies<int>());
    std::vector<double> field1(n_vox,2.0);
    for (int i2 = 0; i2 < nn[2]; ++i2) {
        for (int i1 = 0; i1<nn[1]; ++i1) {
            int i0 = nn[0]/2;
            int idx = i0 + nn[0]*(i1 + nn[1]*i2);
            field1[idx] = 1.0;
        }
    }
    for (int i2 = 0; i2 < nn[2]; ++i2) {
        for (int i0 = 0; i0<nn[0]; ++i0) {
            int i1 = nn[1]/2;
            int idx = i0 + nn[0]*(i1 + nn[1]*i2);
            field1[idx] = 1.0;
        }
    }
    for (int i1 = 0; i1 < nn[1]; ++i1) {
        for (int i0 = 0; i0<nn[0]; ++i0) {
            int i2 = nn[2]/2;
            int idx = i0 + nn[0]*(i1 + nn[1]*i2);
            field1[idx] = 1.0;
        }
    }
    std::vector<double> post1(n_vox);
    median.ApplyFilter(post1.data(),field1.data(),nn);
    for (int i2 = rr[2]; i2<nn[2]-rr[2]; ++i2) {
        for (int i1 = rr[1]; i1<nn[1]-rr[1]; ++i1) {
            for (int i0 = rr[0]; i0<nn[0]-rr[0]; ++i0) {
                int idx = i0 + nn[0]*(i1 + nn[1]*i2);
                ASSERT_DOUBLE_EQ(post1[idx],1.0);
            }
        }
    }
    //
    std::vector<double> field2(n_vox,2.0);
    for (int i2 = 0; i2 < nn[2]; ++i2) {
        int i0 = nn[0]/2;
        int i1 = nn[1]/2;
        int idx = i0 + nn[0]*(i1 + nn[1]*i2);
        field1[idx] = 1.0;
    }
    for (int i1 = 0; i1 < nn[1]; ++i1) {
        int i0 = nn[0]/2;
        int i2 = nn[2]/2;
        int idx = i0 + nn[0]*(i1 + nn[1]*i2);
        field1[idx] = 1.0;
    }
    for (int i0 = 0; i0<nn[0]; ++i0) {
        int i1 = nn[1]/2;
        int i2 = nn[2]/2;
        int idx = i0 + nn[0]*(i1 + nn[1]*i2);
        field1[idx] = 1.0;
    }
    std::vector<double> post2(n_vox);
    median.ApplyFilter(post2.data(),field2.data(),nn);
    for (int i2 = rr[2]; i2<nn[2]-rr[2]; ++i2) {
        for (int i1 = rr[1]; i1<nn[1]-rr[1]; ++i1) {
            for (int i0 = rr[0]; i0<nn[0]-rr[0]; ++i0) {
                int idx = i0 + nn[0]*(i1 + nn[1]*i2);
                ASSERT_DOUBLE_EQ(post2[idx],2.0);
            }
        }
    }
}
