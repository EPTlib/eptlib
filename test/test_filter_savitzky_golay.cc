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
    const std::array<int,N_DIM> nn = {7,5,3};
    const std::array<double,N_DIM> dd = {1.0,2.0,3.0};
    std::array<size_t,N_DIM> ii;
    const int n_vox = std::accumulate(nn.begin(),nn.end(),1,std::multiplies<int>());
    //
    std::array<int,N_DIM> rr = {1,1,1};
    std::array<int,N_DIM> rr2 = {1,0,0};
    Shape cube = shapes::CuboidR(rr[0],rr[1],rr[2]);
    SavitzkyGolay fd_lapl(dd[0],dd[1],dd[2], cube, 2);
    //
    Image<double> q_field(nn[0],nn[1],nn[2]);
    for (int idx = 0; idx<n_vox; ++idx) {
        IdxToIJK(ii[0],ii[1],ii[2],idx,nn[0],nn[1]);
        q_field(idx) = 0.0;
        for (int d = 0; d<N_DIM; ++d) {
            q_field(idx) += ii[d]*ii[d]*dd[d]*dd[d];
        }
    }
    //
    Image<double> q_lapl(nn[0],nn[1],nn[2]);
    EPTlibError error = fd_lapl.Apply(&q_lapl, q_field);
    ASSERT_EQ(error, EPTlibError::Success);
    //
    for (int idx = 0; idx<n_vox; ++idx) {
        IdxToIJK(ii[0],ii[1],ii[2],idx,nn[0],nn[1]);
        bool kernel_in_domain = true;
        for (int d = 0; d<N_DIM; ++d) {
            if (ii[d]==0 || ii[d]==nn[d]-1) {
                kernel_in_domain = false;
                break;
            }
        }
        if (kernel_in_domain) {
            ASSERT_NEAR(q_lapl(idx),2.0*N_DIM,1e-12);
        } else {
            ASSERT_NEAR(q_lapl(idx),0.0,1e-12);
        }
    }
}
