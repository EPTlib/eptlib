/*****************************************************************************
*
*     Program: EPTlib
*     Author: Alessandro Arduino <a.arduino@inrim.it>
*
*  MIT License
*
*  Copyright (c) 2020  Alessandro Arduino
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

#include "eptlib/finite_difference.h"

#include <algorithm>
#include <array>
#include <functional>
#include <numeric>
#include <typeinfo>
#include <vector>

using namespace eptlib;

TEST(FiniteDifferenceGTest,FDLaplacianKernel) {
    const std::array<int,NDIM> nn = {7,5,3};
    const std::array<double,NDIM> dd = {1.0,2.0,3.0};
    std::array<int,NDIM> ii;
    const int n_vox = std::accumulate(nn.begin(),nn.end(),1,std::multiplies<int>());
    //
    std::array<int,NDIM> rr = {1,1,1};
    Shape cross = shapes::Cross(rr);
    ASSERT_TRUE(cross.IsSymmetric());
    FDLaplacianKernel fd_lapl(cross);
    //
    std::vector<double> c_field(n_vox);
    std::vector<double> l_field(n_vox);
    std::vector<double> q_field(n_vox);
    for (int idx = 0; idx<n_vox; ++idx) {
        IdxToMultiIdx(ii,idx,nn);
        c_field[idx] = 1.0;
        l_field[idx] = std::accumulate(ii.begin(),ii.end(),0.0);
        q_field[idx] = 0.0;
        for (int d = 0; d<NDIM; ++d) {
            q_field[idx] += ii[d]*ii[d]*dd[d]*dd[d];
        }
    }
    //
    std::vector<double> c_lapl(n_vox);
    std::vector<double> l_lapl(n_vox);
    std::vector<double> q_lapl(n_vox);
    fd_lapl.ComputeLaplacian(c_lapl.data(),c_field.data(),nn,dd);
    fd_lapl.ComputeLaplacian(l_lapl.data(),l_field.data(),nn,dd);
    fd_lapl.ComputeLaplacian(q_lapl.data(),q_field.data(),nn,dd);
    //
    for (int idx = 0; idx<n_vox; ++idx) {
        ASSERT_NEAR(c_lapl[idx],0.0,1e-14);
        ASSERT_NEAR(l_lapl[idx],0.0,1e-13);
        IdxToMultiIdx(ii,idx,nn);
        bool kernel_in_domain = true;
        for (int d = 0; d<NDIM; ++d) {
            if (ii[d]==0 || ii[d]==nn[d]-1) {
                kernel_in_domain = false;
                break;
            }
        }
        if (kernel_in_domain) {
            ASSERT_NEAR(q_lapl[idx],2.0*NDIM,1e-12);
        } else {
            ASSERT_NEAR(q_lapl[idx],0.0,1e-12);
        }
    }
}

TEST(FiniteDifferenceGTest,AsymmetricFDLaplacianKernel) {
    const std::array<int,NDIM> nn = {7,5,3};
    const std::array<double,NDIM> dd = {1.0,2.0,3.0};
    std::array<int,NDIM> ii;
    const int n_vox = std::accumulate(nn.begin(),nn.end(),1,std::multiplies<int>());
    //
    const std::array<int,NDIM> cross_rr = {2,1,1};
    std::array<int,NDIM> cuboid_nn = {4,2,2};
    Shape cuboid = shapes::Cuboid(cuboid_nn);
    std::array<int,NDIM> l = {1,1,1};
    std::array<int,NDIM> r = {0,0,0};
    cuboid.Pad(l,r);
    Shape cross = cuboid;
    cuboid_nn = {1,3,1};
    cuboid = shapes::Cuboid(cuboid_nn);
    l = {2,0,1};
    r = {2,0,1};
    cuboid.Pad(l,r);
    cross += cuboid;
    cuboid_nn = {1,1,3};
    cuboid = shapes::Cuboid(cuboid_nn);
    l = {2,1,0};
    r = {2,1,0};
    cuboid.Pad(l,r);
    cross += cuboid;
    ASSERT_FALSE(cross.IsSymmetric());
    FDLaplacianKernel fd_lapl(cross);
    //
    std::vector<double> c_field(n_vox);
    std::vector<double> l_field(n_vox);
    std::vector<double> q_field(n_vox);
    for (int idx = 0; idx<n_vox; ++idx) {
        IdxToMultiIdx(ii,idx,nn);
        c_field[idx] = 1.0;
        l_field[idx] = std::accumulate(ii.begin(),ii.end(),0.0);
        q_field[idx] = 0.0;
        for (int d = 0; d<NDIM; ++d) {
            q_field[idx] += ii[d]*ii[d]*dd[d]*dd[d];
        }
    }
    //
    std::vector<double> c_lapl(n_vox);
    std::vector<double> l_lapl(n_vox);
    std::vector<double> q_lapl(n_vox);
    fd_lapl.ComputeLaplacian(c_lapl.data(),c_field.data(),nn,dd);
    fd_lapl.ComputeLaplacian(l_lapl.data(),l_field.data(),nn,dd);
    fd_lapl.ComputeLaplacian(q_lapl.data(),q_field.data(),nn,dd);
    //
    for (int idx = 0; idx<n_vox; ++idx) {
        ASSERT_NEAR(c_lapl[idx],0.0,1e-14);
        ASSERT_NEAR(l_lapl[idx],0.0,1e-13);
        IdxToMultiIdx(ii,idx,nn);
        bool kernel_in_domain = true;
        if (kernel_in_domain) {
            for (int d = 0; d<NDIM; ++d) {
                if (ii[d]<cross_rr[d] || ii[d]>nn[d]-1-cross_rr[d]) {
                    kernel_in_domain = false;
                    break;
                }
            }
        }
        if (kernel_in_domain) {
            ASSERT_NEAR(q_lapl[idx],2.0*NDIM,1e-12);
        } else {
            ASSERT_NEAR(q_lapl[idx],0.0,1e-12);
        }
    }
}
