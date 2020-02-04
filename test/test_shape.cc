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

#include "eptlib/shape.h"

#include <algorithm>
#include <array>
#include <functional>
#include <numeric>
#include <vector>

using namespace eptlib;

TEST(ShapeGTest,Cuboid) {
    constexpr int n_dim = 3;
    const std::array<int,n_dim> nn = {7,3,3};
    std::array<int,n_dim> ii;
    const int n_vox = std::accumulate(nn.begin(),nn.end(),1,std::multiplies<int>());
    int idx;
    //
    Shape cuboid = shapes::Cuboid(nn);
    idx = 0;
    for (int i2 = 0; i2<nn[2]; ++i2) {
        for (int i1 = 0; i1<nn[1]; ++i1) {
            for (int i0 = 0; i0<nn[0]; ++i0) {
                ii = {i0,i1,i2};
                ASSERT_TRUE(cuboid[ii]);
                ASSERT_TRUE(cuboid[idx]);
                ++idx;
            }
        }
    }
    ASSERT_TRUE(cuboid.GetShape().all());
    ASSERT_TRUE(cuboid.IsSymmetric());
    //
    const Shape c_cuboid = shapes::Cuboid(nn);
    idx = 0;
    for (int i2 = 0; i2<nn[2]; ++i2) {
        for (int i1 = 0; i1<nn[1]; ++i1) {
            for (int i0 = 0; i0<nn[0]; ++i0) {
                ii = {i0,i1,i2};
                ASSERT_TRUE((c_cuboid[ii]));
                ASSERT_TRUE(c_cuboid[idx]);
                ++idx;
            }
        }
    }
    ASSERT_TRUE(c_cuboid.GetShape().all());
    ASSERT_TRUE(c_cuboid.IsSymmetric());
}

TEST(ShapeGTest,Ellipsoid) {
    constexpr int n_dim = 3;
    const std::array<int,n_dim> rr = {7,5,3};
    const std::array<int,n_dim> nn = {15,11,7};
    std::array<int,n_dim> ii;
    std::array<real_t,n_dim> xx;
    real_t rho;
    const int n_vox = std::accumulate(nn.begin(),nn.end(),1,std::multiplies<int>());
    int idx;
    //
    Shape ellipsoid = shapes::Ellipsoid(rr);
    idx = 0;
    for (int i2 = 0; i2<nn[2]; ++i2) {
        for (int i1 = 0; i1<nn[1]; ++i1) {
            for (int i0 = 0; i0<nn[0]; ++i0) {
                ii = {i0,i1,i2};
                xx = {static_cast<real_t>(i0-rr[0]),static_cast<real_t>(i1-rr[1]),static_cast<real_t>(i2-rr[2])};
                rho = 0.0;
                for (int d = 0; d<n_dim; ++d) {
                    rho += xx[d]*xx[d]/rr[d]/rr[d];
                }
                if (rho<=1.0) {
                    ASSERT_TRUE(ellipsoid[ii]);
                    ASSERT_TRUE(ellipsoid[idx]);
                } else {
                    ASSERT_FALSE(ellipsoid[ii]);
                    ASSERT_FALSE(ellipsoid[idx]);
                }
                ++idx;
            }
        }
    }
    ASSERT_TRUE(ellipsoid.IsSymmetric());
    return;
}

TEST(ShapeGTest,Operators) {
    constexpr int n_dim = 3;
    std::array<int,n_dim> nn1{3,5,1};
    Shape cuboid1 = shapes::Cuboid(nn1);
    cuboid1.Pad(0,1,1);
    cuboid1.Pad(2,2,2);
    std::array<int,n_dim> nn2{5,3,1};
    Shape cuboid2 = shapes::Cuboid(nn2);
    cuboid2.Pad(1,1,1);
    cuboid2.Pad(2,2,2);
    std::array<int,n_dim> nn3{5,1,5};
    Shape cuboid3 = shapes::Cuboid(nn3);
    cuboid3.Pad(1,2,2);
    // union
    Shape cross = cuboid1+cuboid2+cuboid3;
    std::vector<bool> expected_cross{false,false,false,false,false,
                                     false,false,false,false,false,
                                     true ,true ,true ,true ,true ,
                                     false,false,false,false,false,
                                     false,false,false,false,false,
                                     
                                     false,false,false,false,false,
                                     false,false,false,false,false,
                                     true ,true ,true ,true ,true ,
                                     false,false,false,false,false,
                                     false,false,false,false,false,
                                     
                                     false,true ,true ,true ,false,
                                     true ,true ,true ,true ,true ,
                                     true ,true ,true ,true ,true ,
                                     true ,true ,true ,true ,true ,
                                     false,true ,true ,true ,false,
                                     
                                     false,false,false,false,false,
                                     false,false,false,false,false,
                                     true ,true ,true ,true ,true ,
                                     false,false,false,false,false,
                                     false,false,false,false,false,
                                     
                                     false,false,false,false,false,
                                     false,false,false,false,false,
                                     true ,true ,true ,true ,true ,
                                     false,false,false,false,false,
                                     false,false,false,false,false};
    for (int idx = 0; idx<cross.GetSize(); ++idx) {
        ASSERT_TRUE(cross[idx]==expected_cross[idx]);
    }
    ASSERT_TRUE(cross.IsSymmetric());
    // intersection
    Shape inter = cuboid1&cuboid2&cuboid3;
    std::vector<bool> expected_inter{false,false,false,false,false,
                                     false,false,false,false,false,
                                     false,false,false,false,false,
                                     false,false,false,false,false,
                                     false,false,false,false,false,
                                     
                                     false,false,false,false,false,
                                     false,false,false,false,false,
                                     false,false,false,false,false,
                                     false,false,false,false,false,
                                     false,false,false,false,false,
                                     
                                     false,false,false,false,false,
                                     false,false,false,false,false,
                                     false,true ,true ,true ,false,
                                     false,false,false,false,false,
                                     false,false,false,false,false,
                                     
                                     false,false,false,false,false,
                                     false,false,false,false,false,
                                     false,false,false,false,false,
                                     false,false,false,false,false,
                                     false,false,false,false,false,
                                     
                                     false,false,false,false,false,
                                     false,false,false,false,false,
                                     false,false,false,false,false,
                                     false,false,false,false,false,
                                     false,false,false,false,false};
    for (int idx = 0; idx<inter.GetSize(); ++idx) {
        ASSERT_TRUE(inter[idx]==expected_inter[idx]);
    }
    ASSERT_TRUE(inter.IsSymmetric());
    // difference
    Shape diffe = cross-inter;
    std::vector<bool> expected_diffe{false,false,false,false,false,
                                     false,false,false,false,false,
                                     true ,true ,true ,true ,true ,
                                     false,false,false,false,false,
                                     false,false,false,false,false,
                                     
                                     false,false,false,false,false,
                                     false,false,false,false,false,
                                     true ,true ,true ,true ,true ,
                                     false,false,false,false,false,
                                     false,false,false,false,false,
                                     
                                     false,true ,true ,true ,false,
                                     true ,true ,true ,true ,true ,
                                     true ,false,false,false,true ,
                                     true ,true ,true ,true ,true ,
                                     false,true ,true ,true ,false,
                                     
                                     false,false,false,false,false,
                                     false,false,false,false,false,
                                     true ,true ,true ,true ,true ,
                                     false,false,false,false,false,
                                     false,false,false,false,false,
                                     
                                     false,false,false,false,false,
                                     false,false,false,false,false,
                                     true ,true ,true ,true ,true ,
                                     false,false,false,false,false,
                                     false,false,false,false,false};
    for (int idx = 0; idx<diffe.GetSize(); ++idx) {
        ASSERT_TRUE(diffe[idx]==expected_diffe[idx]);
    }
    ASSERT_TRUE(diffe.IsSymmetric());
}

TEST(ShapeGTest,Pad) {
    constexpr int n_dim1 = 1;
    constexpr int n_dim2 = 2;
    constexpr int n_dim3 = 3;
    const std::array<int,n_dim1> nn1 = {3};
    const std::array<int,n_dim2> nn2 = {3,2};
    const std::array<int,n_dim3> nn3 = {3,2,1};
    Shape cuboid1 = shapes::Cuboid(nn1);
    Shape cuboid2 = shapes::Cuboid(nn2);
    Shape cuboid3 = shapes::Cuboid(nn3);
    const int left = 1;
    const int right = 2;
    //
    cuboid1.Pad(0,left,right);
    std::vector<bool> expected1{false,true,true,true,false,false};
    for (int idx = 0; idx<cuboid1.GetSize(); ++idx) {
        ASSERT_TRUE(cuboid1[idx]==expected1[idx]);
    }
    ASSERT_FALSE(cuboid1.IsSymmetric());
    //
    cuboid2.Pad(0,left,right);
    cuboid2.Pad(1,left,right);
    std::vector<bool> expected2{false,false,false,false,false,false,
                                false,true ,true ,true ,false,false,
                                false,true, true ,true ,false,false,
                                false,false,false,false,false,false,
                                false,false,false,false,false,false};
    for (int idx = 0; idx<cuboid2.GetSize(); ++idx) {
        ASSERT_TRUE(cuboid2[idx]==expected2[idx]);
    }
    ASSERT_FALSE(cuboid1.IsSymmetric());
    //
    cuboid3.Pad(0,left,right);
    cuboid3.Pad(1,left,right);
    cuboid3.Pad(2,left,right);
    std::vector<bool> expected3{false,false,false,false,false,false,
                                false,false,false,false,false,false,
                                false,false,false,false,false,false,
                                false,false,false,false,false,false,
                                false,false,false,false,false,false,
                                
                                false,false,false,false,false,false,
                                false,true ,true ,true ,false,false,
                                false,true, true ,true ,false,false,
                                false,false,false,false,false,false,
                                false,false,false,false,false,false,
                                
                                false,false,false,false,false,false,
                                false,false,false,false,false,false,
                                false,false,false,false,false,false,
                                false,false,false,false,false,false,
                                false,false,false,false,false,false,
                                
                                false,false,false,false,false,false,
                                false,false,false,false,false,false,
                                false,false,false,false,false,false,
                                false,false,false,false,false,false,
                                false,false,false,false,false,false};
    for (int idx = 0; idx<cuboid3.GetSize(); ++idx) {
        ASSERT_TRUE(cuboid3[idx]==expected3[idx]);
    }
    ASSERT_FALSE(cuboid1.IsSymmetric());
    return;
}

TEST(ShapeGTest,Shrink) {
    constexpr int n_dim1 = 1;
    constexpr int n_dim2 = 2;
    constexpr int n_dim3 = 3;
    const std::array<int,n_dim1> nn1 = {3};
    const std::array<int,n_dim2> nn2 = {3,2};
    const std::array<int,n_dim3> nn3 = {3,2,1};
    Shape cuboid1 = shapes::Cuboid(nn1);
    cuboid1.Pad(0,1,2);
    Shape cuboid2 = shapes::Cuboid(nn2);
    cuboid2.Pad(0,1,2);
    cuboid2.Pad(1,1,2);
    Shape cuboid3 = shapes::Cuboid(nn3);
    cuboid3.Pad(0,1,2);
    cuboid3.Pad(1,1,2);
    cuboid3.Pad(2,1,2);
    const int left = 2;
    const int right = 1;
    //
    cuboid1.Shrink(0,left,right);
    std::vector<bool> expected1{true,true,false};
    for (int idx = 0; idx<cuboid1.GetSize(); ++idx) {
        ASSERT_TRUE(cuboid1[idx]==expected1[idx]);
    }
    ASSERT_FALSE(cuboid1.IsSymmetric());
    //
    cuboid2.Shrink(0,left,right);
    cuboid2.Shrink(1,left,right);
    std::vector<bool> expected2{true ,true ,false,
                                false,false,false};
    for (int idx = 0; idx<cuboid2.GetSize(); ++idx) {
        ASSERT_TRUE(cuboid2[idx]==expected2[idx]);
    }
    ASSERT_FALSE(cuboid1.IsSymmetric());
    //
    cuboid3.Shrink(0,left,right);
    cuboid3.Shrink(1,left,right);
    cuboid3.Shrink(2,left,right);
    const std::array<int,n_dim3> xnn3 = {nn3[0]+left+right,nn3[1]+left+right,nn3[2]+left+right};
    std::vector<bool> expected3{false,false,false,
                                false,false,false};
    for (int idx = 0; idx<cuboid3.GetSize(); ++idx) {
        ASSERT_TRUE(cuboid3[idx]==expected3[idx]);
    }
    ASSERT_FALSE(cuboid1.IsSymmetric());
    return;
}
