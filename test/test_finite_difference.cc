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

#include "eptlib/finite_difference.h"

#include <algorithm>
#include <array>
#include <functional>
#include <numeric>
#include <typeinfo>
#include <vector>

using namespace eptlib;

TEST(FiniteDifferenceGTest,SavitzkyGolayLaplacian) {
    const std::array<int,NDIM> nn = {7,5,3};
    const std::array<double,NDIM> dd = {1.0,2.0,3.0};
    std::array<int,NDIM> ii;
    const int n_vox = std::accumulate(nn.begin(),nn.end(),1,std::multiplies<int>());
    //
    std::array<int,NDIM> rr = {1,1,1};
    std::array<int,NDIM> rr2 = {1,0,0};
    Shape cross = shapes::Cross(rr);
    Shape cross2 = shapes::Cross(rr2);
    ASSERT_TRUE(cross.IsSymmetric());
    FDSavitzkyGolayFilter fd_lapl(cross);
    FDSavitzkyGolayFilter fd_lapl2(cross2);
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
    std::vector<double> c_lapl2(n_vox);
    std::vector<double> l_lapl2(n_vox);
    std::vector<double> q_lapl2(n_vox);
    DifferentialOperator diff_op = DifferentialOperator::Laplacian;
    fd_lapl.Apply(diff_op,c_lapl.data(),c_field.data(),nn,dd);
    fd_lapl.Apply(diff_op,l_lapl.data(),l_field.data(),nn,dd);
    fd_lapl.Apply(diff_op,q_lapl.data(),q_field.data(),nn,dd);
    fd_lapl2.Apply(diff_op,c_lapl2.data(),c_field.data(),nn,dd);
    fd_lapl2.Apply(diff_op,l_lapl2.data(),l_field.data(),nn,dd);
    fd_lapl2.Apply(diff_op,q_lapl2.data(),q_field.data(),nn,dd);
    //
    std::vector<double> c_lapl_chi2(n_vox);
    std::vector<double> l_lapl_chi2(n_vox);
    std::vector<double> q_lapl_chi2(n_vox);
    std::vector<double> c_lapl2_chi2(n_vox);
    std::vector<double> l_lapl2_chi2(n_vox);
    std::vector<double> q_lapl2_chi2(n_vox);
    std::vector<double> chi2(n_vox);
    fd_lapl.Apply(diff_op,c_lapl_chi2.data(),chi2.data(),c_field.data(),nn,dd);
    fd_lapl.Apply(diff_op,l_lapl_chi2.data(),chi2.data(),l_field.data(),nn,dd);
    fd_lapl.Apply(diff_op,q_lapl_chi2.data(),chi2.data(),q_field.data(),nn,dd);
    EPTlibError error = fd_lapl2.Apply(diff_op,c_lapl2_chi2.data(),chi2.data(),c_field.data(),nn,dd);
    ASSERT_TRUE(error==EPTlibError::Success);
    error = fd_lapl2.Apply(diff_op,l_lapl2_chi2.data(),chi2.data(),l_field.data(),nn,dd);
    ASSERT_TRUE(error==EPTlibError::Success);
    error = fd_lapl2.Apply(diff_op,q_lapl2_chi2.data(),chi2.data(),q_field.data(),nn,dd);
    ASSERT_TRUE(error==EPTlibError::Success);
    //
    for (int idx = 0; idx<n_vox; ++idx) {
        ASSERT_NEAR(c_lapl[idx],0.0,1e-14);
        ASSERT_NEAR(l_lapl[idx],0.0,1e-13);
        ASSERT_NEAR(c_lapl_chi2[idx],0.0,1e-14);
        ASSERT_NEAR(l_lapl_chi2[idx],0.0,1e-13);
        ASSERT_NEAR(c_lapl2[idx],0.0,1e-14);
        ASSERT_NEAR(l_lapl2[idx],0.0,1e-13);
        ASSERT_NEAR(c_lapl2_chi2[idx],0.0,1e-14);
        ASSERT_NEAR(l_lapl2_chi2[idx],0.0,1e-13);
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
            ASSERT_NEAR(q_lapl_chi2[idx],2.0*NDIM,1e-12);
            ASSERT_NEAR(q_lapl2[idx],2.0,1e-12);
            ASSERT_NEAR(q_lapl2_chi2[idx],2.0,1e-12);
        } else {
            ASSERT_NEAR(q_lapl[idx],0.0,1e-12);
            ASSERT_NEAR(q_lapl_chi2[idx],0.0,1e-12);
        }
    }
}
TEST(FiniteDifferenceGTest,SavitzkyGolayGradient) {
    const std::array<int,NDIM> nn = {7,5,3};
    const std::array<double,NDIM> dd = {1.0,2.0,3.0};
    std::array<int,NDIM> ii;
    const int n_vox = std::accumulate(nn.begin(),nn.end(),1,std::multiplies<int>());
    //
    std::array<int,NDIM> rr = {1,1,1};
    std::array<int,NDIM> rr2 = {1,0,0};
    Shape cross = shapes::Cross(rr);
    Shape cross2 = shapes::Cross(rr2);
    ASSERT_TRUE(cross.IsSymmetric());
    FDSavitzkyGolayFilter fd_grad(cross);
    FDSavitzkyGolayFilter fd_grad2(cross2);
    //
    std::vector<double> c_field(n_vox);
    std::vector<double> l_field(n_vox);
    std::vector<double> q_field(n_vox);
    {
        int idx = 0;
        for (int i2 = 0; i2<nn[2]; ++i2) {
            for (int i1 = 0; i1<nn[1]; ++i1) {
                for (int i0 = 0; i0<nn[0]; ++i0) {
                    c_field[idx] = 1.0;
                    l_field[idx] = i0*dd[0] + i1*dd[1] + i2*dd[2];
                    q_field[idx] = i0*i0*dd[0]*dd[0] + i1*i1*dd[1]*dd[1] + i2*i2*dd[2]*dd[2];
                    ++idx;
                }
            }
        }
    }
    //
    std::array<std::vector<double>,NDIM> c_grad;
    std::array<std::vector<double>,NDIM> l_grad;
    std::array<std::vector<double>,NDIM> q_grad;
    std::array<std::vector<double>,NDIM> c_grad2;
    std::array<std::vector<double>,NDIM> l_grad2;
    std::array<std::vector<double>,NDIM> q_grad2;
    for (int d = 0; d<NDIM; ++d) {
        c_grad[d].resize(n_vox);
        l_grad[d].resize(n_vox);
        q_grad[d].resize(n_vox);
        c_grad2[d].resize(n_vox);
        l_grad2[d].resize(n_vox);
        q_grad2[d].resize(n_vox);
        DifferentialOperator diff_op = static_cast<DifferentialOperator>(d+1);
        fd_grad.Apply(diff_op,c_grad[d].data(),c_field.data(),nn,dd);
        fd_grad.Apply(diff_op,l_grad[d].data(),l_field.data(),nn,dd);
        fd_grad.Apply(diff_op,q_grad[d].data(),q_field.data(),nn,dd);
        fd_grad2.Apply(diff_op,c_grad2[d].data(),c_field.data(),nn,dd);
        fd_grad2.Apply(diff_op,l_grad2[d].data(),l_field.data(),nn,dd);
        fd_grad2.Apply(diff_op,q_grad2[d].data(),q_field.data(),nn,dd);
    }
    //
    std::array<std::vector<double>,NDIM> c_grad_chi2;
    std::array<std::vector<double>,NDIM> l_grad_chi2;
    std::array<std::vector<double>,NDIM> q_grad_chi2;
    std::array<std::vector<double>,NDIM> c_grad2_chi2;
    std::array<std::vector<double>,NDIM> l_grad2_chi2;
    std::array<std::vector<double>,NDIM> q_grad2_chi2;
    std::vector<double> chi2(n_vox);
    for (int d = 0; d<NDIM; ++d) {
        c_grad_chi2[d].resize(n_vox);
        l_grad_chi2[d].resize(n_vox);
        q_grad_chi2[d].resize(n_vox);
        c_grad2_chi2[d].resize(n_vox);
        l_grad2_chi2[d].resize(n_vox);
        q_grad2_chi2[d].resize(n_vox);
        DifferentialOperator diff_op = static_cast<DifferentialOperator>(d+1);
        fd_grad.Apply(diff_op,c_grad_chi2[d].data(),chi2.data(),c_field.data(),nn,dd);
        fd_grad.Apply(diff_op,l_grad_chi2[d].data(),chi2.data(),l_field.data(),nn,dd);
        fd_grad.Apply(diff_op,q_grad_chi2[d].data(),chi2.data(),q_field.data(),nn,dd);
        EPTlibError error = fd_grad2.Apply(diff_op,c_grad2_chi2[d].data(),chi2.data(),c_field.data(),nn,dd);
        if (d==0) {
            ASSERT_TRUE(error==EPTlibError::Success);
        } else {
            ASSERT_TRUE(error==EPTlibError::MissingData);
        }
        error = fd_grad2.Apply(diff_op,l_grad2_chi2[d].data(),chi2.data(),l_field.data(),nn,dd);
        if (d==0) {
            ASSERT_TRUE(error==EPTlibError::Success);
        } else {
            ASSERT_TRUE(error==EPTlibError::MissingData);
        }
        error = fd_grad2.Apply(diff_op,q_grad2_chi2[d].data(),chi2.data(),q_field.data(),nn,dd);
        if (d==0) {
            ASSERT_TRUE(error==EPTlibError::Success);
        } else {
            ASSERT_TRUE(error==EPTlibError::MissingData);
        }
    }
    //
    for (int d = 0; d<NDIM; ++d) {
        for (int idx = 0; idx<n_vox; ++idx) {
            ASSERT_NEAR(c_grad[d][idx],0.0,1e-14);
            ASSERT_NEAR(c_grad_chi2[d][idx],0.0,1e-14);
            ASSERT_NEAR(c_grad2[d][idx],0.0,1e-14);
            ASSERT_NEAR(c_grad2_chi2[d][idx],0.0,1e-14);
            IdxToMultiIdx(ii,idx,nn);
            bool kernel_in_domain = true;
            for (int d = 0; d<NDIM; ++d) {
                if (ii[d]==0 || ii[d]==nn[d]-1) {
                    kernel_in_domain = false;
                    break;
                }
            }
            if (kernel_in_domain) {
                ASSERT_NEAR(l_grad[d][idx],1.0,1e-13);
                ASSERT_NEAR(l_grad_chi2[d][idx],1.0,1e-13);
                ASSERT_NEAR(q_grad[d][idx],2.0*ii[d]*dd[d],1e-12);
                ASSERT_NEAR(q_grad_chi2[d][idx],2.0*ii[d]*dd[d],1e-12);
                if (d==0) {
                    ASSERT_NEAR(l_grad2[d][idx],1.0,1e-13);
                    ASSERT_NEAR(l_grad2_chi2[d][idx],1.0,1e-13);
                    ASSERT_NEAR(q_grad2[d][idx],2.0*ii[d]*dd[d],1e-12);
                    ASSERT_NEAR(q_grad2_chi2[d][idx],2.0*ii[d]*dd[d],1e-12);
                }
            } else {
                ASSERT_NEAR(l_grad[d][idx],0.0,1e-13);
                ASSERT_NEAR(l_grad_chi2[d][idx],0.0,1e-13);
                ASSERT_NEAR(q_grad[d][idx],0.0,1e-12);
                ASSERT_NEAR(q_grad_chi2[d][idx],0.0,1e-12);
            }
        }
    }
}

TEST(FiniteDifferenceGTest,AsymmetricFDSavitzkyGolayLaplacian) {
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
    FDSavitzkyGolayFilter fd_lapl(cross);
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
    DifferentialOperator diff_op = DifferentialOperator::Laplacian;
    fd_lapl.Apply(diff_op,c_lapl.data(),c_field.data(),nn,dd);
    fd_lapl.Apply(diff_op,l_lapl.data(),l_field.data(),nn,dd);
    fd_lapl.Apply(diff_op,q_lapl.data(),q_field.data(),nn,dd);
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
TEST(FiniteDifferenceGTest,AsymmetricFDSavitzkyGolayGradient) {
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
    FDSavitzkyGolayFilter fd_grad(cross);
    //
    std::vector<double> c_field(n_vox);
    std::vector<double> l_field(n_vox);
    std::vector<double> q_field(n_vox);
    for (int idx = 0; idx<n_vox; ++idx) {
        IdxToMultiIdx(ii,idx,nn);
        c_field[idx] = 1.0;
        q_field[idx] = 0.0;
        for (int d = 0; d<NDIM; ++d) {
            l_field[idx] += ii[d]*dd[d];
            q_field[idx] += ii[d]*ii[d]*dd[d]*dd[d];
        }
    }
    //
    std::array<std::vector<double>,NDIM> c_grad;
    std::array<std::vector<double>,NDIM> l_grad;
    std::array<std::vector<double>,NDIM> q_grad;
    for (int d = 0; d<NDIM; ++d) {
        c_grad[d].resize(n_vox);
        l_grad[d].resize(n_vox);
        q_grad[d].resize(n_vox);
        DifferentialOperator diff_op = static_cast<DifferentialOperator>(d+1);
        fd_grad.Apply(diff_op,c_grad[d].data(),c_field.data(),nn,dd);
        fd_grad.Apply(diff_op,l_grad[d].data(),l_field.data(),nn,dd);
        fd_grad.Apply(diff_op,q_grad[d].data(),q_field.data(),nn,dd);
    }
    //
    for (int d = 0; d<NDIM; ++d) {
        for (int idx = 0; idx<n_vox; ++idx) {
            ASSERT_NEAR(c_grad[d][idx],0.0,1e-14);
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
                ASSERT_NEAR(l_grad[d][idx],1.0,1e-13);
                ASSERT_NEAR(q_grad[d][idx],2.0*ii[d]*dd[d],1e-12);
            } else {
                ASSERT_NEAR(l_grad[d][idx],0.0,1e-13);
                ASSERT_NEAR(q_grad[d][idx],0.0,1e-12);
            }
        }
    }
}

TEST(FiniteDifferenceGTest,HigherOrderInterpolation) {
    const std::array<int,NDIM> nn = {7,5,5};
    const std::array<double,NDIM> dd = {1.0,2.0,3.0};
    const int degree = 4;
    std::array<int,NDIM> ii;
    const int n_vox = std::accumulate(nn.begin(),nn.end(),1,std::multiplies<int>());
    //
    std::array<int,NDIM> rr = {2,2,2};
    Shape cube = shapes::CuboidR(rr);
    ASSERT_TRUE(cube.IsSymmetric());
    FDSavitzkyGolayFilter fd_lapl(cube,degree);
    //
    std::vector<double> field(n_vox);
    for (int idx = 0; idx<n_vox; ++idx) {
        IdxToMultiIdx(ii,idx,nn);
        field[idx] = 0.0;
        for (int d = 0; d<NDIM; ++d) {
            field[idx] += ii[d]*ii[d]*dd[d]*dd[d];
            field[idx] += ii[d]*ii[d]*ii[d]*dd[d]*dd[d]*dd[d];
        }
    }
    //
    std::vector<double> lapl(n_vox);
    DifferentialOperator diff_op = DifferentialOperator::Laplacian;
    fd_lapl.Apply(diff_op,lapl.data(),field.data(),nn,dd);
    //
    std::vector<double> lapl_chi2(n_vox);
    std::vector<double> chi2(n_vox);
    EPTlibError error = fd_lapl.Apply(diff_op,lapl_chi2.data(),chi2.data(),field.data(),nn,dd);
    ASSERT_TRUE(error==EPTlibError::Success);
    //
    for (int idx = 0; idx<n_vox; ++idx) {
        IdxToMultiIdx(ii,idx,nn);
        bool kernel_in_domain = true;
        for (int d = 0; d<NDIM; ++d) {
            if (ii[d]-rr[d]<0 || ii[d]+rr[d]>=nn[d]) {
                kernel_in_domain = false;
                break;
            }
        }
        if (kernel_in_domain) {
            double tmp = 0.0;
            for (int d = 0; d<NDIM; ++d) {
                tmp += 2.0 + 6.0*ii[d]*dd[d];
            }
            ASSERT_NEAR(lapl[idx],tmp,1e-12);
            ASSERT_NEAR(lapl_chi2[idx],tmp,1e-12);
        } else {
            ASSERT_NEAR(lapl[idx],0.0,1e-12);
            ASSERT_NEAR(lapl_chi2[idx],0.0,1e-12);
        }
    }
}
