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

TEST(FilterSavitzkyGolayGTest,GetWindow) {
    eptlib::Shape window = eptlib::shapes::Ellipsoid(1,2,3);
    eptlib::filter::SavitzkyGolay sg_filter(1.0,1.0,1.0, window, 2);
    auto& sg_window = sg_filter.GetWindow();
    ASSERT_EQ(sg_window.GetVolume(), window.GetVolume());
    ASSERT_EQ(sg_window.GetSize(0), window.GetSize(0));
    ASSERT_EQ(sg_window.GetSize(1), window.GetSize(1));
    ASSERT_EQ(sg_window.GetSize(2), window.GetSize(2));
    ASSERT_EQ(sg_window.GetNVox(), window.GetNVox());
}

TEST(FilterSavitzkyGolayGTest,GetZeroOrderDerivativeKernel) {
    eptlib::Shape window = eptlib::shapes::Cross(1,1,1);
    eptlib::filter::SavitzkyGolay sg_filter(1.0,1.0,1.0, window, 2);
    auto& zero_order_derivative = sg_filter.GetZeroOrderDerivativeKernel();
    for (size_t idx = 0; idx<zero_order_derivative.size(); ++idx) {
        if (idx == zero_order_derivative.size()/2) {
            ASSERT_NEAR(zero_order_derivative[idx], 1.0, 1e-15);
        } else {
            ASSERT_NEAR(zero_order_derivative[idx], 0.0, 1e-15);
        }
    }
}

TEST(FilterSavitzkyGolayGTest,GetFirstOrderDerivativeKernel) {
    eptlib::Shape window = eptlib::shapes::Cross(1,1,1);
    eptlib::filter::SavitzkyGolay sg_filter(1.0,1.0,1.0, window, 2);
    auto& first_order_derivative_x = sg_filter.GetFirstOrderDerivativeKernel(0);
    auto& first_order_derivative_y = sg_filter.GetFirstOrderDerivativeKernel(1);
    auto& first_order_derivative_z = sg_filter.GetFirstOrderDerivativeKernel(2);
    size_t idx = 0;
    for (size_t k = 0; k<window.GetSize(2); ++k) {
        for (size_t j = 0; j<window.GetSize(1); ++j) {
            for (size_t i = 0; i<window.GetSize(0); ++i) {
                if (window(i,j,k)) {
                    if (i == 0) {
                        ASSERT_NEAR(first_order_derivative_x[idx], -0.5, 1e-15);
                    } else if (i == 2) {
                        ASSERT_NEAR(first_order_derivative_x[idx],  0.5, 1e-15);
                    } else {
                        ASSERT_NEAR(first_order_derivative_x[idx],  0.0, 1e-15);
                    }

                    if (j == 0) {
                        ASSERT_NEAR(first_order_derivative_y[idx], -0.5, 1e-15);
                    } else if (j == 2) {
                        ASSERT_NEAR(first_order_derivative_y[idx],  0.5, 1e-15);
                    } else {
                        ASSERT_NEAR(first_order_derivative_y[idx],  0.0, 1e-15);
                    }

                    if (k == 0) {
                        ASSERT_NEAR(first_order_derivative_z[idx], -0.5, 1e-15);
                    } else if (k == 2) {
                        ASSERT_NEAR(first_order_derivative_z[idx],  0.5, 1e-15);
                    } else {
                        ASSERT_NEAR(first_order_derivative_z[idx],  0.0, 1e-15);
                    }

                    ++idx;
                }
            }
        }
    }
}

TEST(FilterSavitzkyGolayGTest,GetSecondOrderDerivativeKernel) {
    eptlib::Shape window = eptlib::shapes::Cross(1,1,1);
    eptlib::filter::SavitzkyGolay sg_filter(1.0,1.0,1.0, window, 2);
    auto& second_order_derivative_x = sg_filter.GetSecondOrderDerivativeKernel(0);
    auto& second_order_derivative_y = sg_filter.GetSecondOrderDerivativeKernel(1);
    auto& second_order_derivative_z = sg_filter.GetSecondOrderDerivativeKernel(2);
    size_t idx = 0;
    for (size_t k = 0; k<window.GetSize(2); ++k) {
        for (size_t j = 0; j<window.GetSize(1); ++j) {
            for (size_t i = 0; i<window.GetSize(0); ++i) {
                if (window(i,j,k)) {
                    if (j == 1 && k == 1) {
                        ASSERT_NEAR(second_order_derivative_x[idx], i==1 ? -2.0 : 1.0, 1e-14);
                    } else {
                        ASSERT_NEAR(second_order_derivative_x[idx], 0.0, 1e-14);
                    }

                    if (i == 1 && k == 1) {
                        ASSERT_NEAR(second_order_derivative_y[idx], j==1 ? -2.0 : 1.0, 1e-14);
                    } else {
                        ASSERT_NEAR(second_order_derivative_y[idx], 0.0, 1e-14);
                    }

                    if (i == 1 && j == 1) {
                        ASSERT_NEAR(second_order_derivative_z[idx], k==1 ? -2.0 : 1.0, 1e-14);
                    } else {
                        ASSERT_NEAR(second_order_derivative_z[idx], 0.0, 1e-14);
                    }
                    ++idx;
                }
            }
        }
    }
}

TEST(FilterSavitzkyGolayGTest,ZeroOrderDerivative) {
    eptlib::Shape window = eptlib::shapes::Cross(1,1,1);
    eptlib::filter::SavitzkyGolay sg_filter(1.0,1.0,1.0, window, 2);
    std::vector<double> crop(window.GetVolume());
    std::iota(crop.begin(), crop.end(), 1.0);
    ASSERT_NEAR(sg_filter.ZeroOrderDerivative(crop), 4.0, 1e-15);
    ASSERT_NEAR(sg_filter.GetFilter<double>(eptlib::filter::DifferentialOperator::Field)(crop), 4.0, 1e-15);
}

TEST(FilterSavitzkyGolayGTest,FirstOrderDerivative) {
    eptlib::Shape window = eptlib::shapes::Cross(1,1,1);
    eptlib::filter::SavitzkyGolay sg_filter(1.0,1.0,1.0, window, 2);
    std::vector<double> crop(window.GetVolume());
    std::iota(crop.begin(), crop.end(), 1.0);
    ASSERT_NEAR(sg_filter.FirstOrderDerivative(0, crop), 1.0, 1e-15);
    ASSERT_NEAR(sg_filter.FirstOrderDerivative(1, crop), 2.0, 1e-15);
    ASSERT_NEAR(sg_filter.FirstOrderDerivative(2, crop), 3.0, 1e-15);
    ASSERT_NEAR(sg_filter.GetFilter<double>(eptlib::filter::DifferentialOperator::GradientX)(crop), 1.0, 1e-15);
    ASSERT_NEAR(sg_filter.GetFilter<double>(eptlib::filter::DifferentialOperator::GradientY)(crop), 2.0, 1e-15);
    ASSERT_NEAR(sg_filter.GetFilter<double>(eptlib::filter::DifferentialOperator::GradientZ)(crop), 3.0, 1e-15);
}

TEST(FilterSavitzkyGolayGTest,SecondOrderDerivative) {
    eptlib::Shape window = eptlib::shapes::Cross(1,1,1);
    eptlib::filter::SavitzkyGolay sg_filter(1.0,1.0,1.0, window, 2);
    std::vector<double> crop(window.GetVolume());
    size_t idx = 0;
    for (size_t k = 0; k<window.GetSize(2); ++k) {
        for (size_t j = 0; j<window.GetSize(1); ++j) {
            for (size_t i = 0; i<window.GetSize(0); ++i) {
                if (window(i,j,k)) {
                    crop[idx] = i*i + j*j + k*k;
                    ++idx;
                }
            }
        }
    }
    ASSERT_NEAR(sg_filter.SecondOrderDerivative(0, crop), 2.0, 1e-14);
    ASSERT_NEAR(sg_filter.SecondOrderDerivative(1, crop), 2.0, 1e-14);
    ASSERT_NEAR(sg_filter.SecondOrderDerivative(2, crop), 2.0, 1e-14);
    ASSERT_NEAR(sg_filter.Laplacian(crop), 6.0, 1e-14);
    ASSERT_NEAR(sg_filter.GetFilter<double>(eptlib::filter::DifferentialOperator::GradientXX)(crop), 2.0, 1e-14);
    ASSERT_NEAR(sg_filter.GetFilter<double>(eptlib::filter::DifferentialOperator::GradientYY)(crop), 2.0, 1e-14);
    ASSERT_NEAR(sg_filter.GetFilter<double>(eptlib::filter::DifferentialOperator::GradientZZ)(crop), 2.0, 1e-14);
    ASSERT_NEAR(sg_filter.GetFilter<double>(eptlib::filter::DifferentialOperator::Laplacian)(crop), 6.0, 1e-14);
}

TEST(FilterSavitzkyGolayGTest,SavitzkyGolayApply) {
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

    const eptlib::Shape window = eptlib::shapes::Cross(2,2,2);
    const size_t degree = 3;
    eptlib::filter::SavitzkyGolay sg_filter(d0,d1,d2, window, degree);

    eptlib::Image<double> lapl_constant_field (n0,n1,n2);
    eptlib::Image<double> lapl_linear_field   (n0,n1,n2);
    eptlib::Image<double> lapl_quadratic_field(n0,n1,n2);

    eptlib::EPTlibError error;

    error = sg_filter.Apply(eptlib::filter::DifferentialOperator::Laplacian, &lapl_constant_field, constant_field);
    ASSERT_EQ(error, eptlib::EPTlibError::Success);

    error = sg_filter.Apply(eptlib::filter::DifferentialOperator::Laplacian, &lapl_linear_field, linear_field);
    ASSERT_EQ(error, eptlib::EPTlibError::Success);

    error = sg_filter.Apply(eptlib::filter::DifferentialOperator::Laplacian, &lapl_quadratic_field, quadratic_field);
    ASSERT_EQ(error, eptlib::EPTlibError::Success);

    for (int k = 0; k<n2; ++k) {
        for (int j = 0; j<n1; ++j) {
            for (int i = 0; i<n0; ++i) {
                ASSERT_NEAR(lapl_constant_field(i,j,k), 0.0, 1e-13);
                ASSERT_NEAR(lapl_linear_field  (i,j,k), 0.0, 1e-13);
                if (i==0 || i==n0-1 || j==0 || j==n1-1 || k==0 || k==n2-1 ||
                    i==1 || i==n0-2 || j==1 || j==n1-2 || k==1 || k==n2-2) {
                    ASSERT_NEAR(lapl_quadratic_field(i,j,k), 0.0, 1e-13);
                } else {
                    ASSERT_NEAR(lapl_quadratic_field(i,j,k), 6.0, 1e-13);
                }
            }
        }
    }
}

TEST(FilterSavitzkyGolayGTest,WrappedPhase) {
    const double d0 = 1.0e-3;
    const double d1 = 1.0e-3;
    const double d2 = 1.0e-3;

    eptlib::Shape window = eptlib::shapes::Cross(1,1,1);
    eptlib::filter::SavitzkyGolay sg_filter(d0,d1,d2, window, 2);
    std::vector<double> crop(window.GetVolume());
    
    size_t idx = 0;
    for (size_t k = 0; k<window.GetSize(2); ++k) {
        for (size_t j = 0; j<window.GetSize(1); ++j) {
            for (size_t i = 0; i<window.GetSize(0); ++i) {
                if (window(i, j, k)) {
                    double x = i*d0;
                    double y = j*d1;
                    double z = k*d2;
                    crop[idx] = std::arg(std::exp(std::complex<double>(0.0, x*x + y*y + z*z)));
                    ++idx;
                }
            }
        }
    }

    ASSERT_NEAR(sg_filter.GetFilterWrappedPhase(eptlib::filter::DifferentialOperator::Field)(crop), crop[3], 1e-15);
    
    ASSERT_NEAR(sg_filter.GetFilterWrappedPhase(eptlib::filter::DifferentialOperator::GradientX)(crop), 2.0*d0, 1e-14);
    ASSERT_NEAR(sg_filter.GetFilterWrappedPhase(eptlib::filter::DifferentialOperator::GradientY)(crop), 2.0*d1, 1e-14);
    ASSERT_NEAR(sg_filter.GetFilterWrappedPhase(eptlib::filter::DifferentialOperator::GradientZ)(crop), 2.0*d2, 1e-14);

    ASSERT_NEAR(sg_filter.GetFilterWrappedPhase(eptlib::filter::DifferentialOperator::GradientXX)(crop), 2.0, 1e-10);
    ASSERT_NEAR(sg_filter.GetFilterWrappedPhase(eptlib::filter::DifferentialOperator::GradientYY)(crop), 2.0, 1e-10);
    ASSERT_NEAR(sg_filter.GetFilterWrappedPhase(eptlib::filter::DifferentialOperator::GradientZZ)(crop), 2.0, 1e-10);
    ASSERT_NEAR(sg_filter.GetFilterWrappedPhase(eptlib::filter::DifferentialOperator::Laplacian)(crop),  6.0, 1e-10);
}

TEST(FilterSavitzkyGolayGTest,SavitzkyGolayApplyWrappedPhase) {
    const size_t n0 = 10;
    const size_t n1 = 10;
    const size_t n2 = 10;

    const double d0 = 1.0e-3;
    const double d1 = 1.0e-3;
    const double d2 = 1.0e-3;

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

    const eptlib::Shape window = eptlib::shapes::Cross(2,2,2);
    const size_t degree = 3;
    eptlib::filter::SavitzkyGolay sg_filter(d0,d1,d2, window, degree);

    eptlib::Image<double> lapl_constant_field (n0,n1,n2);
    eptlib::Image<double> lapl_linear_field   (n0,n1,n2);
    eptlib::Image<double> lapl_quadratic_field(n0,n1,n2);

    eptlib::EPTlibError error;

    error = sg_filter.ApplyWrappedPhase(eptlib::filter::DifferentialOperator::Laplacian, &lapl_constant_field, constant_field);
    ASSERT_EQ(error, eptlib::EPTlibError::Success);

    error = sg_filter.ApplyWrappedPhase(eptlib::filter::DifferentialOperator::Laplacian, &lapl_linear_field, linear_field);
    ASSERT_EQ(error, eptlib::EPTlibError::Success);

    error = sg_filter.ApplyWrappedPhase(eptlib::filter::DifferentialOperator::Laplacian, &lapl_quadratic_field, quadratic_field);
    ASSERT_EQ(error, eptlib::EPTlibError::Success);

    for (int k = 0; k<n2; ++k) {
        for (int j = 0; j<n1; ++j) {
            for (int i = 0; i<n0; ++i) {
                ASSERT_NEAR(lapl_constant_field(i,j,k), 0.0, 1e-10);
                ASSERT_NEAR(lapl_linear_field  (i,j,k), 0.0, 1e-10);
                if (i==0 || i==n0-1 || j==0 || j==n1-1 || k==0 || k==n2-1 ||
                    i==1 || i==n0-2 || j==1 || j==n1-2 || k==1 || k==n2-2) {
                    ASSERT_NEAR(lapl_quadratic_field(i,j,k), 0.0, 1e-8);
                } else {
                    ASSERT_NEAR(lapl_quadratic_field(i,j,k), 6.0, 1e-8);
                }
            }
        }
    }
}
