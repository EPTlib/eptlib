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

#ifndef EPTLIB_FILTER_SAVITZKY_GOLAY_H_
#define EPTLIB_FILTER_SAVITZKY_GOLAY_H_

#include <array>
#include <functional>
#include <numeric>
#include <tuple>
#include <vector>

#include "eptlib/differential_operator.h"
#include "eptlib/image.h"
#include "eptlib/shape.h"
#include "eptlib/util.h"

#include "eptlib/filter/moving_window.h"

#include "eptlib/linalg/linalg_util.h"

namespace eptlib {

namespace filter {

    /**
     * @brief Class implementing the Savitzky-Golay filter.
     * 
     * Savitzky-Golay filter consists in a local polynomial approximation of
     * which the derivatives are computed analytically.
     */
    class SavitzkyGolay {
        public:
            /**
             * @brief Construct a new Savitzky-Golay object.
             * 
             * @param d0 resolution in meter along direction x.
             * @param d1 resolution in meter along direction y.
             * @param d2 resolution in meter along direction z.
             * @param window mask over which apply the filter.
             * @param degree degree of the fitting polynomial (default: 2).
             */
            SavitzkyGolay(const double d0, const double d1, const double d2,
                const Shape &window, const size_t degree = 2);

            /**
             * @brief Destroy the Savitzky Golay object.
             */
            virtual ~SavitzkyGolay();

            /**
             * @brief Get a constant reference to the filter window.
             * 
             * @return a constant reference to the filter window.
             */
            inline const Shape& GetWindow() const {
                return window_;
            }

            /**
             * @brief Get a constant reference to the filter resolution.
             * 
             * @return a constant reference to the resolution.
             */
            inline const std::array<double,N_DIM>& GetResolution() const {
                return dd_;
            }

            /**
             * @brief Get the filter resolution along dimension `d'.
             * 
             * @param d dimension of interest.
             * 
             * @return the resolution along dimension `d'.
             */
            inline double GetResolution(const size_t d) const {
                return dd_[d];
            }

            /**
             * @brief Get the degree of the fitting polynomial.
             * 
             * @return the degree of the fitting polynomial.
             */
            inline size_t GetDegree() const {
                return degree_;
            }

            /**
             * @brief Get a constant reference to the coefficients of the zero order derivative (field approximation).
             * 
             * @return a constant reference to the coefficients of the zero order derivative (field approximation).
             */
            inline const std::vector<double>& GetZeroOrderDerivativeKernel() const {
                return zero_order_derivative_;
            }

            /**
             * @brief Get a constant reference to the coefficients of the first order derivative approximation.
             * 
             * @param d derivative direction.
             * 
             * @return a constant reference to the coefficients of the first order derivative approximation.
             */
            inline const std::vector<double>& GetFirstOrderDerivativeKernel(const size_t d) const {
                return first_order_derivative_[d];
            }

            /**
             * @brief Get a constant reference to the coefficients of the second order derivative approximation.
             * 
             * @param d derivative direction.
             * 
             * @return a constant reference to the coefficients of the second order derivative approximation.
             */
            inline const std::vector<double>& GetSecondOrderDerivativeKernel(const size_t d) const {
                return second_order_derivative_[d];
            }

            /**
             * @brief Compute the zero order derivative (field approximation) by applying the filter in a crop.
             * 
             * @tparam Scalar numerical type of the data.
             * 
             * @param crop data values to which apply the filter.
             * 
             * @return the approximated zero order derivative (field approximation).
             */
            template <typename Scalar>
            inline Scalar ZeroOrderDerivative(const std::vector<Scalar> &crop) const {
                return std::inner_product(zero_order_derivative_.begin(), zero_order_derivative_.end(), crop.begin(), static_cast<Scalar>(0.0));
            }

            /**
             * @brief Compute the first order derivative along direction d by applying the filter in a crop.
             * 
             * @tparam Scalar numerical type of the data.
             * 
             * @param d derivative direction.
             * @param crop data values to which apply the filter.
             * 
             * @return the approximated first order derivative.
             */
            template <typename Scalar>
            inline Scalar FirstOrderDerivative(const size_t d, const std::vector<Scalar> &crop) const {
                return std::inner_product(first_order_derivative_[d].begin(), first_order_derivative_[d].end(), crop.begin(), static_cast<Scalar>(0.0));
            }

            /**
             * @brief Compute the second order derivative along direction d by applying the filter in a crop.
             * 
             * @tparam Scalar numerical type of the data.
             * 
             * @param d derivative direction.
             * @param crop data values to which apply the filter.
             * 
             * @return the approximated second order derivative.
             */
            template <typename Scalar>
            inline Scalar SecondOrderDerivative(const size_t d, const std::vector<Scalar> &crop) const {
                return std::inner_product(second_order_derivative_[d].begin(), second_order_derivative_[d].end(), crop.begin(), static_cast<Scalar>(0.0));
            }

            /**
             * @brief Compute the Laplacian by applying the filter in a crop.
             * 
             * @tparam Scalar numerical type of the data.
             * 
             * @param crop data values to which apply the filter.
             * 
             * @return the approximated Laplacian.
             */
            template <typename Scalar>
            inline Scalar Laplacian(const std::vector<Scalar> &crop) const {
                Scalar laplacian = 0.0;
                laplacian += SecondOrderDerivative(0, crop);
                laplacian += SecondOrderDerivative(1, crop);
                laplacian += SecondOrderDerivative(2, crop);
                return laplacian;
            }

            /**
             * @brief Get the filter function corresponding to the selected differential operator.
             * 
             * @tparam Scalar numerical type of the data.
             * 
             * @param differential_operator differential operator to be applied.
             * 
             * @return the filter function corresponding to the selected differential operator.
             */
            template <typename Scalar>
            inline std::function<Scalar(const std::vector<Scalar>&)> GetFilter(const DifferentialOperator differential_operator) const {
                switch (differential_operator) {
                case DifferentialOperator::Field:
                    return [&](const std::vector<Scalar> &crop) -> Scalar {
                        return this->ZeroOrderDerivative(crop);
                    };
                case DifferentialOperator::GradientX:
                    return [&](const std::vector<Scalar> &crop) -> Scalar {
                        return this->FirstOrderDerivative(0, crop);
                    };
                case DifferentialOperator::GradientY:
                    return [&](const std::vector<Scalar> &crop) -> Scalar {
                        return this->FirstOrderDerivative(1, crop);
                    };
                case DifferentialOperator::GradientZ:
                    return [&](const std::vector<Scalar> &crop) -> Scalar {
                        return this->FirstOrderDerivative(2, crop);
                    };
                case DifferentialOperator::GradientXX:
                    return [&](const std::vector<Scalar> &crop) -> Scalar {
                        return this->SecondOrderDerivative(0, crop);
                    };
                case DifferentialOperator::GradientYY:
                    return [&](const std::vector<Scalar> &crop) -> Scalar {
                        return this->SecondOrderDerivative(1, crop);
                    };
                case DifferentialOperator::GradientZZ:
                    return [&](const std::vector<Scalar> &crop) -> Scalar {
                        return this->SecondOrderDerivative(2, crop);
                    };
                case DifferentialOperator::Laplacian:
                    return [&](const std::vector<Scalar> &crop) -> Scalar {
                        return this->Laplacian(crop);
                    };
                default:
                    return [ ](const std::vector<Scalar> &crop) -> Scalar {
                        return 0.0;
                    };
                };
            }

            /**
             * @brief Compute the variance with which the differential operator is estimated on the crop.
             * 
             * @tparam Scalar numerical type of the data.
             * 
             * @param differential_operator differential operator to be applied.
             * @param crop data values to which the filter is applied.
             * 
             * @return the variance with which the differential operator is estimated.
             */
            template <typename Scalar>
            inline double ComputeVariance(const DifferentialOperator differential_operator, const std::vector<Scalar> &crop) const {
                const size_t n_row = residuals_.size();
                std::vector<Scalar> residual(n_row, 0.0);
                for (size_t k = 0; k<n_row; ++k) {
                    for (size_t row = 0; row<n_row; ++row) {
                        residual[row] += crop[k] * residuals_[k][row];
                    }
                }
                return variance_coefficients_[static_cast<size_t>(differential_operator)] * eptlib::linalg::Norm2(residual.data(), n_row);
            }

            /**
             * @brief Apply the Savitzky-Golay filter to estimate the selected differential operator.
             * 
             * @tparam Scalar numerical type of the data.
             * 
             * @param differential_operator differential operator to be applied.
             * @param dst image where the filter result is written.
             * @param src image to which the filter is applied.
             * 
             * @return a Success or a WrongDataFormat if the argument sizes are inconsistent.
             */
            template <typename Scalar>
            EPTlibError Apply(const DifferentialOperator differential_operator, Image<Scalar> *dst, const Image<Scalar> &src) const {
                auto filter = GetFilter<Scalar>(differential_operator);
                return MovingWindow(dst, src, window_, filter);
            }

            /**
             * @brief Apply the Savitzky-Golay filter to estimate the selected differential operator
             *     and compute the output variance.
             * 
             * @tparam Scalar numerical type of the data.
             * 
             * @param differential_operator differential operator to be applied.
             * @param dst image where the filter result is written.
             * @param variance image where the output variance is written
             * @param src image to which the filter is applied.
             * 
             * @return a Success or a WrongDataFormat if the argument sizes are incosistent.
             */
            template <typename Scalar>
            EPTlibError Apply(const DifferentialOperator differential_operator, Image<Scalar> *dst,
                Image<double> *variance, const Image<Scalar> &src) const {
                auto filter = [&](const std::vector<Scalar> &crop) -> std::tuple<Scalar, double> {
                    auto compute_derivative = this->GetFilter<Scalar>(differential_operator);
                    Scalar derivative = compute_derivative(crop);
                    double variance = this->ComputeVariance(differential_operator, crop);
                    return {derivative, variance};
                };
                return MovingWindow(dst, src, window_, filter, variance);
            }

            /**
             * @brief Get the filter function corresponding to the selected differential operator
             *     for a wrapped phase map.
             * 
             * @param differential_operator differential operator to be applied.
             * 
             * @return the filter function corresponding to the selected differential operator for
             *     a wrapped phase map.
             */
            std::function<double(const std::vector<double>&)> GetFilterWrappedPhase(const DifferentialOperator differential_operator) const;

            /**
             * @brief Apply the Savitzky-Golay filter to estimate the selected differential operator
             *     for a wrapped phase map.
             * 
             * @param differential_operator differential operator to be applied.
             * @param dst image where the filter result is written.
             * @param src image to which the filter is applied.
             * 
             * @return a Success or a WrongDataFormat if the argument sizes are inconsistent.
             */
            EPTlibError ApplyWrappedPhase(const DifferentialOperator differential_operator, Image<double> *dst, const Image<double> &src) const;
        private:
            /// Resolution in meter of the voxels along each direction.
            std::array<double,N_DIM> dd_;
            /// Mask over which apply the filter.
            Shape window_;
            /// Degree of the fitting polynomial.
            size_t degree_;

            /// Coefficients of the zero order derivative (map approximation).
            std::vector<double> zero_order_derivative_;
            /// Coefficients of the first order derivatives.
            std::array<std::vector<double>, N_DIM> first_order_derivative_;
            /// Coefficients of the second order derivatives.
            std::array<std::vector<double>, N_DIM> second_order_derivative_;

            /// Residuals of the fitting, used to compute the variance.
            linalg::MatrixReal residuals_;
            /// Coefficients to rescale the output variance.
            std::array<double, static_cast<size_t>(DifferentialOperator::END)> variance_coefficients_;
    };

}  // namespace filter

}  // namespace eptlib

#endif  // EPTLIB_FILTER_SAVITZKY_GOLAY_H_
