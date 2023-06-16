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

#ifndef EPTLIB_FILTER_ANATOMICAL_SAVITZKY_GOLAY_H_
#define EPTLIB_FILTER_ANATOMICAL_SAVITZKY_GOLAY_H_

#include <array>
#include <functional>
#include <tuple>
#include <vector>

#include "eptlib/differential_operator.h"
#include "eptlib/image.h"
#include "eptlib/shape.h"

#include "eptlib/linalg/matrix.h"
#include "eptlib/linalg/qr.h"
#include "eptlib/linalg/regression.h"
#include "eptlib/linalg/vector.h"

#include "eptlib/filter/moving_window.h"

namespace eptlib {

namespace filter {

    /**
     * @brief Class implementing an Savitzky-Golay filter whose kernel is
     *     automatically adapted to the sample anatomy.
     * 
     * Savitzky-Golay filter consists in a local polynomial approximation of
     * which the derivatives are computed analytically. The kernel of the
     * filter is selected on the basis of a reference image.
     */
    class AnatomicalSavitzkyGolay {
        public:
            /**
             * @brief Construct a new anatomical Savitzky-Golay object.
             * 
             * @param d0 resolution in meter along direction x.
             * @param d1 resolution in meter along direction y.
             * @param d2 resolution in meter along direction z.
             * @param window mask over which apply the filter.
             * @param degree degree of the fitting polynomial (default: 2).
             */
            AnatomicalSavitzkyGolay(const double d0, const double d1, const double d2,
                const Shape &window, const size_t degree = 2);

            /**
             * @brief Destroy the anatomical Savitzky-Golay object.
             */
            virtual ~AnatomicalSavitzkyGolay();

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
             * @brief Get the design matrix.
             * 
             * @return the design matrix.
             */
            inline const eptlib::linalg::Matrix<double>& GetDesignMatrix() const {
                return design_matrix_;
            }

            /**
             * @brief Compute the weights for the polynomial fitting according to the relative
             *     contrasts in the reference image.
             * 
             * @param ref_img_crop reference image whose contrast is used to define the weights.
             * 
             * @return a vector of weights for the polynomial fitting.
             */
            std::vector<double> ComputeWeights(const std::vector<double> &ref_img_crop) const;

            /**
             * @brief Apply the weights to the design matrix and the forcing term of the polynomial fitting.
             * 
             * @tparam Scalar numerical type of the data.
             * 
             * @param weights vector of weights for the polynomial fitting.
             * @param src_crop forcint term that will be altered according to the weights.
             * @param design_matrix design matrix that will be altered according to the weights.
             */
            template <typename Scalar>
            void ApplyWeights(const std::vector<double> &weights, std::vector<Scalar> *src_crop, 
                eptlib::linalg::Matrix<double> *design_matrix) const {
                for (size_t idx = 0; idx < weights.size(); ++idx) {
                    (*src_crop)[idx] *= weights[idx];
                }
                for (size_t col = 0; col < design_matrix->GetNCol(); ++col) {
                    for (size_t idx = 0; idx < weights.size(); ++idx) {
                        (*design_matrix)(idx, col) *= weights[idx];
                    }
                }
                return;
            }

            /**
             * @brief Extract the zero order derivative (field approximation) from the fitted polynomial coefficients.
             * 
             * @tparam Scalar numerical type of the data.
             * 
             * @param a coefficients of the fitted polynomial.
             * 
             * @return the approximated zero order derivative (field approximation).
             */
            template <typename Scalar>
            Scalar ExtractZeroOrderDerivative(const std::vector<Scalar> &a) const {
                return a[derivative_indices_[0]];
            }

            /**
             * @brief Extract the first order derivative along direction d from the fitted polynomial coefficients.
             * 
             * @tparam Scalar numerical type of the data.
             * 
             * @param d derivative direction.
             * @param a coefficients of the fitted polynomial.
             * 
             * @return the approximated first order derivative.
             */
            template <typename Scalar>
            Scalar ExtractFirstOrderDerivative(const size_t d, const std::vector<Scalar> &a) const {
                return a[derivative_indices_[1+d]];
            }

            /**
             * @brief Extract the second order derivative along direction d from the fitted polynomial coefficients.
             * 
             * @tparam Scalar numerical type of the data.
             * 
             * @param d derivative direction.
             * @param a coefficients of the fitted polynomial.
             * 
             * @return the approximated second order derivative.
             */
            template <typename Scalar>
            Scalar ExtractSecondOrderDerivative(const size_t d, const std::vector<Scalar> &a) const {
                return a[derivative_indices_[4+d]];
            }

            /**
             * @brief Extract the Laplacian from the fitted polynomial coefficients.
             * 
             * @tparam Scalar numerical type of the data.
             * 
             * @param a coefficients of the fitted polynomial.
             * 
             * @return the approximated Laplacian.
             */
            template <typename Scalar>
            Scalar ExtractLaplacian(const std::vector<Scalar> &a) const {
                return 2.0 * (a[derivative_indices_[4]] + a[derivative_indices_[5]] + a[derivative_indices_[6]]);
            }

            /**
             * @brief Get the extractor function corresponding to the selected differential operator.
             * 
             * @tparam Scalar numerical type of the data.
             * 
             * @param differential_operator differential operator to be applied.
             * 
             * @return the extractor function corresponding to the selected differential operator.
             */
            template <typename Scalar>
            inline std::function<Scalar(const std::vector<Scalar>&)> GetExtractor(const DifferentialOperator differential_operator) const {
                switch (differential_operator) {
                    case DifferentialOperator::Field:
                        return [&](const std::vector<Scalar> &a) -> Scalar {
                            return this->ExtractZeroOrderDerivative(a);
                        };
                    case DifferentialOperator::GradientX:
                        return [&](const std::vector<Scalar> &a) -> Scalar {
                            return this->ExtractFirstOrderDerivative(0, a);
                        };
                    case DifferentialOperator::GradientY:
                        return [&](const std::vector<Scalar> &a) -> Scalar {
                            return this->ExtractFirstOrderDerivative(1, a);
                        };
                    case DifferentialOperator::GradientZ:
                        return [&](const std::vector<Scalar> &a) -> Scalar {
                            return this->ExtractFirstOrderDerivative(2, a);
                        };
                    case DifferentialOperator::GradientXX:
                        return [&](const std::vector<Scalar> &a) -> Scalar {
                            return this->ExtractSecondOrderDerivative(0, a);
                        };
                    case DifferentialOperator::GradientYY:
                        return [&](const std::vector<Scalar> &a) -> Scalar {
                            return this->ExtractSecondOrderDerivative(1, a);
                        };
                    case DifferentialOperator::GradientZZ:
                        return [&](const std::vector<Scalar> &a) -> Scalar {
                            return this->ExtractSecondOrderDerivative(2, a);
                        };
                    case DifferentialOperator::Laplacian:
                        return [&](const std::vector<Scalar> &a) -> Scalar {
                            return this->ExtractLaplacian(a);
                        };
                    default:
                        return [ ](const std::vector<Scalar> &a) -> Scalar {
                            return 0.0;
                        };
                };
            }

            /**
             * @brief Apply the Savitzky-Golay filter to estimate the selected differential operator.
             * 
             * @tparam Scalar numerical type of the data.
             * 
             * @param differential_operator differential operator to be applied.
             * @param dst image where the filter result is written.
             * @param src image to which the filter is applied.
             * @param ref_img reference image used for an improved filter.
             * 
             * @return a Success or a WrongDataFormat if the argument sizes are inconsistent.
             */
            template <typename Scalar>
            EPTlibError Apply(const DifferentialOperator differential_operator, Image<Scalar> *dst,
                const Image<Scalar> &src, const Image<double> &ref_img) const {
                auto filter = [&](std::vector<Scalar> src_crop, const std::vector<double> &ref_img_crop) -> Scalar {
                    eptlib::linalg::Matrix<double> design_matrix = this->design_matrix_;
                    std::vector<double> weights = this->ComputeWeights(ref_img_crop);
                    this->ApplyWeights(weights, &src_crop, &design_matrix);
                    auto [a, chi2n] = eptlib::linalg::LinearRegression(design_matrix, src_crop);
                    return this->GetExtractor<Scalar>(differential_operator)(a);
                };
                return MovingWindow(dst, src, window_, filter, nullptr, &ref_img);
            }

            /**
             * @brief Apply the Savitzky-Golay filter ot estimate the selected differential operator
             *     and compute the output variance.
             * 
             * @tparam Scalar numerical type of the data.
             * 
             * @param differential_operator differential operator to be applied.
             * @param dst image where the filter result is written.
             * @param variance image where the output variance is written.
             * @param src image to which the filter is applied.
             * @param ref_img reference image used for an improved filter.
             * 
             * @return a Success or a WrongDataFormat if the argument sizes are inconsistent.
             */
            template <typename Scalar>
            EPTlibError Apply(const DifferentialOperator differential_operator, Image<Scalar> *dst,
                Image<double> *variance, const Image<Scalar> &src, const Image<double> &ref_img) const {
                const size_t n_col = design_matrix_.GetNCol();
                std::vector<double> w(n_col, 0.0);
                switch (differential_operator) {
                    case DifferentialOperator::Field:
                        w[derivative_indices_[0]] = 1.0;
                        break;
                    case DifferentialOperator::GradientX:
                        w[derivative_indices_[1]] = 1.0;
                        break;
                    case DifferentialOperator::GradientY:
                        w[derivative_indices_[2]] = 1.0;
                        break;
                    case DifferentialOperator::GradientZ:
                        w[derivative_indices_[3]] = 1.0;
                        break;
                    case DifferentialOperator::GradientXX:
                        w[derivative_indices_[4]] = 2.0;
                        break;
                    case DifferentialOperator::GradientYY:
                        w[derivative_indices_[5]] = 2.0;
                        break;
                    case DifferentialOperator::GradientZZ:
                        w[derivative_indices_[6]] = 2.0;
                        break;
                    case DifferentialOperator::Laplacian:
                        w[derivative_indices_[4]] = 2.0;
                        w[derivative_indices_[5]] = 2.0;
                        w[derivative_indices_[6]] = 2.0;
                        break;
                    default:
                        break;
                }   
                auto filter = [&](std::vector<Scalar> src_crop, const std::vector<double> &ref_img_crop) -> std::tuple<Scalar, double> {
                    // get the derivative
                    eptlib::linalg::Matrix<double> design_matrix = this->design_matrix_;
                    std::vector<double> weights = this->ComputeWeights(ref_img_crop);
                    this->ApplyWeights(weights, &src_crop, &design_matrix);
                    auto [QR, p] = eptlib::linalg::QRDecomposition(design_matrix);
                    auto [a, chi] = eptlib::linalg::QRSolve(QR, src_crop);
                    eptlib::linalg::Permute(a.begin(), a.end(), p);
                    Scalar derivative = this->GetExtractor<Scalar>(differential_operator)(a);
                    // get the variance
                    size_t m = design_matrix.GetNRow() - std::count(weights.begin(), weights.end(), 0.0);
                    size_t rank = eptlib::linalg::QRGetRank(QR);
                    double chi2n = chi*chi/(m-rank);
                    std::vector<double> w_p(n_col);
                    for (size_t col = 0; col < n_col; ++col) {
                        w_p[col] = w[p[col]];
                    }
                    std::vector<double> r = eptlib::linalg::RTransposeSolve(QR, w_p);
                    double var = eptlib::linalg::Norm2(r.begin(), r.end());
                    var *= var*chi2n;
                    return {derivative, var};
                };
                return MovingWindow(dst, src, window_, filter, variance, &ref_img);
            }
        private:
            /// Resolution in meter of the voxels along each direction.
            std::array<double,N_DIM> dd_;
            /// Mask over which apply the filter.
            Shape window_;
            /// Degree of the fitting polynomial.
            size_t degree_;

            /// Complete design matrix.
            eptlib::linalg::Matrix<double> design_matrix_;

            /// Where the derivative contributions are found in the design matrix columns.
            static constexpr std::array<size_t, 7> derivative_indices_ = {0, 1, 2, 3, 4, 6, 9};
    };

}  // namespace filter

}  // namespace eptlib

#endif  // EPTLIB_FILTER_ANATOMICAL_SAVITZKY_GOLAY_H_
