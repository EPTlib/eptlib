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

#include <numeric>
#include <vector>

#include "eptlib/image.h"
#include "eptlib/shape.h"
#include "eptlib/util.h"

#include "eptlib/filter/moving_window.h"

namespace eptlib {

    class SavitzkyGolay {
        public:
            /**
             * @brief Construct a new Savitzky Golay object.
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
             * @brief Compute the derivative by applying the filter in a crop.
             * 
             * @tparam Scalar numerical type of the data.
             * 
             * @param crop data values to which apply the filter.
             * 
             * @return Scalar the approximated derivative.
             */
            template <typename Scalar>
            Scalar Filter(const std::vector<Scalar> &crop) {
                return std::inner_product(lapl_.begin(), lapl_.end(), crop.begin(), 0.0);
            }

            template <typename Scalar>
            EPTlibError Apply(Image<Scalar> *dst, const Image<Scalar> &src) {
                auto filter = [&](const std::vector<Scalar> &crop) -> Scalar {
                    return this->Filter(crop);
                };
                return MovingWindow(dst, src, window_, filter);
            }
        private:
            /// Mask over which apply the filter.
            Shape window_;
            /// Coefficients of the Laplacian approximation.
            std::vector<double> lapl_;
    };

}

#endif  // EPTLIB_FILTER_SAVITZKY_GOLAY_H_
