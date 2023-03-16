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

#include "eptlib/shape.h"

namespace eptlib {

namespace filter {

    /**
     * @brief Class implementing the Savitzky-Golay filter in an anatomically
     *     aware fashion.
     * 
     * The anatomical Savitzky-Golay filter consists in a local polynomial
     * approximation of the image using only voxels reasonably belonging to
     * the same anatomical tissue. Then, the derivatives of the polynomial
     * approximation are computed analytically.
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
        private:
            /// Resolution in meter of the voxels along each direction.
            std::array<double,N_DIM> dd_;
            /// Mask over which apply the filter.
            Shape window_;
            /// Degree of the fitting polynomial.
            size_t degree_;
    };

}  // namespace filter

}  // namespace eptlib

#endif  // EPTLIB_FILTER_ANATOMICAL_SAVITZKY_GOLAY_H_
