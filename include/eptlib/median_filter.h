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

#ifndef EPTLIB_MEDIAN_FILTER_H_
#define EPTLIB_MEDIAN_FILTER_H_

#include <vector>

#include "eptlib/shape.h"
#include "eptlib/util.h"

namespace eptlib {

/**
 * Class for the application of median filters.
 */
class MedianFilter {
    public:
        /**
         * Constructor.
         * 
         * @param shape mask over which apply the median filter.
         */
        MedianFilter(const Shape &shape);

        /**
         * Apply the median filter to an input field.
         * 
         * @param[out] dst pointer to the output destination.
         * @param[in] src pointer to the input source.
         * @param[in] nn number of voxels in each direction.
         * @param[in] img pointer to the reference image. If it is not
         *     nullptr, then it is used as a reference.
         * 
         * @return a Success or Unknown error.
         */
        EPTlibError ApplyFilter(double *dst, const double *src,
            const std::array<int,NDIM> &nn, const double *img = nullptr);
    private:
        /// Shape of the kernel for median filter application.
        Shape shape_;
        /// Number of voxels in the shape.
        int m_vol_;
};

}  // namespace eptlib

#endif  // EPTLIB_MEDIAN_FILTER_H_
