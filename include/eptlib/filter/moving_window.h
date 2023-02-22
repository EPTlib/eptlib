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

#ifndef EPTLIB_FILTER_MOVING_WINDOW_H_
#define EPTLIB_FILTER_MOVING_WINDOW_H_

#include <vector>

#include "eptlib/image.h"
#include "eptlib/shape.h"
#include "eptlib/util.h"

namespace eptlib {

    template <typename Scalar, typename Filter>
    EPTlibError MovingWindow(Image<Scalar> *dst, const Image<Scalar> &src, const Shape &window, const Filter &filter) {
        const size_t n0 = src.GetSize(0);
        const size_t n1 = src.GetSize(1);
        const size_t n2 = src.GetSize(2);
        const size_t m0 = window.GetSize(0);
        const size_t m1 = window.GetSize(1);
        const size_t m2 = window.GetSize(2);
        const size_t r0 = m0/2;
        const size_t r1 = m1/2;
        const size_t r2 = m2/2;
        // loop over the destination voxels
        for (size_t i2 = r2; i2<n2-r2; ++i2) {
            for (size_t i1 = r1; i1<n1-r1; ++i1) {
                for (size_t i0 = r0; i0<n0-r0; ++i0) {
                    std::vector<Scalar> src_crop(window.GetVolume());
                    // loop over the window voxels
                    size_t idx_crop = 0;
                    for (size_t iw2 = 0; iw2<m2; ++iw2) {
                        for (size_t iw1 = 0; iw1<m1; ++iw1) {
                            for (size_t iw0 = 0; iw0<m0; ++iw0) {
                                if (window(iw0, iw1, iw2)) {
                                    src_crop[idx_crop] = src(i0-r0+iw0, i1-r1+iw1, i2-r2+iw2);
                                    ++idx_crop;
                                }
                            }
                        }
                    }
                    // apply the filter
                    dst->At(i0, i1, i2) = filter(src_crop, window);
                }
            }
        }
        return EPTlibError::Success;
    }

}  // eptlib

#endif  // EPTLIB_FILTER_MOVING_WINDOW_H_
