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

#include <functional>
#include <tuple>
#include <type_traits>
#include <vector>

#include "eptlib/image.h"
#include "eptlib/shape.h"
#include "eptlib/util.h"

namespace eptlib {

namespace filter {

    /**
     * @brief Apply a moving window filter on a three-dimensional image.
     * 
     * @tparam Scalar numerical type of the image data.
     * @tparam Filter functional type of the filter to be applied.
     * 
     * @param dst image where the filter result is written.
     * @param src image to which the filter is applied.
     * @param window mask over which apply the filter.
     * @param filter filter to be applied.
     * @param variance optional image where the second filter result (if it exists) is written (default: `nullptr').
     * @param supporting_images optional supporting images to be used for an improved filter (default: `{}', empty list).
     * 
     * @return a Success or a WrongDataFormat if the argument sizes are inconsistent.
     * 
     * `Filter' must be a callable function with one of the following signatures:
     * 
     *  - `Scalar filter(const std::vector<Scalar>&)';
     * 
     *  - `std::tuple<Scalar,double> filter(const std::vector<Scalar>&)`, with `variance!=nullptr';
     * 
     *  - `Scalar filter(const std::vector<Scalar>&, const std::vector<double>&)', with `supporting_images.size()!=0';
     * 
     *  - `std::tuple<Scalar,double> filter(const std::vector<Scalar>&, const std::vector<double>&)', with `variance!=nullptr' and `supporting_images.size()!=0'.
     * 
     * The filter always receives in input the image `src' cropped in the filter window
     * and provides in output the value to be assigned to image `dst'. Additionally,
     * it can receive the `supporting_images' (cropped and listed in a single vector by
     * alternating the entries of each image) and can provide a second output (`variance').
     */
    template <typename Scalar, typename Filter>
    EPTlibError MovingWindow(Image<Scalar> *dst, const Image<Scalar> &src, const Shape &window, const Filter &filter,
        Image<Scalar> *variance = nullptr,
        const std::initializer_list<const Image<double> *> supporting_images = {}) {
        // define compile-time variables about the Filter function
        constexpr bool filter_with_variance = std::is_same_v<FunctionReturnType<Filter>, std::tuple<Scalar, Scalar> >;
        constexpr bool filter_with_supporting_images  = std::is_invocable_v<Filter, const std::vector<Scalar>&, const std::vector<double>&>;
        if ((filter_with_variance && !variance) || (filter_with_supporting_images && supporting_images.size()==0)) {
            return EPTlibError::WrongDataFormat;
        }
        // initialize the runtime variables
        const size_t n0 = src.GetSize(0);
        const size_t n1 = src.GetSize(1);
        const size_t n2 = src.GetSize(2);
        if (dst->GetSize(0)!=n0 || dst->GetSize(1)!=n1 || dst->GetSize(2)!=n2) {
            return EPTlibError::WrongDataFormat;
        }
        const size_t m0 = window.GetSize(0);
        const size_t m1 = window.GetSize(1);
        const size_t m2 = window.GetSize(2);
        const size_t r0 = m0/2;
        const size_t r1 = m1/2;
        const size_t r2 = m2/2;
        const size_t num_supporting_images = filter_with_supporting_images ? supporting_images.size() : 0;
        // loop over the destination voxels
        #pragma omp parallel for
        for (int idx = 0; idx < src.GetNVox(); ++idx) {
            std::vector<Scalar> src_crop;
            std::vector<double> supporting_images_crop;
            src_crop.reserve(window.GetVolume());
            supporting_images_crop.reserve(window.GetVolume() * num_supporting_images);
            // get the sub-indices from the linear index
            size_t i0 = idx % n0;
            size_t tmp = idx / n0;
            size_t i1 = tmp % n1;
            size_t i2 = tmp / n1;
            // loop over the window voxels
            size_t iw2;
            ptrdiff_t ic2;
            for (iw2 = 0, ic2 = static_cast<ptrdiff_t>(i2)-r2; iw2<m2; ++iw2, ++ic2) {
                size_t iw1;
                ptrdiff_t ic1;
                for (iw1 = 0, ic1 = static_cast<ptrdiff_t>(i1)-r1; iw1<m1; ++iw1, ++ic1) {
                    size_t iw0;
                    ptrdiff_t ic0;
                    for (iw0 = 0, ic0 = static_cast<ptrdiff_t>(i0)-r0; iw0<m0; ++iw0, ++ic0) {
                        if (window(iw0, iw1, iw2)) {
                            if (ic0 < 0 || ic1 < 0 || ic2 < 0 || ic0 >= n0 || ic1 >= n1 || ic2 >= n2) {
                                src_crop.push_back(static_cast<Scalar>(0.0));
                                if constexpr (filter_with_supporting_images) {
                                    for (size_t idx_supporting_image = 0; idx_supporting_image < num_supporting_images; ++idx_supporting_image) {
                                        supporting_images_crop.push_back(-1.0);
                                    }
                                }
                            } else {
                                src_crop.push_back(src(ic0, ic1, ic2));
                                if constexpr (filter_with_supporting_images) {
                                    for (auto supporting_image = supporting_images.begin(); supporting_image != supporting_images.end(); ++supporting_image) {
                                        supporting_images_crop.push_back((*supporting_image)->At(ic0, ic1, ic2));
                                    }
                                }
                            }
                        }
                    }
                }
            }
            // apply the filter
            //   the filter application is obtained in two steps
            //   1) the reference to where to write the result is initialised by assigning an rvalue to output
            //   2) the result is written by assigning an lvalue to output
            using output_t = std::conditional_t<filter_with_variance, std::tuple<Scalar&,Scalar&>, Scalar&>;
            output_t output = ConstexprIf<filter_with_variance>(std::tie(dst->At(i0,i1,i2),variance->At(i0,i1,i2)), std::ref(dst->At(i0,i1,i2)));
            if constexpr (filter_with_supporting_images) {
                output = std::invoke(filter, src_crop, supporting_images_crop);
            } else {
                output = std::invoke(filter, src_crop);
            }
        }
        return EPTlibError::Success;
    }

}  // namespace filter

}  // namespace eptlib

#endif  // EPTLIB_FILTER_MOVING_WINDOW_H_
