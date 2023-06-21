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
     * @param variance optional image where the second filter result (if it exists) is written. (Default: `nullptr')
     * @param ref_img optional reference image to be used for an improved filter. (Default: `nullptr')
     * 
     * @return a Success or a WrongDataFormat if the argument sizes are inconsistent.
     * 
     * `Filter' must be a callable function with one of the following signatures:
     * 
     *  - `Scalar filter(const std::vector<Scalar>&)';
     * 
     *  - `std::tuple<Scalar,double> filter(const std::vector<Scalar>&)`, with `variance!=nullptr';
     * 
     *  - `Scalar filter(const std::vector<Scalar>&, const std::vector<double>&)', with `ref_img!=nullptr';
     * 
     *  - `std::tuple<Scalar,double> filter(const std::vector<Scalar>&, const std::vector<double>&)', with `variance!=nullptr' and `ref_img!=nullptr'.
     * 
     * The filter always receives in input the image `src' cropped in the filter window
     * and provides in output the value to be assigned to image `dst'. Additionally,
     * it can receive a reference image (`ref_img') and can provide a second
     * output (`variance').
     */
    template <typename Scalar, typename Filter>
    EPTlibError MovingWindow(Image<Scalar> *dst, const Image<Scalar> &src, const Shape &window, const Filter &filter,
        Image<double> *variance = nullptr,
        const Image<double> *ref_img = nullptr) {
        // define compile-time variables about the Filter function
        constexpr bool filter_with_variance = std::is_same_v<FunctionReturnType<Filter>, std::tuple<Scalar, double> >;
        constexpr bool filter_with_ref_img  = std::is_invocable_v<Filter, const std::vector<Scalar>&, const std::vector<double>&>;
        if ((filter_with_variance && !variance) || (filter_with_ref_img && !ref_img)) {
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
        const size_t num_ref_imgs = filter_with_ref_img ? 1 : 0;
        // loop over the destination voxels
        #pragma omp parallel for collapse(3)
        for (size_t i2 = r2; i2<n2-r2; ++i2) {
            for (size_t i1 = r1; i1<n1-r1; ++i1) {
                for (size_t i0 = r0; i0<n0-r0; ++i0) {
                    std::vector<Scalar> src_crop(window.GetVolume());
                    std::vector<double> ref_img_crop(window.GetVolume() * num_ref_imgs);
                    double ref_img0 = 0.0;
                    if constexpr (filter_with_ref_img) {
                        ref_img0 = ref_img->At(i0,i1,i2);
                    }
                    // loop over the window voxels
                    size_t idx_crop = 0;
                    for (size_t iw2 = 0; iw2<m2; ++iw2) {
                        for (size_t iw1 = 0; iw1<m1; ++iw1) {
                            for (size_t iw0 = 0; iw0<m0; ++iw0) {
                                if (window(iw0, iw1, iw2)) {
                                    src_crop[idx_crop] = src(i0-r0+iw0, i1-r1+iw1, i2-r2+iw2);
                                    if constexpr (filter_with_ref_img) {
                                        ref_img_crop[idx_crop] = ref_img->At(i0-r0+iw0, i1-r1+iw1, i2-r2+iw2) / ref_img0;
                                    }
                                    ++idx_crop;
                                }
                            }
                        }
                    }
                    // apply the filter
                    //   the filter application is obtained in two steps
                    //   1) the reference to where to write the result is initialised by assigning an rvalue to output
                    //   2) the result is written by assigning an lvalue to output
                    using output_t = std::conditional_t<filter_with_variance, std::tuple<Scalar&,double&>, Scalar&>;
                    output_t output = ConstexprIf<filter_with_variance>(std::tie(dst->At(i0,i1,i2),variance->At(i0,i1,i2)), std::ref(dst->At(i0,i1,i2)));
                    if constexpr (filter_with_ref_img) {
                        output = std::invoke(filter, src_crop, ref_img_crop);
                    } else {
                        output = std::invoke(filter, src_crop);
                    }
                }
            }
        }
        return EPTlibError::Success;
    }

}  // namespace filter

}  // namespace eptlib

#endif  // EPTLIB_FILTER_MOVING_WINDOW_H_
