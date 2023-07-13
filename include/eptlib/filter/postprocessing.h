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

#ifndef EPTLIB_FILTER_POSTPROCESSING_H_
#define EPTLIB_FILTER_POSTPROCESSING_H_

#include <algorithm>

#include "eptlib/image.h"
#include "eptlib/shape.h"
#include "eptlib/util.h"

#include "eptlib/filter/weight_functions.h"

namespace eptlib {

namespace filter {

    template <typename Scalar>
    EPTlibError Postprocessing(Image<Scalar> *dst, const Image<Scalar> &src, const Shape &window,
        const Image<double> &variance, const Image<double> &ref_img, const double max) {
        auto filter = [&](const std::vector<Scalar> &src_crop, const std::vector<double> &ref_img_crop, const std::vector<double> &variance_crop) -> Scalar {
            size_t idx0 = ref_img_crop.size() / 2;
            double ref0 = ref_img_crop[idx0];
            // remove data from the computation
            std::vector<Scalar> src_crop_tmp(0);
            std::vector<double> var_crop_tmp(0);
            src_crop_tmp.reserve(src_crop.size());
            for (size_t idx = 0; idx < src_crop.size(); ++idx) {
                const double ref = ref_img_crop[idx];
                const Scalar src = src_crop[idx];
                const double var = variance_crop[idx];
                if (HardThreshold(std::abs(ref-ref0)/(ref+ref0)*2.0, 0.10) == 0.0) {
                    continue;
                }
                if (std::isnan(src)) {
                    continue;
                }
                if (src < 0.0 || src > max) {
                    continue;
                }
                src_crop_tmp.push_back(src);
                var_crop_tmp.push_back(var);
            }
            // sort the values based on variance
            size_t n = var_crop_tmp.size() / 10;
            std::vector<size_t> indices(var_crop_tmp.size());
            std::iota(indices.begin(), indices.end(), 0);
            std::partial_sort(indices.begin(), indices.begin()+n, indices.end(), [&](int a, int b) -> bool { return var_crop_tmp[a] < var_crop_tmp[b]; });
            // compute the weighted mean
            return std::accumulate(indices.begin(), indices.begin()+n, 0.0, [&](double a, int b) -> double { return a + src_crop_tmp[b]/n; });
//            // compute the weights based on the ref_img
//            size_t idx0 = ref_img_crop.size() / 2;
//            double ref0 = ref_img_crop[idx0];
//            std::vector<double> weights(ref_img_crop.size());
//            std::transform(ref_img_crop.begin(), ref_img_crop.end(), weights.begin(),
//                [ref0](const double x) -> double {
//                    return Gaussian(x-ref0, 0.05);
//                }
//            );
//            // combine the weights with the standard deviations
//            std::transform(variance_crop.begin(), variance_crop.end(), weights.begin(), weights.begin(),
//                [](const double s, const double x) -> double {
//                    if (s > 0.0) {
//                        return x/std::sqrt(s);
//                    }
//                    return 0.0;
//                }
//            );
//            // cut the negative values
//            std::transform(src_crop.begin(), src_crop.end(), weights.begin(), weights.begin(),
//                [max](const double s, const double x) -> double {
//                    if (std::isnan(s)) {
//                        return 0.0;
//                    }
//                    if (s < 0.0 || s > max) {
//                        return 0.0;
//                    }
//                    return x;
//                }
//            );
//            // remove nan from computation
//            std::vector<double> src_crop_tmp(src_crop.size());
//            std::transform(src_crop.begin(), src_crop.end(), src_crop_tmp.begin(),
//                [](const double x) -> double {
//                    if (std::isnan(x)) {
//                        return 0.0;
//                    }
//                    return x;
//                }
//            );
//            // normalize the weights
//            double mass = Sum(weights);
//            std::transform(weights.begin(), weights.end(), weights.begin(),
//                [mass](const double x) -> double {
//                    return x/mass;
//                }
//            );
//            // apply the filter
//            return std::inner_product(src_crop_tmp.begin(), src_crop_tmp.end(), weights.begin(), 0.0);
        };

        //

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
        // loop over the destination voxels
        #pragma omp parallel for collapse(3)
        for (size_t i2 = r2; i2<n2-r2; ++i2) {
            for (size_t i1 = r1; i1<n1-r1; ++i1) {
                for (size_t i0 = r0; i0<n0-r0; ++i0) {
                    std::vector<Scalar> src_crop(window.GetVolume());
                    std::vector<double> ref_img_crop(window.GetVolume());
                    std::vector<double> variance_crop(window.GetVolume());
                    double ref_img0 = ref_img(i0,i1,i2);
                    // loop over the window voxels
                    size_t idx_crop = 0;
                    for (size_t iw2 = 0; iw2<m2; ++iw2) {
                        for (size_t iw1 = 0; iw1<m1; ++iw1) {
                            for (size_t iw0 = 0; iw0<m0; ++iw0) {
                                if (window(iw0, iw1, iw2)) {
                                    src_crop[idx_crop] = src(i0-r0+iw0, i1-r1+iw1, i2-r2+iw2);
                                    ref_img_crop[idx_crop] = ref_img(i0-r0+iw0, i1-r1+iw1, i2-r2+iw2);
                                    variance_crop[idx_crop] = variance(i0-r0+iw0, i1-r1+iw1, i2-r2+iw2);
                                    ++idx_crop;
                                }
                            }
                        }
                    }
                    // apply the filter
                    dst->At(i0,i1,i2) = std::invoke(filter, src_crop, ref_img_crop, variance_crop);
                }
            }
        }
        return EPTlibError::Success;
    }

}  // namespace filter

}  // namespace eptlib

#endif  // EPTLIB_FILTER_POSTPROCESSING_H_
