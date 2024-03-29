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

#include "eptlib/filter/postprocessing.h"

#include <algorithm>

#include "eptlib/filter/moving_window.h"
#include "eptlib/filter/weight_functions.h"

namespace {

    // Compute the median.
    double Median(std::vector<double> v) {
        auto m = v.size() / 2;
        std::nth_element(v.begin(), v.begin()+m, v.end());
        return v[m];
    }

}  //

// Apply the median filter.
double eptlib::filter::MedianFilter(const std::vector<double> &src_crop) {
    // remove data from the computation
    std::vector<double> src_selected(0);
    src_selected.reserve(src_crop.size());
    for (size_t idx = 0; idx < src_crop.size(); ++idx) {
        const double src = src_crop[idx];
        if (!std::isfinite(src)) {
            continue;
        }
        src_selected.push_back(src);
    }
    // compute the median
    return src_selected.size()>0 ? ::Median(src_selected) : 0.0;
}

// Apply the median filter to selected values based on the anatomy of a reference image.
double eptlib::filter::AnatomicalMedianFilter(const std::vector<double> &src_crop, const std::vector<double> &ref_img_crop, const double weight_param) {
    // identify the central reference value
    const size_t idx0 = ref_img_crop.size() / 2;
    const double ref0 = ref_img_crop[idx0];
    // remove data from the computation
    std::vector<double> src_selected(0);
    src_selected.reserve(src_crop.size());
    for (size_t idx = 0; idx < src_crop.size(); ++idx) {
        const double ref = ref_img_crop[idx];
        const double src = src_crop[idx];
        if (eptlib::filter::HardThreshold(std::abs((ref - ref0) / std::min(ref, ref0)), weight_param) == 0) {
            continue;
        }
        if (!std::isfinite(src)) {
            continue;
        }
        src_selected.push_back(src);
    }
    // compute the median
    return src_selected.size()>0 ? ::Median(src_selected) : nand;
}

// Filter the elements of a three-dimensional window based on an uncertainty index.
double eptlib::filter::UncertainFilter(const std::vector<double> &src_crop, const std::vector<double> &uncertainty_crop) {
    // remove data from the computation
    std::vector<double> src_selected(0);
    std::vector<double> unc_selected(0);
    src_selected.reserve(src_crop.size());
    unc_selected.reserve(src_crop.size());
    for (size_t idx = 0; idx < src_crop.size(); ++idx) {
        const double src = src_crop[idx];
        const double unc = uncertainty_crop[idx];
        if (!std::isfinite(src) || !std::isfinite(unc)) {
            continue;
        }
        src_selected.push_back(src);
        unc_selected.push_back(unc);
    }
    // compute the mean
    double mass = std::accumulate(unc_selected.begin(), unc_selected.end(), 0.0, [](double a, double b) -> double { return a + 1 / b; });
    return std::inner_product(src_selected.begin(), src_selected.end(), unc_selected.begin(), 0.0, [](double a, double b) -> double { return a + b; }, [&](double a, double b) -> double { return a / b / mass; });
}

// Filter the elements selected from a three-dimensional window based on an uncertainty index.
double eptlib::filter::AnatomicalUncertainFilter(const std::vector<double> &src_crop, const std::vector<double> &supporting_crop, const double weight_param) {
    // identify the central reference value
    const size_t idx0 = src_crop.size() / 2;
    const double ref0 = supporting_crop[idx0 * 2];
    // remove data from the computation
    std::vector<double> src_selected(0);
    std::vector<double> unc_selected(0);
    src_selected.reserve(src_crop.size());
    unc_selected.reserve(src_crop.size());
    for (size_t idx = 0; idx < src_crop.size(); ++idx) {
        const double src = src_crop[idx];
        const double ref = supporting_crop[idx * 2];
        const double unc = supporting_crop[idx * 2 + 1];
        if (eptlib::filter::HardThreshold(std::abs((ref - ref0) / std::min(ref, ref0)), weight_param) == 0) {
            continue;
        }
        if (!std::isfinite(src) || !std::isfinite(unc)) {
            continue;
        }
        src_selected.push_back(src);
        unc_selected.push_back(unc);
    }
    // select the best (less uncertain) values
    double mass = std::accumulate(unc_selected.begin(), unc_selected.end(), 0.0, [](double a, double b) -> double { return a + 1 / b; });
    return std::inner_product(src_selected.begin(), src_selected.end(), unc_selected.begin(), 0.0, [](double a, double b) -> double { return a + b; }, [&](double a, double b) -> double { return a / b / mass; });
}

// Postprocess the input image depending on the amount of available information.
eptlib::EPTlibError eptlib::filter::Postprocessing(eptlib::Image<double> *dst, const eptlib::Image<double> &src,
    const eptlib::Shape &window, const Image<double> *reference_image, const Image<double> *uncertainty,
    const double weight_param) {
    // no additional information available
    if (!reference_image && !uncertainty) {
        return eptlib::filter::MovingWindow(dst, src, window, eptlib::filter::MedianFilter);
    }
    // only reference image available
    if (reference_image && !uncertainty) {
        auto filter = [&] (const std::vector<double> &src_crop, const std::vector<double> &ref_img_crop) -> double {
            return eptlib::filter::AnatomicalMedianFilter(src_crop, ref_img_crop, weight_param);
        };
        return eptlib::filter::MovingWindow<double>(dst, src, window, filter, nullptr, {reference_image});
    }
    // only uncertainty available
    if (!reference_image && uncertainty) {
        auto filter = [&] (const std::vector<double> &src_crop, const std::vector<double> &uncertainty_crop) -> double {
            return eptlib::filter::AnatomicalMedianFilter(src_crop, uncertainty_crop);
        };
        return eptlib::filter::MovingWindow<double>(dst, src, window, filter, nullptr, {uncertainty});
    }
    // both additional information available
    if (reference_image && uncertainty) {
        auto filter = [&] (const std::vector<double> &src_crop, const std::vector<double> &supporting_crop) -> double {
            return eptlib::filter::AnatomicalUncertainFilter(src_crop, supporting_crop, weight_param);
        };
        return eptlib::filter::MovingWindow<double>(dst, src, window, filter, nullptr, {reference_image, uncertainty});
    }
    return EPTlibError::Unknown;
}
