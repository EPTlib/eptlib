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

    /**
     * @brief Apply the median filter to the elements of a three-dimensional window.
     * 
     * @param src_crop data values to which apply the filter.
     * 
     * @return the median value.
     */
    double MedianFilter(const std::vector<double> &src_crop);

    /**
     * @brief Apply the median filter to the elements selected from a three-dimensional
     *     window according to a reference image.
     * 
     * @param src_crop data values to which apply the filter.
     * @param ref_img_crop reference data used for selecting the data values.
     * @param weight_param parameter of the thresholding function (default: 0.05).
     * 
     * @return the median value of selected values.
     */
    double AnatomicalMedianFilter(const std::vector<double> &src_crop, const std::vector<double> &ref_img_crop, const double weight_param = 0.05);

    /**
     * @brief Filter the elements of a three-dimensional window based on an uncertainty index.
     * 
     * @param src_crop data values to which apply the filter.
     * @param uncertainty_crop uncertainty values associated to each data value.
     * @param ratio_of_best_values ratio of best values to be selected for filtering (default: 0.10).
     * 
     * @return the filtered value.
     */
    double UncertainFilter(const std::vector<double> &src_crop, const std::vector<double> &uncertainty_crop, const double ratio_of_best_values = 0.10);

    /**
     * @brief Filter the elements selected from a three-dimensional window based on an uncertainty
     *     index.
     * 
     * @param src_crop data values to which apply the filter.
     * @param supporting_crop supporting values for selecting the data values (reference and
     *     uncertainty values).
     * @param weight_param parameter of the thresholding function (default: 0.05).
     * @param ratio_of_best_values ratio of best values to be selected for filter (default: 0.10).
     * 
     * @return the filtered value.
     */
    double AnatomicalUncertainFilter(const std::vector<double> &src_crop, const std::vector<double> &supporting_crop, const double weight_param = 0.05, const double ratio_of_best_values = 0.10);

    /**
     * @brief Postprocess the input image depending on the amount of available information.
     * 
     * @param dst image where the filter result is written.
     * @param src image to which the filter is applied.
     * @param window mask over which apply the filter.
     * @param reference_image optional reference image used for an improved filter (default:
     *     `nullptr').
     * @param uncertainty optional uncertainty image used for an improved filter (default:
     *     `nullptr').
     * @param weight_param parameter of the thresholding function, used only in combination with
     *     `reference_image' (default: 0.05).
     * @param ratio_of_best_values ratio of best values to be selected for filter, used only in
     *     combination with `uncertainty' (default: 0.10).
     * 
     * @return a Success, a WrongDataFormat if the argument sizes are inconsistent, or a Unknown.
     * 
     * Based on the availability of `reference_image' and `uncertainty', four combinations are
     * possible:
     * 
     * - if both are not available, the MedianFilter is used for postprocessing the data;
     * 
     * - if only the `reference_image' is available, the AnatomicalMedianFilter is used;
     * 
     * - if only the `uncertainty' is available, the UncertainFilter is used;
     * 
     * - if both are available, the AnatomicalUncertainFilter is used.
     */
    EPTlibError Postprocessing(Image<double> *dst, const Image<double> &src, const Shape &window, const Image<double> *reference_image = nullptr, const Image<double> *uncertainty = nullptr, const double weight_param = 0.05, const double ratio_of_best_values = 0.10);

}  // namespace filter

}  // namespace eptlib

#endif  // EPTLIB_FILTER_POSTPROCESSING_H_
