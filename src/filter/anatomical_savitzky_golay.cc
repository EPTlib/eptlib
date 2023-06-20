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

#include "eptlib/filter/anatomical_savitzky_golay.h"

#include "eptlib/polynomial/fitting.h"

// AnatomicalSavitzkyGolay constructor
eptlib::filter::AnatomicalSavitzkyGolay::
AnatomicalSavitzkyGolay(const double d0, const double d1, const double d2,
    const eptlib::Shape &window, const size_t degree) :
    window_(window),
    degree_(degree) {
    dd_[0] = d0;
    dd_[1] = d1;
    dd_[2] = d2;
    // get coordinates of window voxels
    std::vector<double> x(window_.GetVolume());
    std::vector<double> y(window_.GetVolume());
    std::vector<double> z(window_.GetVolume());
    const double x0 = window_.GetSize(0)/2*dd_[0];
    const double y0 = window_.GetSize(1)/2*dd_[1];
    const double z0 = window_.GetSize(2)/2*dd_[2];
    size_t idx = 0;
    for (size_t k = 0; k < window_.GetSize(2); ++k) {
        for (size_t j = 0; j < window_.GetSize(1); ++j) {
            for (size_t i = 0; i < window_.GetSize(0); ++i) {
                if (window_(i,j,k)) {
                    x[idx] = i*dd_[0] - x0;
                    y[idx] = j*dd_[1] - y0;
                    z[idx] = k*dd_[2] - z0;
                    ++idx;
                }
            }
        }
    }
    // compute the design matrix
    design_matrix_ = eptlib::polynomial::DesignMatrixWithMonomialBasis(x, y, z, degree_);
    return;
}

// AnatomicalSavitzkyGolay virtual destructor
eptlib::filter::AnatomicalSavitzkyGolay::
~AnatomicalSavitzkyGolay() {
    return;
}

// AnatomicalSavitzkyGolay compute the weights for the polynomial fitting
std::vector<double> eptlib::filter::AnatomicalSavitzkyGolay::
ComputeWeights(const std::vector<double> &ref_img_crop) const {
//    return WeightHardThreshold(ref_img_crop, 0.1);
    return WeightGaussian(ref_img_crop, 0.05);
}

//
std::vector<double> eptlib::filter::AnatomicalSavitzkyGolay::
WeightHardThreshold(const std::vector<double> &ref_img_crop, const double threshold) const {
    std::vector<double> weights(ref_img_crop.size());
    size_t idx0 = ref_img_crop.size()/2;
    for (size_t idx = 0; idx < ref_img_crop.size(); ++idx) {
        double relative_contrast = std::abs((ref_img_crop[idx] - ref_img_crop[idx0]) / (ref_img_crop[idx] + ref_img_crop[idx0]) * 2.0);
        if (relative_contrast > threshold) {
            weights[idx] = 0.0;
        } else {
            weights[idx] = 1.0;
        }
    }
    return weights;
}

//
std::vector<double> eptlib::filter::AnatomicalSavitzkyGolay::
WeightGaussian(const std::vector<double> &ref_img_crop, const double sigma) const {
    std::vector<double> weights(ref_img_crop.size());
    size_t idx0 = ref_img_crop.size()/2;
    for (size_t idx = 0; idx < ref_img_crop.size(); ++idx) {
        double delta = ref_img_crop[idx] - ref_img_crop[idx0];
        weights[idx] = std::exp(-delta*delta/2/sigma/sigma);
    }
    return weights;
}
