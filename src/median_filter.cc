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

#include "eptlib/median_filter.h"

using namespace eptlib;

// MedianFilter constructor
MedianFilter::
MedianFilter(const Shape &shape) :
    shape_(shape), m_vol_(shape_.GetVolume()) {
    // check to have odd shape size
    int m_vox = std::accumulate(shape.GetSize().begin(),shape.GetSize().end(),1,std::multiplies<int>());
    assert(m_vox%2);
    return;
}

// MedianFilter apply
EPTlibError MedianFilter::
ApplyFilter(double *dst, const double *src, const std::array<int,NDIM> &nn,
    const double *img) {
    const int n_vox = std::accumulate(nn.begin(),nn.end(),1,std::multiplies<int>());
    std::array<int,NDIM> ii;
    std::array<int,NDIM> rr;
    std::array<int,NDIM> inc;
    for (int d = 0; d<NDIM; ++d) {
        rr[d] = shape_.GetSize()[d]/2;
    }
    inc[0] = 1;
    inc[1] = nn[0]-shape_.GetSize()[0];
    inc[2] = nn[0]*(nn[1]-shape_.GetSize()[1]);
    // loop over field voxels
    for (ii[2] = rr[2]; ii[2]<nn[2]-rr[2]; ++ii[2]) {
        for (ii[1] = rr[1]; ii[1]<nn[1]-rr[1]; ++ii[1]) {
            for (ii[0] = rr[0]; ii[0]<nn[0]-rr[0]; ++ii[0]) {
                std::array<int,NDIM> ii_l;
                std::vector<double> field_crop(m_vol_);
                int idx = ii[0] + nn[0]*(ii[1] + nn[1]*ii[2]);
                double img0 = img ? img[idx] : 0.0;
                double tol = 0.1*img0;
                // inner loop over kernel voxels
                int idx_l = 0;
                int idx_s = 0;
                int idx_g = ii[0]-rr[0] + nn[0]*(ii[1]-rr[1] + nn[1]*(ii[2]-rr[2]));
                for (ii_l[2] = -rr[2]; ii_l[2]<=rr[2]; ++ii_l[2]) {
                    for (ii_l[1] = -rr[1]; ii_l[1]<=rr[1]; ++ii_l[1]) {
                        for (ii_l[0] = -rr[0]; ii_l[0]<=rr[0]; ++ii_l[0]) {
                            if (shape_[idx_l] && (!img || (img && abs(img[idx_g]-img0)<tol))) {
                                // set the field_crop to the field value
                                field_crop[idx_s] = src[idx_g];
                                ++idx_s;
                            }
                            ++idx_l;
                            idx_g += inc[0];
                        }
                        idx_g += inc[1];
                    }
                    idx_g += inc[2];
                }
                // compute the median
                double *ptr_begin = field_crop.data();
                double *ptr_end = ptr_begin+idx_s;
                double *ptr_middle = ptr_begin + (ptr_end-ptr_begin)/2;
                std::nth_element(ptr_begin,ptr_middle,ptr_end);
                dst[idx] = *ptr_middle;
            }
        }
    }
    return EPTlibError::Success;
}
