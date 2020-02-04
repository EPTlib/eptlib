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

#include "eptlib/shape.h"

#include <utility>

using namespace eptlib;

// Shape getters
int Shape::
GetNDim() const {
    return n_dim_;
}
int Shape::
GetSize(const int d) const {
    return nn_[d];
}
int Shape::
GetSize() const {
    return n_vox_;
}
boost::dynamic_bitset<>& Shape::
GetShape() {
    return shape_;
}
const boost::dynamic_bitset<>& Shape::
GetShape() const {
    return shape_;
}
int Shape::
GetVolume() const {
    return static_cast<int>(shape_.count());
}
bool Shape::
IsSymmetric() const {
    return is_symmetric_;
}

// Shape operator overload
Shape& Shape::
operator+=(const Shape &rhs) {
    // check shape dimensions
    assert(n_dim_==rhs.n_dim_);
    for (int d = 0; d<n_dim_; ++d) {
        assert(nn_[d]==rhs.nn_[d]);
    }
    // unite
    shape_ |= rhs.shape_;
    // check for symmetry
    if (!(is_symmetric_&&rhs.is_symmetric_)) {
        CheckSymmetry();
    }
    return *this;
}
Shape& Shape::
operator-=(const Shape &rhs) {
    // check shape dimensions
    assert(n_dim_==rhs.n_dim_);
    for (int d = 0; d<n_dim_; ++d) {
        assert(nn_[d]==rhs.nn_[d]);
    }
    // subtract
    shape_ -= rhs.shape_;
    // check for symmetry
    if (!(is_symmetric_&&rhs.is_symmetric_)) {
        CheckSymmetry();
    }
    return *this;
}
Shape& Shape::
operator&=(const Shape &rhs) {
    // check shape dimensions
    assert(n_dim_==rhs.n_dim_);
    for (int d = 0; d<n_dim_; ++d) {
        assert(nn_[d]==rhs.nn_[d]);
    }
    // subtract
    shape_ &= rhs.shape_;
    // check for symmetry
    if (!(is_symmetric_&&rhs.is_symmetric_)) {
        CheckSymmetry();
    }
    return *this;
}

// Shape check symmetry
bool Shape::
CheckSymmetry() {
    int step = 1;
    std::vector<int> ii0(n_dim_);
    for (int d = 0; d<n_dim_; ++d) {
        // initialise line length and half-length
        int l = nn_[d];
        int m = l/2;
        for (int i0 = 0; i0<n_vox_; ++i0) {
            // check if i0 starts a line parallel to d
            IdxToMultiIdx(ii0,i0,nn_);
            if (ii0[d]==0) {
                // check the mirror symmetry
                for (int i = 0; i<m; ++i) {
                    if (shape_[i0+i*step]!=shape_[i0+(l-1-i)*step]) {
                        is_symmetric_ = false;
                        return is_symmetric_;
                    }
                }
            }
        }
        // update the step length
        step *= nn_[d];
    }
    // all tests are passed
    is_symmetric_ = true;
    return is_symmetric_;
}

// Shape padding
void Shape::
Pad(const int d, const int l, const int r) {
    std::vector<int> ii0(n_dim_);
    int n = nn_[d];
    // initialise the new shape
    std::vector<int> xnn(n_dim_);
    std::copy(nn_.begin(),nn_.end(),xnn.begin());
    xnn[d] += l+r;
    int xn_vox = std::accumulate(xnn.begin(),xnn.end(),1,std::multiplies<int>());
    boost::dynamic_bitset<> xshape(xn_vox,false);
    // compute the step lengths
    int step = std::accumulate(nn_.begin(),std::next(nn_.begin(),d),1,std::multiplies<int>());
    int xstep = std::accumulate(xnn.begin(),std::next(xnn.begin(),d),1,std::multiplies<int>());
    for (int i0 = 0; i0<n_vox_; ++i0) {
        // check if i0 starts a line parallel to d
        IdxToMultiIdx(ii0,i0,nn_);
        if (ii0[d]==0) {
            int xi0 = MultiIdxToIdx(ii0,xnn);
            // copy the original line in the padded shape
            for (int i = 0; i<n; ++i) {
                xshape[xi0+(l+i)*xstep] = shape_[i0+i*step];
            }
        }
    }
    // update the attributes
    nn_ = xnn;
    n_vox_ = xn_vox;
    shape_ = std::move(xshape);
    CheckSymmetry();
    return;
}

// Shape shrinking
void Shape::
Shrink(const int d, const int l, const int r) {
    std::vector<int> xii0(n_dim_);
    // initialise the new shape
    std::vector<int> xnn(n_dim_);
    std::copy(nn_.begin(),nn_.end(),xnn.begin());
    xnn[d] -= l+r;
    int xn_vox = std::accumulate(xnn.begin(),xnn.end(),1,std::multiplies<int>());
    boost::dynamic_bitset<> xshape(xn_vox,false);
    // compute the step lengths
    int n = xnn[d];
    int step = std::accumulate(nn_.begin(),std::next(nn_.begin(),d),1,std::multiplies<int>());
    int xstep = std::accumulate(xnn.begin(),std::next(xnn.begin(),d),1,std::multiplies<int>());
    for (int xi0 = 0; xi0<xn_vox; ++xi0) {
        // check if i0 starts a line parallel to d
        IdxToMultiIdx(xii0,xi0,xnn);
        if (xii0[d]==0) {
            int i0 = MultiIdxToIdx(xii0,nn_);
            // copy the original line in the shrinked shape
            for (int i = 0; i<n; ++i) {
                xshape[xi0+i*xstep] = shape_[i0+(l+i)*step];
            }
        }
    }
    // update the attributes
    nn_ = xnn;
    n_vox_ = xn_vox;
    shape_ = std::move(xshape);
    CheckSymmetry();
    return;
}
