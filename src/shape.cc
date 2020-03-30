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

#include <iostream>

#include <utility>

using namespace eptlib;

// Shape constructor
Shape::
Shape(const std::array<int,NDIM> &nn) :
    nn_(nn),
    n_vox_(std::accumulate(nn.begin(),nn.end(),1,std::multiplies<int>())),
    shape_(n_vox_,false), is_symmetric_(false) {
    return;
}

// Shape getters
const std::array<int,NDIM>& Shape::
GetSize() const {
    return nn_;
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
int Shape::
GetBoxVolume() const {
    return n_vox_;
}
bool Shape::
IsSymmetric() const {
    return is_symmetric_;
}

// Shape operator[] overload
boost::dynamic_bitset<>::reference Shape::
operator[](const std::array<int,NDIM> &ii) {
    int idx = ii[0] + nn_[0]*(ii[1] + nn_[1]*ii[2]);
    return shape_[idx];
}
boost::dynamic_bitset<>::reference Shape::
operator[](const int &idx) {
    return shape_[idx];
}
// Shape operator[] const overload
bool Shape::
operator[](const std::array<int,NDIM> &ii) const {
    int idx = ii[0] + nn_[0]*(ii[1] + nn_[1]*ii[2]);
    return shape_[idx];
}
bool Shape::
operator[](const int &idx) const {
    return shape_[idx];
}

// Shape operator+,-,& overload
Shape& Shape::
operator+=(const Shape &rhs) {
    shape_ |= rhs.shape_;
    if (!(is_symmetric_&&rhs.is_symmetric_)) {
        CheckSymmetry();
    }
    return *this;
}
Shape& Shape::
operator-=(const Shape &rhs) {
    shape_ -= rhs.shape_;
    if (!(is_symmetric_&&rhs.is_symmetric_)) {
        CheckSymmetry();
    }
    return *this;
}
Shape& Shape::
operator&=(const Shape &rhs) {
    shape_ &= rhs.shape_;
    if (!(is_symmetric_&&rhs.is_symmetric_)) {
        CheckSymmetry();
    }
    return *this;
}

// Shape check symmetry
bool Shape::
CheckSymmetry() {
    std::array<int,NDIM> ii;
    std::array<int,NDIM> jj;
    std::array<int,NDIM> mm;
    for (int d = 0; d<NDIM; ++d) {
        std::copy(nn_.begin(),nn_.end(),mm.begin());
        mm[d] = mm[d]/2;
        for (ii[2] = 0; ii[2]<mm[2]; ++ii[2]) {
            for (ii[1] = 0; ii[1]<mm[1]; ++ii[1]) {
                for (ii[0] = 0; ii[0]<mm[0]; ++ii[0]) {
                    std::copy(ii.begin(),ii.end(),jj.begin());
                    jj[d] = nn_[d]-1-ii[d];
                    if ((*this)[ii]!=(*this)[jj]) {
                        is_symmetric_ = false;
                        return is_symmetric_;
                    }
                }
            }
        }
    }
    is_symmetric_ = true;
    return is_symmetric_;
}

// Shape padding
void Shape::
Pad(const std::array<int,NDIM> &l, const std::array<int,NDIM> &r) {
    // initialise the new shape
    std::array<int,NDIM> xnn;
    for (int d = 0; d<NDIM; ++d) {
        xnn[d] = nn_[d]+l[d]+r[d];
    }
    int xn_vox = std::accumulate(xnn.begin(),xnn.end(),1,std::multiplies<int>());
    boost::dynamic_bitset<> xshape(xn_vox,false);
    // fill the new shape
    for (int i2 = 0; i2<nn_[2]; ++i2) {
        for (int i1 = 0; i1<nn_[1]; ++i1) {
            for (int i0 = 0; i0<nn_[0]; ++i0) {
                int idx = i0 + nn_[0]*(i1 + nn_[1]*i2);
                int xidx = i0+l[0] + xnn[0]*(i1+l[1] + xnn[1]*(i2+l[2]));
                xshape[xidx] = shape_[idx];
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
Shrink(const std::array<int,NDIM> &l, const std::array<int,NDIM> &r) {
    // initialise the new shape
    std::array<int,NDIM> xnn;
    for (int d = 0; d<NDIM; ++d) {
        xnn[d] = nn_[d]-l[d]-r[d];
    }
    int xn_vox = std::accumulate(xnn.begin(),xnn.end(),1,std::multiplies<int>());
    boost::dynamic_bitset<> xshape(xn_vox,false);
    // fill the new shape
    for (int xi2 = 0; xi2<xnn[2]; ++xi2) {
        for (int xi1 = 0; xi1<xnn[1]; ++xi1) {
            for (int xi0 = 0; xi0<xnn[0]; ++xi0) {
                int idx = xi0+l[0] + nn_[0]*(xi1+l[1] + nn_[1]*(xi2+l[2]));
                int xidx = xi0 + xnn[0]*(xi1 + xnn[1]*xi2);
                xshape[xidx] = shape_[idx];
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


// Collection of shapes
namespace eptlib::shapes {

    // Cuboid
    Shape Cuboid(const std::array<int,NDIM> &nn) {
        Shape cuboid(nn);
        cuboid.GetShape().set();
        cuboid.CheckSymmetry();
        return cuboid;
    }

    // Ellipsoid
    Shape Ellipsoid(const std::array<int,NDIM> &rr) {
        // compute the dimension of the fitting grid
        std::array<int,NDIM> nn;
        for (int d = 0; d<NDIM; ++d) {
            nn[d] = 2*rr[d]+1;
        }
        // create the ellipsoid
        Shape ellipsoid(nn);
        std::array<int,NDIM> xx;
        int idx = 0;
        for (xx[2] = -rr[2]; xx[2]<=rr[2]; ++xx[2]) {
            for (xx[1] = -rr[1]; xx[1]<=rr[1]; ++xx[1]) {
                for (xx[0] = -rr[0]; xx[0]<=rr[0]; ++xx[0]) {
                    double rho = 0.0;
                    for (int d = 0; d<NDIM; ++d) {
                        rho += static_cast<double>(xx[d])*xx[d]/rr[d]/rr[d];
                    }
                    if (rho<=1.0) {
                        ellipsoid.GetShape().set(idx);
                    }
                    ++idx;
                }
            }
        }
        ellipsoid.CheckSymmetry();
        return ellipsoid;
    }

    // Cross
    Shape Cross(const std::array<int,NDIM> &rr) {
        // compute the dimension of the fitting grid
        std::array<int,NDIM> nn;
        for (int d = 0; d<NDIM; ++d) {
            nn[d] = 2*rr[d]+1;
        }
        int n_vox = std::accumulate(nn.begin(),nn.end(),1,std::multiplies<int>());
        int idx0 = n_vox/2;
        // create the cross
        Shape cross(nn);
        for (int d = 0; d<NDIM; ++d) {
            int step = 1;
            for (int d2 = 0; d2<d; ++d2) {
                step *= nn[d2];
            }
            int idx = idx0-step*rr[d];
            for (int i = 0; i<nn[d]; ++i) {
                cross.GetShape().set(idx);
                idx += step;
            }
        }
        cross.CheckSymmetry();
        return cross;
    }

}  // shapes
