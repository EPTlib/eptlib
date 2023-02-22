/*****************************************************************************
*
*     Program: EPTlib
*     Author: Alessandro Arduino <a.arduino@inrim.it>
*
*  MIT License
*
*  Copyright (c) 2020-2023  Alessandro Arduino
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

// Shape default constructor
eptlib::Shape::
Shape() :
    Shape(0,0,0) {
    return;
}

// Shape constructor
eptlib::Shape::
Shape(const size_t n0, const size_t n1, const size_t n2) :
    data_(n0*n1*n2) {
    nn_[0] = n0;
    nn_[1] = n1;
    nn_[2] = n2;
    return;
}

eptlib::Shape::
~Shape() {
    return;
}

// Shape operator+,-,& overload
eptlib::Shape& eptlib::Shape::
operator+=(const eptlib::Shape &rhs) {
    this->data_ |= rhs.data_;
    return *this;
}

eptlib::Shape& eptlib::Shape::
operator-=(const eptlib::Shape &rhs) {
    this->data_ -= rhs.data_;
    return *this;
}

eptlib::Shape& eptlib::Shape::
operator&=(const eptlib::Shape &rhs) {
    this->data_ &= rhs.data_;
    return *this;
}

// Shape padding
void eptlib::Shape::
Pad(const size_t left0, const size_t left1, const size_t left2,
    const size_t right0, const size_t right1, const size_t right2) {
    std::array<size_t,N_DIM> nn_old = nn_;
    boost::dynamic_bitset<> data_old = data_;
    // update shape size
    nn_[0] += left0+right0;
    nn_[1] += left1+right1;
    nn_[2] += left2+right2;
    // update shape data
    data_.resize(nn_[0]*nn_[1]*nn_[2]);
    data_.reset();
    for (size_t i2 = 0; i2<nn_old[2]; ++i2) {
        for (size_t i1 = 0; i1<nn_old[1]; ++i1) {
            for (size_t i0 = 0; i0<nn_old[0]; ++i0) {
                size_t idx = IJKToIdx(i0+left0,i1+left1,i2+left2, nn_[0],nn_[1]);
                size_t idx_old = IJKToIdx(i0,i1,i2, nn_old[0],nn_old[1]);
                data_[idx] = data_old[idx_old];
            }
        }
    }
    return;
}

// Shape shrinking
void eptlib::Shape::
Shrink(const size_t left0, const size_t left1, const size_t left2,
    const size_t right0, const size_t right1, const size_t right2) {
    std::array<size_t,N_DIM> nn_old = nn_;
    boost::dynamic_bitset<> data_old = data_;
    // update shape size
    nn_[0] -= left0+right0;
    nn_[1] -= left1+right1;
    nn_[2] -= left2+right2;
    // update shape data
    data_.resize(nn_[0]*nn_[1]*nn_[2]);
    data_.reset();
    for (size_t i2 = 0; i2<nn_[2]; ++i2) {
        for (size_t i1 = 0; i1<nn_[1]; ++i1) {
            for (size_t i0 = 0; i0<nn_[0]; ++i0) {
                size_t idx = IJKToIdx(i0,i1,i2, nn_[0],nn_[1]);
                size_t idx_old = IJKToIdx(i0+left0,i1+left1,i2+left2, nn_old[0],nn_old[1]);
                data_[idx] = data_old[idx_old];
            }
        }
    }
    return;
}

// Collection of shapes
namespace eptlib {

namespace shapes {

    // Cuboid
    Shape Cuboid(const size_t n0, const size_t n1, const size_t n2) {
        Shape cuboid(n0,n1,n2);
        cuboid.GetData().set();
        return cuboid;
    }

    // CuboidR
    Shape CuboidR(const size_t r0, const size_t r1, const size_t r2) {
        Shape cuboid(2*r0+1,2*r1+1,2*r2+1);
        cuboid.GetData().set();
        return cuboid;
    }

    // Ellipsoid
    Shape Ellipsoid(const size_t r0, const size_t r1, const size_t r2) {
        size_t n0 = r0*2+1;
        size_t n1 = r1*2+1;
        size_t n2 = r2*2+1;
        Shape ellipsoid(n0,n1,n2);
        int idx = 0;
        for (size_t i2 = 0; i2<n2; ++i2) {
            for (size_t i1 = 0; i1<n1; ++i1) {
                for (size_t i0 = 0; i0<n0; ++i0) {
                    double x0 = static_cast<double>(i0)-r0;
                    double x1 = static_cast<double>(i1)-r1;
                    double x2 = static_cast<double>(i2)-r2;
                    double rho = 0.0;
                    rho += r0>0 ? x0*x0/r0/r0 : 0.0;
                    rho += r1>0 ? x1*x1/r1/r1 : 0.0;
                    rho += r2>0 ? x2*x2/r2/r2 : 0.0;
                    if (rho<=1.0) {
                        ellipsoid.GetData().set(idx);
                    }
                    ++idx;
                }
            }
        }
        return ellipsoid;
    }

    // Cross
    Shape Cross(const size_t r0, const size_t r1, const size_t r2) {
        size_t n0 = r0*2+1;
        size_t n1 = r1*2+1;
        size_t n2 = r2*2+1;
        Shape cross(n0,n1,n2);
        for (size_t i = 0; i<n0; ++i) {
            cross.GetData().set(IJKToIdx(i,r1,r2, n0,n1));
        }
        for (size_t j = 0; j<n1; ++j) {
            cross.GetData().set(IJKToIdx(r0,j,r2, n0,n1));
        }
        for (size_t k = 0; k<n2; ++k) {
            cross.GetData().set(IJKToIdx(r0,r1,k, n0,n1));
        }
        return cross;
    }

}  // shapes

}  // eptlib
