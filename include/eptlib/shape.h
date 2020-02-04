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

#ifndef EPTLIB_SHAPE_H_
#define EPTLIB_SHAPE_H_

#include <algorithm>
#include <cassert>
#include <functional>
#include <numeric>

#include <boost/dynamic_bitset.hpp>
#include <boost/operators.hpp>

#include "eptlib/config.h"
#include "eptlib/util.h"

namespace eptlib {

/**
 * Class describing a generic shape in the space. Can be used in the
 * definition of domains, masks, kernels, patches, and so on.
 */
class Shape : public boost::addable<Shape>, boost::subtractable<Shape>,
    boost::andable<Shape> {
    public:
        /**
         * Constructor.
         * 
         * @tparam T iterator typename.
         * 
         * @param nn number of voxels in each direction.
         * 
         * Argument `nn' must have the methods `begin', `end', `size' and
         * `operator[]' (any sequence container from STL works fine).
         */
        template <typename T>
        Shape(const T &nn);

        /**
         * Get the number of spatial dimensions.
         * 
         * @return the number of spatial dimensions.
         */
        int GetNDim() const;
        /**
         * Get the number of grid voxels in a direction.
         * 
         * @param d Cartesian direction id.
         * 
         * @return the number of grid voxels in direction `d'.
         */
        int GetSize(const int d) const;
        /**
         * Get the total number of grid voxels.
         * 
         * @return the total number of grid voxels.
         */
        int GetSize() const;
        /**
         * Get a reference to the shape descriptor.
         * 
         * @return a reference to the shape descriptor.
         */
        boost::dynamic_bitset<>& GetShape();
        /**
         * Get a read-only reference to the shape descriptor.
         * 
         * @return a read-only reference to the shape descriptor.
         */
        const boost::dynamic_bitset<>& GetShape() const;
        /**
         * Get the number of voxels in the shape.
         * 
         * @return the number of voxels in the shape.
         */
        int GetVolume() const;
        /**
         * Get if the shape has mirror symmetries.
         * 
         * @return if the shape has mirror symmetries.
         */
        bool IsSymmetric() const;

        /**
         * Get a reference to the ii-th element in the grid.
         * 
         * @tparam T iterator typename.
         * 
         * @param ii multi-index of the element.
         * 
         * @return a reference to the ii-th element in the grid.
         */
        template <typename T>
        boost::dynamic_bitset<>::reference operator[](const T &ii);
        /**
         * Get a reference to the idx-th element in the grid.
         * 
         * @param idx single index of the element.
         * 
         * @return a reference to the idx-th element in the grid.
         */
        template <>
        boost::dynamic_bitset<>::reference operator[](const int &idx);
        /**
         * Get a copy of the ii-th element in the grid.
         * 
         * @tparam T iterator typename.
         * 
         * @param ii multi-index of the element.
         * 
         * @return a copy of the ii-th element in the grid.
         */
        template <typename T>
        bool operator[](const T &ii) const;
        /**
         * Get a copy to the idx-th element in the grid.
         * 
         * @param idx single index of the element.
         * 
         * @return a copy to the idx-th element in the grid.
         */
        template <>
        bool operator[](const int &idx) const;

        /**
         * Set union of two shapes.
         * 
         * @param rhs shape to be united.
         * 
         * @return a reference to _this_.
         */
        Shape& operator+=(const Shape &rhs);
        /**
         * Set difference of two shapes.
         * 
         * @param rhs shape to be subtracted.
         * 
         * @return a reference to _this_.
         */
        Shape& operator-=(const Shape &rhs);
        /**
         * Set intersection of two shapes.
         * 
         * @param rhs shape to be intersected.
         * 
         * @return a reference to _this_.
         */
        Shape& operator&=(const Shape &rhs);

        /**
         * Check and set if the shape has mirror symmetries.
         * 
         * @return if the shape has mirror symmetries.
         */
        bool CheckSymmetry();

        /**
         * Add void layers around the grid in a certain direction.
         * 
         * @param d Cartesian direction along which add the layers.
         * @param left number of layers before the existing grid.
         * @param right number of layers after the existing grid.
         * 
         * Updates the symmetry of the shape.
         */
        void Pad(const int d, const int left, const int right);
        /**
         * Remove external layers from the grid in a certain direction.
         * 
         * @param d Cartesian direction along which remove the layers.
         * @param left number of layers removed before the new grid.
         * @param right number of layers removed after the new grid.
         * 
         * Updates the symmetry of the shape.
         */
        void Shrink(const int d, const int left, const int right);
    private:
        /// Number of spatial dimensions.
        int n_dim_;
        /// Number of voxels in each direction.
        std::vector<int> nn_;
        /// Total number of voxels.
        int n_vox_;
        /// "Black and white" shape descriptor.
        boost::dynamic_bitset<> shape_;
        /// Has the shape got mirror symmetries?
        bool is_symmetric_;
};

/// Collection of methods to generate elementary shapes.
namespace shapes {

/**
 * Create a cuboid shape that fill the input grid.
 * 
 * @tparam T iterator typename.
 * 
 * @param nn number of voxels in each direction.
 * 
 * @return a cuboid shape.
 */
template <typename T>
Shape Cuboid(const T &nn);
/**
 * Create an ellipsoid shape fitted within a grid.
 * 
 * @tparam T iterator typename
 * 
 * @param rr semi-axes (in voxels) of the ellipsoid.
 */
template <typename T>
Shape Ellipsoid(const T &rr);
/**
 * Create a cross shape fitted within a grid.
 * 
 * @tparam T iterator typename
 * 
 * @param rr semi-length (in voxels) of the cross lines.
 */
template <typename T>
Shape Cross(const T &rr);

}  // shapes


// ---------------------------------------------------------------------------
// -------------------------  Implementation detail  -------------------------
// ---------------------------------------------------------------------------

// Shape constructor
template <typename T>
Shape::
Shape(const T &nn) :
    n_dim_(static_cast<int>(nn.size())) {
    // copy the content of nn
    nn_.resize(n_dim_);
    std::copy(nn.begin(),nn.end(),nn_.begin());
    n_vox_ = std::accumulate(nn_.begin(),nn_.end(),1,std::multiplies<int>());
    // initialise a void shape
    shape_.resize(n_vox_,false);
    // initialise the data flags
    is_symmetric_ = false;
    return;
}
// Shape operator[] overload
template <typename T>
boost::dynamic_bitset<>::reference Shape::
operator[](const T &ii) {
    assert(ii.size()==n_dim_);
    int idx = MultiIdxToIdx(ii,nn_);
    return shape_[idx];
}
template <>
boost::dynamic_bitset<>::reference Shape::
operator[](const int &idx) {
    return shape_[idx];
}
// Shape operator[] const overload
template <typename T>
bool Shape::
operator[](const T &ii) const {
    assert(ii.size()==n_dim_);
    int idx = MultiIdxToIdx(ii,nn_);
    return shape_[idx];
}
template<>
bool Shape::
operator[](const int &idx) const {
    return shape_[idx];
}

// Collection of shapes
namespace shapes {

    // Cuboid
    template <typename T>
    Shape Cuboid(const T &nn) {
        Shape cuboid(nn);
        cuboid.GetShape().set();
        cuboid.CheckSymmetry();
        return cuboid;
    }
    // Ellipsoid
    template <typename T>
    Shape Ellipsoid(const T &rr) {
        int n_dim = static_cast<int>(rr.size());
        // compute the dimension of the fitting grid
        std::vector<int> nn(n_dim);
        std::transform(rr.begin(),rr.end(),nn.begin(),
            [](const int &r)->int{return 2*r+1;});
        // compute the coordinates of voxel's barycenter in the grid
        std::vector<std::vector<int>> xx(n_dim);
        for (int d = 0; d<n_dim; ++d) {
            xx[d].resize(nn[d]);
            std::iota(xx[d].begin(),xx[d].end(),-rr[d]);
        }
        // create the ellipsoid
        Shape ellipsoid(nn);
        std::vector<int> ii(n_dim);
        real_t rho;
        for (int idx = 0; idx<ellipsoid.GetSize(); ++idx) {
            IdxToMultiIdx(ii,idx,nn);
            rho = 0.0;
            for (int d = 0; d<n_dim; ++d) {
                real_t x = xx[d][ii[d]];
                real_t r = rr[d];
                rho += x*x/r/r;
            }
            if (rho <= 1.0) {
                ellipsoid.GetShape().set(idx);
            }
        }
        ellipsoid.CheckSymmetry();
        return ellipsoid;
    }
    // Cross
    template <typename T>
    Shape Cross(const T &rr) {
        int n_dim = static_cast<int>(rr.size());
        // compute the dimension of the fitting grid
        std::vector<int> nn(n_dim);
        std::transform(rr.begin(),rr.end(),nn.begin(),
            [](const int &r)->int{return 2*r+1;});
        // create the cross
        Shape cross(nn);
        std::vector<int> ii(n_dim);
        std::copy(rr.begin(),rr.end(),ii.begin());
        for (int d = 0; d<n_dim; ++d) {
            for (int i = 0; i<nn[d]; ++i) {
                ii[d] = i;
                int idx = MultiIdxToIdx(ii,nn);
                cross.GetShape().set(idx);
            }
            ii[d] = rr[d];
        }
        cross.CheckSymmetry();
        return cross;
    }

}  // shapes

}  // eptlib

#endif  // EPTLIB_SHAPE_H_
