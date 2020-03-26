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

#include <array>
#include <numeric>

#include <boost/dynamic_bitset.hpp>
#include <boost/operators.hpp>

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
         * @param nn number of voxels in each direction.
         */
        Shape(const std::array<int,NDIM> &nn);

        /**
         * Get a read-only reference to the number of voxels in each direction.
         * 
         * @return a read-only reference to the number of voxels.
         */
        const std::array<int,NDIM>& GetSize() const;
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
         * Get the number of voxels inside the shape.
         * 
         * @return the number of voxels inside the shape.
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
         * @param ii multi-index of the element.
         * 
         * @return a reference to the ii-th element in the grid.
         */
        boost::dynamic_bitset<>::reference operator[](const std::array<int,NDIM> &ii);
        /**
         * Get a reference to the idx-th element in the grid.
         * 
         * @param idx single index of the element.
         * 
         * @return a reference to the idx-th element in the grid.
         */
        boost::dynamic_bitset<>::reference operator[](const int &idx);
        /**
         * Get a copy of the ii-th element in the grid.
         * 
         * @param ii multi-index of the element.
         * 
         * @return a copy of the ii-th element in the grid.
         */
        bool operator[](const std::array<int,NDIM> &ii) const;
        /**
         * Get a copy to the idx-th element in the grid.
         * 
         * @param idx single index of the element.
         * 
         * @return a copy to the idx-th element in the grid.
         */
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
         * @param l number of layers before the existing grid.
         * @param r number of layers after the existing grid.
         * 
         * Updates the symmetry of the shape.
         */
        void Pad(const std::array<int,NDIM> &l, const std::array<int,NDIM> &r);
        /**
         * Remove external layers from the grid in a certain direction.
         * 
         * @param l number of layers removed before the new grid.
         * @param r number of layers removed after the new grid.
         * 
         * Updates the symmetry of the shape.
         */
        void Shrink(const std::array<int,NDIM> &l, const std::array<int,NDIM> &r);
    private:
        /// Number of voxels in each direction.
        std::array<int,NDIM> nn_;
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
     * @param nn number of voxels in each direction.
     * 
     * @return a cuboid shape.
     */
    Shape Cuboid(const std::array<int,NDIM> &nn);
    /**
     * Create an ellipsoid shape fitted within a grid.
     * 
     * @param rr semi-axes (in voxels) of the ellipsoid.
     */
    Shape Ellipsoid(const std::array<int,NDIM> &rr);
    /**
     * Create a cross shape fitted within a grid.
     * 
     * @param rr semi-length (in voxels) of the cross lines.
     */
    Shape Cross(const std::array<int,NDIM> &rr);

}  // shapes

}  // eptlib

#endif  // EPTLIB_SHAPE_H_
