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
class Shape : boost::addable<Shape>, boost::subtractable<Shape>, boost::andable<Shape> {
    public:
        /**
         * Default constructor.
         */
        Shape();

        /**
         * 3-D constructor.
         * 
         * @param n0 number of voxels along the x-direction.
         * @param n1 number of voxels along the y-direction.
         * @param n2 number of voxels along the z-direction.
         */
        Shape(const size_t n0, const size_t n1, const size_t n2);

        /**
         * Destructor.
         */
        virtual ~Shape();

        /**
         * Get a constant reference to the shape dimensions.
         * 
         * @return a constant reference to the dimensions.
         */
        inline const std::array<size_t,N_DIM>& GetSize() const {
            return nn_;
        }

        /**
         * Get the number of voxels along dimensions `d'.
         * 
         * @param d dimension of interest.
         * 
         * @return the number of voxels along dimension `d'.
         */
        inline size_t GetSize(const size_t d) const {
            return nn_[d];
        }

        /**
         * Get the number of voxels of the shape's bounding box.
         * 
         * @return the number of voxels of the shape's bounding box.
         */
        inline size_t GetNVox() const {
            return data_.size();
        }
        
        /**
         * Get the number of voxels inside the shape.
         * 
         * @return the number of voxels inside the shape.
         */
        inline size_t GetVolume() const {
            return data_.count();
        }

        /**
         * Get a reference to the data bitset.
         * 
         * @return a reference to the data bitset.
         */
        inline boost::dynamic_bitset<>& GetData() {
            return data_;
        }

        /**
         * Get a constant reference to the data bitset.
         * 
         * @return a constant reference to the data bitset.
         */
        inline const boost::dynamic_bitset<>& GetData() const {
            return data_;
        }

        /**
         * Get a reference to the idx-th voxel in the shape.
         * 
         * @param idx single index of the voxel.
         * 
         * @return a reference to the idx-th voxel.
         */
        inline boost::dynamic_bitset<>::reference operator()(const size_t idx) {
            return data_[idx];
        }

        /**
         * Get a copy of the idx-th voxel in the shape.
         * 
         * @param idx single index of the voxel.
         * 
         * @return a copy of the idx-th voxel.
         */
        inline bool operator()(const size_t idx) const {
            return data_[idx];
        };

        /**
         * Get a reference to the voxel in the shape of indices i,j,k.
         * 
         * @param i index along the x-direction.
         * @param j index along the y-direction.
         * @param k index along the z-direction.
         * 
         * @return a reference to the voxel.
         */
        inline boost::dynamic_bitset<>::reference operator()(const size_t i, const size_t j, const size_t k) {
            return data_[IJKToIdx(i,j,k,nn_[0],nn_[1])];
        }

        /**
         * Get a copy of the idx-th voxel in the shape.
         * 
         * @param i index along the x-direction.
         * @param j index along the y-direction.
         * @param k index along the z-direction.
         * 
         * @return a copy of the voxel.
         */
        inline bool operator()(const size_t i, const size_t j, const size_t k) const {
            return data_[IJKToIdx(i,j,k,nn_[0],nn_[1])];
        };

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
         * Add empty layers around the shape.
         * 
         * @param left0 number of layers before the x-direction.
         * @param left1 number of layers before the y-direction.
         * @param left2 number of layers before the z-direction.
         * @param right0 number of layers after the x-direction.
         * @param right1 number of layers after the y-direction.
         * @param right2 number of layers after the z-direction.
         */
        void Pad(const size_t left0, const size_t left1, const size_t left2,
            const size_t right0, const size_t right1, const size_t right2);

        /**
         * Remove external layers from the shape.
         * 
         * @param left0 number of layers before the x-direction.
         * @param left1 number of layers before the y-direction.
         * @param left2 number of layers before the z-direction.
         * @param right0 number of layers after the x-direction.
         * @param right1 number of layers after the y-direction.
         * @param right2 number of layers after the z-direction.
         */
        void Shrink(const size_t left0, const size_t left1, const size_t left2,
            const size_t right0, const size_t right1, const size_t right2);
    private:
        /// Number of voxels in each direction.
        std::array<size_t,N_DIM> nn_;
        /// "Black and white" shape descriptor.
        boost::dynamic_bitset<> data_;
};

/// Collection of methods to generate elementary shapes.
namespace shapes {

    /**
     * Create a cuboid shape that fill the bounding box of sizes `n0', `n1' and `n2'.
     * 
     * @param n0 number of voxels along x-direction.
     * @param n1 number of voxels along y-direction.
     * @param n2 number of voxels along z-direction.
     * 
     * @return a cuboid shape.
     */
    Shape Cuboid(const size_t n0, const size_t n1, const size_t n2);

    /**
     * Create a cuboid shape that fill the bounding box of semiaxes `r0', `r1' and `r2'.
     * 
     * @param r0 semi-axes (in voxel) along x-direction.
     * @param r1 semi-axes (in voxel) along y-direction.
     * @param r2 semi-axes (in voxel) along z-direction.
     * 
     * The bounding box sizes are `r0*2+1', `r1*2+1' and `r2*2+1'.
     * 
     * @return a cuboid shape.
     */
    Shape CuboidR(const size_t r0, const size_t r1, const size_t r2);

    /**
     * Create an ellipsoid shape of semiaxes `r0', `r1' and `r2'.
     * 
     * @param r0 semi-axes (in voxel) along x-direction.
     * @param r1 semi-axes (in voxel) along y-direction.
     * @param r2 semi-axes (in voxel) along z-direction.
     * 
     * The bounding box sizes are `r0*2+1', `r1*2+1' and `r2*2+1'.
     * 
     * @return an ellipsoid shape.
     */
    Shape Ellipsoid(const size_t r0, const size_t r1, const size_t r2);

    /**
     * Create a cross shape of semiaxes `r0', `r1' and `r2'.
     * 
     * @param r0 semi-axes (in voxel) along x-direction.
     * @param r1 semi-axes (in voxel) along y-direction.
     * @param r2 semi-axes (in voxel) along z-direction.
     * 
     * The bounding box sizes are `r0*2+1', `r1*2+1' and `r2*2+1'.
     * 
     * @return a cross shape.
     */
    Shape Cross(const size_t r0, const size_t r1, const size_t r2);

}  // shapes

}  // eptlib

#endif  // EPTLIB_SHAPE_H_
