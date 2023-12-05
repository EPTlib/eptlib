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

#ifndef EPTLIB_IMAGE_H_
#define EPTLIB_IMAGE_H_

#include <array>
#include <vector>

#include "eptlib/util.h"

namespace eptlib {

/**
 * Class describing the input and output images, intended as three-dimensional
 * maps.
 * 
 * @tparam Scalar numerical type of the image data.
 */
template <typename Scalar>
class Image {
    public:
        /**
         * Default constructor.
         */
        Image() :
        Image(0,0,0) {
            return;
        }

        /**
         * 3-D constructor.
         * 
         * @param n0 number of voxels along the x-direction.
         * @param n1 number of voxels along the y-direction.
         * @param n2 number of voxels along the z-direction.
         */
        Image(const size_t n0, const size_t n1, const size_t n2) :
        data_(n0*n1*n2) {
            nn_[0] = n0;
            nn_[1] = n1;
            nn_[2] = n2;
            return;
        }

        /**
         * Destructor.
         */
        virtual ~Image() {
            return;
        }

        /**
         * Get a constant reference to the image dimensions.
         * 
         * @return a constant reference to the dimensions.
         */
        inline const std::array<size_t,N_DIM>& GetSize() const {
            return nn_;
        }

        /**
         * Get the number of voxels along dimension `d'.
         * 
         * @param d dimension of interest.
         * 
         * @return the number of voxels along dimension `d'.
         */
        inline size_t GetSize(const size_t d) const {
            return nn_[d];
        }

        /**
         * Get the number of voxels of the image.
         * 
         * @return the number of voxels of the image.
         */
        inline size_t GetNVox() const {
            return data_.size();
        }

        /**
         * Get a reference to the data vector.
         * 
         * @return a reference to the data.
         */
        inline std::vector<Scalar>& GetData() {
            return data_;
        }

        /**
         * Get a constant reference to the data vector.
         * 
         * @return a constant reference to the data.
         */
        inline const std::vector<Scalar>& GetData() const {
            return data_;
        }

        /**
         * Get a reference to the idx-th voxel in the image.
         * 
         * @param idx single index of the voxel.
         * 
         * @return a reference to the idx-th voxel.
         */
        inline Scalar& operator()(const size_t idx) {
            return data_[idx];
        }

        /**
         * Get a copy of the idx-th voxel in the image.
         * 
         * @param idx single index of the voxel.
         * 
         * @return a copy of the idx-th voxel.
         */
        inline Scalar operator()(const size_t idx) const {
            return data_[idx];
        }

        /**
         * Get a reference to the voxel in the image of indices i,j,k.
         * 
         * @param i index along the x-direction.
         * @param j index along the y-direction.
         * @param k index along the z-direction.
         * 
         * @return a reference to the voxel.
         */
        inline Scalar& operator()(const size_t i, const size_t j, const size_t k) {
            return data_[IJKToIdx(i,j,k,nn_[0],nn_[1])];
        }

        /**
         * Get a copy of the voxel in the image of indices i,j,k.
         * 
         * @param i index along the x-direction.
         * @param j index along the y-direction.
         * @param k index along the z-direction.
         * 
         * @return a copy of the voxel.
         */
        inline Scalar operator()(const size_t i, const size_t j, const size_t k) const {
            return data_[IJKToIdx(i,j,k,nn_[0],nn_[1])];
        }

        /**
         * Proxy for operator().
         * 
         * @param idx single index of the voxel.
         * 
         * @return a reference to the idx-th voxel.
         */
        inline Scalar& At(const size_t idx) {
            return this->operator()(idx);
        }

        /**
         * Proxy for operator().
         * 
         * @param idx single index of the voxel.
         * 
         * @return a copy of the idx-th voxel.
         */
        inline Scalar At(const size_t idx) const {
            return this->operator()(idx);
        }

        /**
         * Proxy for operator().
         * 
         * @param i index along the x-direction.
         * @param j index along the y-direction.
         * @param k index along the z-direction.
         * 
         * @return a reference to the voxel.
         */
        inline Scalar& At(const size_t i, const size_t j, const size_t k) {
            return this->operator()(i,j,k);
        }

        /**
         * Proxy for operator().
         * 
         * @param i index along the x-direction.
         * @param j index along the y-direction.
         * @param k index along the z-direction.
         * 
         * @return a copy of the voxel.
         */
        inline Scalar At(const size_t i, const size_t j, const size_t k) const {
            return this->operator()(i,j,k);
        }
    private:
        /// Number of voxels in each direction.
        std::array<size_t,N_DIM> nn_;
        /// Image data.
        std::vector<Scalar> data_;
};

}  // namespace eptlib

#endif  // EPTLIB_IMAGE_H_
