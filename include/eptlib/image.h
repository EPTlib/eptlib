/*****************************************************************************
*
*     Program: EPTlib
*     Author: Alessandro Arduino <a.arduino@inrim.it>
*
*  MIT License
*
*  Copyright (c) 2020-2021  Alessandro Arduino
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

#include <vector>

#include "eptlib/util.h"

namespace eptlib {

/**
 * Class for N-dimensional images.
 * 
 * @tparam NumType numerical typename of the image data.
 */
template <typename NumType>
class Image {
	public:
		/**
		 * Default constructor.
		 */
		Image();
		/**
		 * 1-D constructor.
		 * 
		 * @param n0 Number of voxels along the line.
		 */
		Image(const int n0);
		/**
		 * 2-D constructor.
		 * 
		 * @param n0 Number of voxels along the x-direction.
		 * @param n1 Number of voxels along the y-direction.
		 */
		Image(const int n0, const int n1);
		/**
		 * 3-D constructor.
		 * 
		 * @param n0 Number of voxels along the x-direction.
		 * @param n1 Number of voxels along the y-direction.
		 * @param n2 Number of voxels along the z-direction.
		 */
		Image(const int n0, const int n1, const int n2);
		/**
		 * N-D constructor.
		 * 
		 * @param nn Number of voxels along each direction.
		 */
		Image(const std::vector<int> &nn);

		/**
		 * Number of dimensions of the image.
		 * 
		 * @return the number of dimensions.
		 */
		int GetNDim() const;
		/**
		 * Get a reference to the image dimensions.
		 * 
		 * @return a reference to the dimensions.
		 */
		std::vector<int>& GetSize();
		/**
		 * Get a constant reference to the image dimensions.
		 * 
		 * @return a constant reference to the dimensions.
		 */
		const std::vector<int>& GetSize() const;
		/**
		 * Number of voxels along dimension d.
		 * 
		 * @param d dimension of interest.
		 * 
		 * @return the number of voxels.
		 */
		int GetSize(const int d) const;
		/**
		 * Number of voxels of the image.
		 * 
		 * @return the number of voxels.
		 */
		size_t GetNVox() const;
		/**
		 * Get a reference to the data vector.
		 * 
		 * @return a reference to the data.
		 */
		std::vector<NumType>& GetData();
		/**
		 * Get a constant reference to the data vector.
		 * 
		 * @return a constant reference to the data.
		 */
		const std::vector<NumType>& GetData() const;
		/**
		 * Get a reference to the idx-th voxel in the image.
		 * 
		 * @param idx single index of the voxel.
		 * 
		 * @return a reference to the idx-th voxel.
		 */
		NumType& operator[](const size_t idx);
		/**
		 * Get a copy of the idx-th voxel in the image.
		 * 
		 * @param idx single index of the voxel.
		 * 
		 * @return a copy of the idx-th voxel.
		 */
		NumType operator[](const size_t idx) const;
	private:
		/// Number of voxels in each direction.
		std::vector<int> nn_;
		/// Image data.
		std::vector<NumType> data_;
};

}  // namespace eptlib

#endif  // EPTLIB_IMAGE_H_
