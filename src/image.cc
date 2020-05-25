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

#include "eptlib/image.h"

using namespace eptlib;

// Image constructors
template <typename NumType>
Image<NumType>::
Image() :
	nn_(0), data_(0) {
	return;
}
template <typename NumType>
Image<NumType>::
Image(const int n0) :
	nn_(1,n0), data_(n0) {
	return;
}
template <typename NumType>
Image<NumType>::
Image(const int n0, const int n1) :
	nn_(2), data_(0) {
	nn_[0] = n0;
	nn_[1] = n1;
	data_.resize(Prod(nn_));
	return;
}
template <typename NumType>
Image<NumType>::
Image(const int n0, const int n1, const int n2) :
	nn_(3), data_(0) {
	nn_[0] = n0;
	nn_[1] = n1;
	nn_[2] = n2;
	data_.resize(Prod(nn_));
	return;
}
template <typename NumType>
Image<NumType>::
Image(const std::vector<int> &nn) :
	nn_(nn), data_(Prod(nn_)) {
	return;
}

// Image getters
template <typename NumType>
int Image<NumType>::
GetNDim() const {
	return static_cast<int>(nn_.size());
}
template <typename NumType>
std::vector<int>& Image<NumType>::
GetSize() {
	return nn_;
}
template <typename NumType>
const std::vector<int>& Image<NumType>::
GetSize() const {
	return nn_;
}
template <typename NumType>
int Image<NumType>::
GetSize(const int d) const {
	return nn_[d];
}
template <typename NumType>
size_t Image<NumType>::
GetNVox() const {
	return data_.size();
}
template <typename NumType>
std::vector<NumType>& Image<NumType>::
GetData() {
	return data_;
}
template <typename NumType>
const std::vector<NumType>& Image<NumType>::
GetData() const {
	return data_;
}
template <typename NumType>
NumType& Image<NumType>::
operator[](const size_t idx) {
	return data_[idx];
}
template <typename NumType>
NumType Image<NumType>::
operator[](const size_t idx) const {
	return data_[idx];
}

template class Image<size_t>;
template class Image<float>;
template class Image<double>;
template class Image<int>;
template class Image<long>;
