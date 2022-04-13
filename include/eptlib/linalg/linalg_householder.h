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

#ifndef LINALG_HOUSEHOLDER_H_
#define LINALG_HOUSEHOLDER_H_

#include "eptlib/util.h"

namespace eptlib {

namespace linalg {

    /**
     * @brief Compute the elementary reflector 
     * 
     * @param x vector for which the reflector is computed.
     * @param n vector size.
     * 
     * The elementary reflector u is stored in place of `x', which should be
     * larger than `n' to allow storing the factor pi.
     */
    void HouseholderReflector(double *x,const size_t n);

    /**
     * @brief Compute the product of the elementary reflector with a vector.
     * 
     * @tparam NumType numeric typename.
     * 
     * @param x vector to which the reflector is applied.
     * @param u elementary reflector.
     * @param n vector size.
     * 
     * The result of the application is stored in place of `x'.
     */
    template <typename NumType>
    void HouseholderLeft(NumType *x,const double *u,const size_t n);

}  // namespace linalg

}  // namespace eptlib

#endif  // LINALG_HOUSEHOLDER_H_
