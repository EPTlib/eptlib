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

#ifndef LINALG_UTIL_H_
#define LINALG_UTIL_H_

#include <vector>

#include "eptlib/util.h"

namespace eptlib {

namespace linalg {

    using MatrixReal = std::vector<std::vector<real_t> >;

    /**
     * @brief Compute the quadratic norm of a vector.
     * 
     * @param x vector.
     * @param n vector size.
     * @return the quadratic norm of vector `x'.
     * 
     * The function is robust with respect to possible overflow.
     */
    real_t Norm2(const real_t *x,const size_t n);

    /**
     * @brief Compute the inner product between two vectors.
     * 
     * @param x vector 1.
     * @param y vector 2.
     * @param n vector size.
     * @return the inner product between `x' and `y'.
     */
    real_t Dot(const real_t *x,const real_t *y,const size_t n);

    /**
     * @brief Get the maximum absolute value of a vector.
     * 
     * @param x vector.
     * @param n vector size.
     * @return the maximum absolute value of vector `x'.
     */
    real_t MaxAbs(const real_t *x,const size_t n);

    /**
     * @brief Solve a square diagonal system.
     * 
     * @param x vector of the unknowns.
     * @param A square diagonal matrix of the coefficients.
     * @param b vector of the forcing terms.
     * @param n vector size.
     */
    void SolveDiag(real_t *x,const MatrixReal &A,const real_t *b,const size_t n);

    /**
     * @brief Solve a square upper triangular system.
     * 
     * @param x vector of the unknowns.
     * @param A square upper triangular matrix of the coefficients.
     * @param b vector of the forcing terms.
     * @param n vector size.
     */
    void SolveTriU(real_t *x,const MatrixReal &A,const real_t *b,const size_t n);
    
    /**
     * @brief Solve a square lower triangular system.
     * 
     * @param x vector of the unknowns.
     * @param A square lower triangular matrix of the coefficients.
     * @param b vector of the forcing terms.
     * @param n vector size.
     */
    void SolveTriL(real_t *x,const MatrixReal &A,const real_t *b,const size_t n);

}  // namespace linalg

}  // namespace eptlib

#endif  // LINALG_UTIL_H_
