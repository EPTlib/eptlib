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

#ifndef LINALG_QR_H_
#define LINALG_QR_H_

#include "eptlib/util.h"
#include "eptlib/linalg/linalg_util.h"

namespace eptlib {

namespace linalg {

    /**
     * @brief Perform the QR decomposition of a vertical full-rank matrix.
     * 
     * @param qr compact form of the QR decomposition (column-wise matrix).
     * @param A vertical full-rank matrix to be decomposed (column-wise matrix).
     * @param m number of rows in the matrix.
     * @param n number of columns in the matrix.
     * 
     * The upper triangular part of `qr' stores the matrix R. The lower triangular
     * part of `qr' stores the elementary reflections describing the application
     * of the transponse of matrix Q.
     */
    void HouseholderQR(MatrixReal *qr,const MatrixReal &A,const size_t m,const size_t n);

    /**
     * @brief Solve a linear system using the QR decomposition.
     * 
     * @tparam NumType numeric typename.
     * 
     * @param x vector of the system solution.
     * @param qr compact form of the QR decomposition (output of HouseholderQR).
     * @param b vector of the forcing term.
     * @param m number of rows in the matrix.
     * @param n numer of columns in the matrix.
     * @return square of quadratic norm of the residual.
     */
    template <typename NumType>
    double QRSolve(NumType *x,const MatrixReal &qr,const NumType *b,const size_t m,const size_t n);

}  // namespace linalg

}  // namespace eptlib

#endif  // LINALG_QR_H_
