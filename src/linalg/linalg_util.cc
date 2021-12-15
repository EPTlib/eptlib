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

#include "eptlib/linalg/linalg_util.h"

#include <cmath>
#include <complex>
#include <limits>

using namespace eptlib;
using namespace eptlib::linalg;

// Compute the quadratic norm of a vector
real_t eptlib::linalg::Norm2(const real_t *x,const size_t n) {
    real_t eta = MaxAbs(x,n);
    if (eta==0.0) {
        return 0.0;
    }
    real_t sigma = 0;
    for (int i = 0; i<n; ++i) {
        // check to avoid summing values in the underflow region
        if (std::abs(x[i])>=std::sqrt(std::numeric_limits<real_t>::epsilon())*eta) {
            // normalise the vector components to avoid overflow
            sigma += (x[i]/eta)*(x[i]/eta);
        }
    }
    // correct the result with the factor `eta'
    sigma = std::sqrt(sigma)*eta;
    return sigma;
}

// Compute the dot product between two vectors
template <typename NumType>
NumType eptlib::linalg::Dot(const NumType *x,const real_t *y,const size_t n) {
    NumType tau = 0.0;
    for (int i = 0; i<n; ++i) {
        tau += x[i]*y[i];
    }
    return tau;
}

// Get the maximum absolute value of a vector
real_t eptlib::linalg::MaxAbs(const real_t *x,const size_t n) {
    real_t eta = 0;
    for (int i = 0; i<n; ++i) {
        if (std::abs(x[i])>eta) {
            eta = std::abs(x[i]);
        }
    }
    return eta;
}

// Solve a square diagonal system
template <typename NumType>
void eptlib::linalg::SolveDiag(NumType *x,const MatrixReal &A,const NumType *b,const size_t n) {
    for (int i = 0; i<n; ++i) {
        x[i] = b[i]/A[i][i];
    }
    return;
}

// Solve a square upper triangular system
template <typename NumType>
void eptlib::linalg::SolveTriU(NumType *x,const MatrixReal &A,const NumType *b,const size_t n) {
    x[n-1] = b[n-1]/A[n-1][n-1];
    for (int i = n-2; i>=0; --i) {
        NumType s = 0;
        for (int j = i+1; j<n; ++j) {
            s += x[j]*A[j][i];
        }
        x[i] = (b[i]-s)/A[i][i];
    }
    return;
}

// Solve a square lower triangular system
template <typename NumType>
void eptlib::linalg::SolveTriL(NumType *x,const MatrixReal &A,const NumType *b,const size_t n) {
    x[0] = b[0]/A[0][0];
    for (int i = 1; i<n; ++i) {
        NumType s = 0;
        for (int j = 0; j<i; ++j) {
            s += x[j]*A[j][i];
        }
        x[i] = (b[i]-s)/A[i][i];
    }
    return;
}

template real_t eptlib::linalg::Dot<real_t>(const real_t *x,const real_t *y,const size_t n);
template void eptlib::linalg::SolveDiag<real_t>(real_t *x,const MatrixReal &A,const real_t *b,const size_t n);
template void eptlib::linalg::SolveTriU<real_t>(real_t *x,const MatrixReal &A,const real_t *b,const size_t n);
template void eptlib::linalg::SolveTriL<real_t>(real_t *x,const MatrixReal &A,const real_t *b,const size_t n);

template std::complex<real_t> eptlib::linalg::Dot<std::complex<real_t> >(const std::complex<real_t> *x,const real_t *y,const size_t n);
template void eptlib::linalg::SolveDiag<std::complex<real_t> >(std::complex<real_t> *x,const MatrixReal &A,const std::complex<real_t> *b,const size_t n);
template void eptlib::linalg::SolveTriU<std::complex<real_t> >(std::complex<real_t> *x,const MatrixReal &A,const std::complex<real_t> *b,const size_t n);
template void eptlib::linalg::SolveTriL<std::complex<real_t> >(std::complex<real_t> *x,const MatrixReal &A,const std::complex<real_t> *b,const size_t n);
