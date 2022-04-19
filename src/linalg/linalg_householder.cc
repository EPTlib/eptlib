/*****************************************************************************
*
*     Program: EPTlib
*     Author: Alessandro Arduino <a.arduino@inrim.it>
*
*  MIT License
*
*  Copyright (c) 2020-2022  Alessandro Arduino
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

#include "eptlib/linalg/linalg_householder.h"

#include <cmath>
#include <complex>
#include <limits>

#include "eptlib/linalg/linalg_util.h"

using namespace eptlib;
using namespace eptlib::linalg;

// Compute the elementary reflector
void eptlib::linalg::HouseholderReflector(double *x,const size_t n) {
    double sigma = std::copysign(Norm2(x,n),x[0]);
    // sum avoiding loss of significance
    x[0] += sigma;
    // store the factor pi
    x[n] = sigma*x[0];
    return;
}

// Compute the product of the elementary reflector with a vector
template <typename NumType>
void eptlib::linalg::HouseholderLeft(NumType *x,const double *u,const size_t n) {
    NumType tau = Dot(x,u,n)/u[n];
    for (int i = 0; i<n; ++i) {
        x[i] -= tau*u[i];
    }
    return;
}

template void eptlib::linalg::HouseholderLeft<double>(double *x,const double *u,const size_t n);
template void eptlib::linalg::HouseholderLeft<std::complex<double> >(std::complex<double> *x,const double *u,const size_t n);
