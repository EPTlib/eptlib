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

#include "eptlib/linalg/linalg_qr.h"

#include <cstring>

#include "eptlib/linalg/linalg_householder.h"

using namespace eptlib;
using namespace eptlib::linalg;

// Perform the QR decomposition of a vertical full-rank matrix
void eptlib::linalg::HouseholderQR(MatrixReal *qr,const MatrixReal &A,const size_t m,const size_t n) {
    // initialise matrix qr
    qr->resize(n);
    for (int j = 0; j<n; ++j) {
        qr->at(j).resize(m+2);
        std::copy(A.at(j).begin(),A.at(j).end(),qr->at(j).begin()+1);
    }
    // perform the decomposition
    for (int j = 0; j<n; ++j) {
        real_t *u = qr->at(j).data()+1+j;
        int sub_m = m-j;
        HouseholderReflector(u,sub_m);
        for (int i = 0; i<j; ++i) {
            qr->at(j)[i] = qr->at(j)[i+1];
        }
        qr->at(j)[j] = -u[sub_m]/u[0];
        for (int j2 = j+1; j2<n; ++j2) {
            real_t *x = qr->at(j2).data()+1+j;
            HouseholderLeft(x,u,sub_m);
        }
    }
    return;
}

// Solve a linear system using the QR decomposition.
real_t eptlib::linalg::QRSolve(real_t *x,const MatrixReal &qr,const real_t *b,const size_t m,const size_t n) {
    // rotate the forcing term
    real_t *c = new real_t[m];
    std::memcpy(c,b,m*sizeof(real_t));
    for (int j = 0; j<n; ++j) {
        HouseholderLeft(c+j,qr[j].data()+1+j,m-j);
    }
    // solve the system
    SolveTriU(x,qr,c,n);
    // compute the residual
    real_t chi2 = 0.0;
    for (int i = n; i<m; ++i) {
        chi2 += c[i]*c[i];
    }
    // deallocate
    delete[] c;
    return chi2;
}
