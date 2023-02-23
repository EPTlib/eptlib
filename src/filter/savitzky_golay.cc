/*****************************************************************************
*
*     Program: EPTlib
*     Author: Alessandro Arduino <a.arduino@inrim.it>
*
*  MIT License
*
*  Copyright (c) 2023  Alessandro Arduino
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

#include "eptlib/filter/savitzky_golay.h"

#include "eptlib/linalg/linalg_householder.h"
#include "eptlib/linalg/linalg_qr.h"
#include "eptlib/linalg/linalg_util.h"

namespace {

    // Compute the maximum number of monomials in a three-dimensional polynomial
    // of given maximum degree.
    size_t GetNumberOfMonomials(const size_t degree) {
        size_t n_monomials = 0;
        for (size_t deg = 0; deg<degree; ++deg) {
            n_monomials += (deg+2)*(deg+1)/2;
        }
        return n_monomials;
    }

}  //

// SavitzkyGolay constructor
eptlib::SavitzkyGolay::
SavitzkyGolay(const double d0, const double d1, const double d2,
    const eptlib::Shape &window, const size_t degree) :
    window_(window),
    lapl_(0) {
    double x0[3] {window.GetSize(0)/2*d0, window.GetSize(1)/2*d1, window.GetSize(2)/2*d2};
    // initialise the design matrix
    size_t n_row = window_.GetVolume();
    size_t n_col = ::GetNumberOfMonomials(degree);
    eptlib::linalg::MatrixReal F(n_col, std::vector<double>(n_row));
    // fill the design matrix (1, x^2, y^2, z^2, x, y, z, x*y, x*z, y*z, x^3, x^2*y, x*y^2, ...)
    size_t row = 0;
    for (size_t i2 = 0; i2<window.GetSize(2); ++i2) {
        for (size_t i1 = 0; i1<window.GetSize(1); ++i1) {
            for (size_t i0 = 0; i0<window.GetSize(0); ++i0) {
                if (window_(i0, i1, i2)) {
                    double dx[3] {i0*d0-x0[0], i1*d1-x0[1], i2*d2-x0[2]};
                    // even monomials with deg<=2 (1, x^2, y^2, z^2)
                    F[0][row] = 1.0;
                    for (size_t d = 0; d<3; ++d) {
                        F[1+d][row] = dx[d]*dx[d]/2.0;
                    }
                    // odd monomials with deg<=2 (x, y, z, xy, xz, yz)
                    {
                        size_t col = 0;
                        for (size_t d = 0; d<3; ++d) {
                            F[4+d][row] = dx[d];
                            for (size_t d2 = d+1; d2<3; ++d2) {
                                F[7+col][row] = dx[d]*dx[d2];
                                ++col;
                            }
                        }
                    }
                    // monomials with deg>2 (x^j[0] * y^j[1] * z^j[2])
                    {
                        size_t col = 0;
                        for (size_t deg = 3; deg<=degree; ++deg) {
                            size_t j[3];
                            for (j[2] = 0; j[2]<=deg; ++j[2]) {
                                for (j[1] = 0; j[1]<=deg-j[2]; ++j[1]) {
                                    j[0] = deg-j[1]-j[2];
                                    F[11+col][row] = 1;
                                    for (size_t d = 0; d<N_DIM; ++d) {
                                        for (size_t counter = 0; counter<j[d]; ++counter) {
                                            F[11+col][row] *= dx[d];
                                        }
                                    }
                                    ++col;
                                }
                            }
                        }
                    }
                    ++row;
                }
            }
        }
    }
    // solve the linear system
    eptlib::linalg::MatrixReal QR;
    eptlib::linalg::HouseholderQR(&QR, F ,n_row ,n_col);
    eptlib::linalg::MatrixReal A(n_row, std::vector<double>(n_col));
    for (int row = 0; row<n_row; ++row) {
        std::vector<double> b(n_row,0.0);
        b[row] = 1.0;
        eptlib::linalg::QRSolve(A[row].data(), QR, b.data(), n_row, n_col);
    }
    // compute the Laplacian coefficients
    lapl_.resize(n_row);
    for (int row = 0; row<n_row; ++row) {
        lapl_[row] = 2.0*(A[row][1]+A[row][2]+A[row][3]);
    }
    return;
}

// SavitzkyGolay virtual destructor
eptlib::SavitzkyGolay::
~SavitzkyGolay() {
    return;
}
