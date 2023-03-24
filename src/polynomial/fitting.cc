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

#include "eptlib/polynomial/fitting.h"

// Get the number of three-dimensional monomials of given maximum degree 
size_t eptlib::polynomial::GetNumberOfMonomials(const size_t degree) {
    size_t n_monomials = 0;
    for (size_t n = 0; n<=degree; ++n) {
        n_monomials += (n+2)*(n+1)/2;
    }
    return n_monomials;
}

// Evaluate all the monomials of given degree
std::vector<double> eptlib::polynomial::EvaluateMonomials(const double x, const double y, const double z, const size_t n) {
    std::vector<double> monomials;
    for (size_t k = 0; k<=n; ++k) {
        for (size_t j = 0; j<=n-k; ++j) {
            size_t i = n-j-k;
            double monomial = 1.0;
            for (size_t counter = 0; counter<i; ++counter) {
                monomial *= x;
            }
            for (size_t counter = 0; counter<j; ++counter) {
                monomial *= y;
            }
            for (size_t counter = 0; counter<k; ++counter) {
                monomial *= z;
            }
            monomials.push_back(monomial);
        }
    }
    return monomials;
}

// Fill the design matrix for polynomial fitting with a basis of monomials
eptlib::linalg::MatrixReal eptlib::polynomial::DesignMatrixWithMonomialsBasis(const double d0, const double d1, const double d2,
    const eptlib::Shape &window, const size_t degree) {
    // initialize the matrix
    size_t n_row = window.GetVolume();
    size_t n_col = eptlib::polynomial::GetNumberOfMonomials(degree);
    eptlib::linalg::MatrixReal F(n_col, std::vector<double>(n_row));
    // fill the matrix
    double x0 = window.GetSize(0)/2*d0;
    double y0 = window.GetSize(1)/2*d1;
    double z0 = window.GetSize(2)/2*d2;
    size_t row = 0;
    for (size_t i2 = 0; i2<window.GetSize(2); ++i2) {
        for (size_t i1 = 0; i1<window.GetSize(1); ++i1) {
            for (size_t i0 = 0; i0<window.GetSize(0); ++i0) {
                if (window(i0, i1, i2)) {
                    double dx = i0*d0-x0;
                    double dy = i1*d1-y0;
                    double dz = i2*d2-z0;
                    size_t col = 0;
                    for (size_t n = 0; n<=degree; ++n) {
                        std::vector<double> monomials = eptlib::polynomial::EvaluateMonomials(dx,dy,dz, n);
                        for (double monomial : monomials) {
                            F[col][row] = monomial;
                            ++col;
                        }
                    }
                    ++row;
                }
            }
        }
    }
    return F;
}

// Permute the columns of a matrix to have all the null columns at the end
std::tuple<std::vector<size_t>, size_t> eptlib::polynomial::PermuteColumns(eptlib::linalg::MatrixReal *A, const size_t n_row, const size_t n_col) {
    std::vector<size_t> p(n_col);
    std::iota(p.begin(), p.end(), 0);
    size_t n_col_shrinked = n_col;
    for (size_t col = 0; col<n_col_shrinked; ++col) {
        if (eptlib::linalg::Norm2(A->at(col).data(), n_row) == 0.0) {
            --n_col_shrinked;
            std::swap(A->at(col), A->at(n_col_shrinked));
            std::swap(p[col], p[n_col_shrinked]);
            --col;
        }
    }
    return {p, n_col_shrinked};
}
