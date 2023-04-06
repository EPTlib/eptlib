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
eptlib::linalg::Matrix<double> eptlib::polynomial::DesignMatrixWithMonomialsBasis(const std::vector<double> &x,
    const std::vector<double> &y, const std::vector<double> &z, const size_t degree) {
    size_t n_row = x.size();
    size_t n_col = eptlib::polynomial::GetNumberOfMonomials(degree);
    eptlib::linalg::Matrix<double> F(n_row, n_col);
    for (size_t row = 0; row < n_row; ++row) {
        size_t col = 0; 
        for (size_t n = 0; n <= degree; ++n) {
            std::vector<double> monomials = eptlib::polynomial::EvaluateMonomials(x[row], y[row], z[row], n);
            for (double monomial : monomials) {
                F(row, col) = monomial;
                ++col;
            }
        }
    }
    return F;
}
