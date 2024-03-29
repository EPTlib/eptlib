/*****************************************************************************
*
*     Program: EPTlib
*     Author: Alessandro Arduino <a.arduino@inrim.it>
*
*  MIT License
*
*  Copyright (c) 2020-2023  Alessandro Arduino
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

#include "eptlib/linalg/qr.h"

#include <algorithm>
#include <cmath>

// Perform the QR decomposition with column pivoting
std::tuple<eptlib::linalg::Matrix<double>, std::vector<size_t> > eptlib::linalg::QRDecomposition(const eptlib::linalg::Matrix<double> &A) {
    using std::swap;
    const size_t n_row = A.GetNRow();
    const size_t n_col = A.GetNCol();
    const size_t rank_max = std::min(n_row, n_col);
    // initialise matrix QR
    eptlib::linalg::Matrix<double> QR(n_row+2, n_col);
    for (size_t col = 0; col < n_col; ++col) {
        std::copy(A.begin(col), A.end(col), QR.begin(col)+1);
    }
    // initialise the column norms
    std::vector<double> column_norm(n_col);
    for (size_t col = 0; col < n_col; ++col) {
        column_norm[col] = eptlib::linalg::Norm2(A.begin(col), A.end(col));
        column_norm[col] *= column_norm[col];
    }
    // initialise the permutation
    std::vector<size_t> p(n_col);
    std::iota(p.begin(), p.end(), 0);
    // perform the decomposition
    double tol;
    for (size_t col = 0; col < rank_max; ++col) {
        // permute the columns
        auto max_norm = std::max_element(column_norm.begin()+col, column_norm.end());
        auto max_norm_col = std::distance(column_norm.begin(), max_norm);
        if (max_norm_col != col) {
            std::swap_ranges(QR.begin(col), QR.end(col), QR.begin(max_norm_col));
            swap(column_norm[col], column_norm[max_norm_col]);
            swap(p[col], p[max_norm_col]);
        }
        // check if the matrix rank has been reached
        if (col == 0) {
            const double tmp = std::sqrt(column_norm[col]);
            tol = (std::nextafter(tmp, +INFINITY) - tmp) * std::max(n_row, n_col);
            tol *= tol;
        } else {
            if (column_norm[col] < tol) {
                break;
            }
        }
        // compute and apply the reflector
        HouseholderReflector(QR.begin(col)+col+1, QR.end(col));
        std::copy(QR.begin(col)+1, QR.begin(col)+col+1, QR.begin(col));
        QR(col, col) = -QR(n_row+1, col) / QR(col+1, col);
        for (size_t col2 = col+1; col2 < n_col; ++col2) {
            HouseholderLeft(QR.begin(col2)+col+1, QR.end(col2)-1, QR.begin(col)+col+1, QR(n_row+1,col));
        }
        // update the column_norm vector
        for (size_t col2 = col+1; col2 < n_col; ++col2) {
            column_norm[col2] -= QR(col,col2)*QR(col,col2);
        }
    }
    return {QR, p};
}

// Compute the numerical rank of a matrix given its rank-revealing QR decomposition
size_t eptlib::linalg::QRGetRank(const eptlib::linalg::Matrix<double> &QR) {
    const size_t n_row = QR.GetNRow()-2;
    const size_t n_col = QR.GetNCol();
    const size_t max_rank = std::min(n_row, n_col);
    const double tol = (std::nextafter(QR(0,0), +INFINITY) - QR(0,0)) * std::max(n_row, n_col);
    for (size_t rank = 0; rank < max_rank; ++rank) {
        if (std::abs(QR(rank,rank)) < tol) {
            return rank;
        }
    }
    return max_rank;
}
