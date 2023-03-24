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

#include "gtest/gtest.h"

#include "eptlib/linalg/matrix.h"

#include <algorithm>
#include <complex>
#include <vector>

TEST(LinalgMatrixGTest,SolveDiag) {
    const size_t n_row = 10;
    const size_t n_col = 10;

    eptlib::linalg::Matrix<double> A(n_row, n_col);
    std::vector<double> x_ref(n_col);
    for (size_t col = 0; col < n_col; ++col) {
        for (size_t row = 0; row < n_row; ++row) {
            A(row, col) = col + row*n_col + 1.0;
        }
        x_ref[col] = col + 1.0;
    }
    
    std::vector<double> b(n_row);
    for (size_t row = 0; row < n_row; ++row) {
        b[row] = A(row, row) * x_ref[row];
    }

    std::vector<std::complex<double> > b_c(n_row);
    std::copy(b.begin(), b.end(), b_c.begin());

    std::vector<double> x = eptlib::linalg::SolveDiag(A, b);
    std::vector<std::complex<double> > x_c = eptlib::linalg::SolveDiag(A, b_c);

    for (size_t col = 0; col < n_col; ++col) {
        ASSERT_NEAR(x[col], x_ref[col], 1e-12);
        ASSERT_NEAR(x_c[col].real(), x_ref[col], 1e-12);
        ASSERT_NEAR(x_c[col].imag(), 0.0, 1e-12);
    }
}

TEST(LinalgMatrixGTest,SolveTriU) {
    const size_t n_row = 10;
    const size_t n_col = 10;

    eptlib::linalg::Matrix<double> A(n_row, n_col);
    std::vector<double> x_ref(n_col);
    for (size_t col = 0; col < n_col; ++col) {
        for (size_t row = 0; row < n_row; ++row) {
            A(row, col) = col + row*n_col + 1.0;
        }
        x_ref[col] = col + 1.0;
    }

    std::vector<double> b(n_row);
    for (size_t row = 0; row < n_row; ++row) {
        b[row] = 0.0;
        for (size_t col = row; col < n_col; ++col) {
            b[row] += A(row, col) * x_ref[col];
        }
    }

    std::vector<std::complex<double> > b_c(n_row);
    std::copy(b.begin(), b.end(), b_c.begin());

    std::vector<double> x = eptlib::linalg::SolveTriU(A, b);
    std::vector<std::complex<double> > x_c = eptlib::linalg::SolveTriU(A, b_c);

    for (size_t col = 0; col < n_col; ++col) {
        ASSERT_NEAR(x[col], x_ref[col], 1e-12);
        ASSERT_NEAR(x_c[col].real(), x_ref[col], 1e-12);
        ASSERT_NEAR(x_c[col].imag(), 0.0, 1e-12);
    }
}

TEST(LinalgMatrixGTest,SolveTriL) {
    const size_t n_row = 10;
    const size_t n_col = 10;

    eptlib::linalg::Matrix<double> A(n_row, n_col);
    std::vector<double> x_ref(n_col);
    for (size_t col = 0; col < n_col; ++col) {
        for (size_t row = 0; row < n_row; ++row) {
            A(row, col) = col + row*n_col + 1.0;
        }
        x_ref[col] = col + 1.0;
    }

    std::vector<double> b(n_row);
    for (size_t row = 0; row < n_row; ++row) {
        b[row] = 0.0;
        for (size_t col = 0; col <= row; ++col) {
            b[row] += A(row, col) * x_ref[col];
        }
    }

    std::vector<std::complex<double> > b_c(n_row);
    std::copy(b.begin(), b.end(), b_c.begin());

    std::vector<double> x = eptlib::linalg::SolveTriL(A, b);
    std::vector<std::complex<double> > x_c = eptlib::linalg::SolveTriL(A, b_c);

    for (size_t col = 0; col < n_col; ++col) {
        ASSERT_NEAR(x[col], x_ref[col], 1e-12);
        ASSERT_NEAR(x_c[col].real(), x_ref[col], 1e-12);
        ASSERT_NEAR(x_c[col].imag(), 0.0, 1e-12);
    }
}
