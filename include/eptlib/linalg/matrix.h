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

#ifndef EPTLIB_LINALG_MATRIX_H_
#define EPTLIB_LINALG_MATRIX_H_

#include <functional>
#include <type_traits>
#include <vector>

namespace eptlib {

namespace linalg {

    /**
     * Class describing a matrix in column-major order.
     * 
     * @tparam Scalar numerical type of the matrix entries. 
     */
    template <typename Scalar>
    class Matrix {
        public:
            /**
             * Construct a matrix with given rows and columns.
             * 
             * @param n_row number of rows.
             * @param n_col number of columns.
             */
            Matrix(const size_t n_row, const size_t n_col) :
            data_(n_row*n_col),
            column_(n_col),
            n_row_(n_row) {
                for (size_t col = 0; col<n_col; ++col) {
                    column_[col] = data_.data() + col*n_row;
                }
                return;
            }

            /**
             * Destructor
             */
            virtual ~Matrix() {
                return;
            }

            /**
             * Get the number of rows in the matrix.
             * 
             * @return the number of rows in the matrix.
             */
            inline size_t GetNRow() const {
                return n_row_;
            }

            /**
             * Get the number of columns in the matrix.
             * 
             * @return the number of columns in the matrix.
             */
            inline size_t GetNCol() const {
                return column_.size();
            }

            /**
             * Get a reference to the idx-th entry of the matrix.
             * 
             * @param idx single index of the entry.
             * 
             * @return a reference to the idx-th entry.
             */
            inline Scalar& operator()(const size_t idx) {
                return data_[idx];
            }

            /**
             * Get a constant reference to the idx-th entry of the matrix.
             * 
             * @param idx single index of the entry.
             * 
             * @return a constant reference to the idx-th entry.
             */
            inline const Scalar& operator()(const size_t idx) const {
                return data_[idx];
            }

            /**
             * Get a reference to the entry of the matrix at given row and column.
             * 
             * @param row row of the entry.
             * @param col column of the entry.
             * 
             * @return a reference to the entry.
             */
            inline Scalar& operator()(const size_t row, const size_t col) {
                return column_[col][row];
            }

            /**
             * Get a constant reference to the entry of the matrix at given row and column.
             * 
             * @param row row of the entry.
             * @param col column of the entry.
             * 
             * @return a constant reference to the entry.
             */
            inline Scalar operator()(const size_t row, const size_t col) const {
                return column_[col][row];
            }
        private:
            /// Matrix entries ordered by columns.
            std::vector<Scalar> data_;
            /// Pointers to the beginning of each column.
            std::vector<Scalar*> column_;
            /// Number of rows.
            size_t n_row_;
    };

    template <typename ScalarMatrix, typename ScalarVector>
    auto SolveDiag(const Matrix<ScalarMatrix> &A, const std::vector<ScalarVector> &b) {
        using ScalarOutput = decltype(A(0,0) * b[0]);
        size_t n = b.size();
        std::vector<ScalarOutput> x(n);
        for (size_t i = 0; i<n; ++i) {
            x[i] = b[i] / A(i,i);
        }
        return x;
    }

    template <typename ScalarMatrix, typename ScalarVector>
    auto SolveTriU(const Matrix<ScalarMatrix> &A, const std::vector<ScalarVector> &b) {
        using ScalarOutput = decltype(A(0,0) * b[0]);
        size_t n = b.size();
        std::vector<ScalarOutput> x(n);
        x[n-1] = b[n-1] / A(n-1,n-1);
        size_t i = n-2;
        for (size_t i_ = 0; i_<n-1; ++i_) {
            ScalarOutput s = 0.0;
            for (size_t j = i+1; j<n; ++j) {
                s += A(i, j) * x[j];
            }
            x[i] = (b[i] - s) / A(i,i);
            --i;
        }
        return x;
    }

    template <typename ScalarMatrix, typename ScalarVector>
    auto SolveTriL(const Matrix<ScalarMatrix> &A, const std::vector<ScalarVector> &b) {
        using ScalarOutput = decltype(A(0,0) * b[0]);
        size_t n = b.size();
        std::vector<ScalarOutput> x(n);
        x[0] = b[0] / A(0,0);
        for (size_t i = 1; i<n; ++i) {
            ScalarOutput s = 0.0;
            for (size_t j = 0; j<i; ++j) {
                s += A(i, j) * x[j];
            }
            x[i] = (b[i] - s) / A(i,i);
        }
        return x;
    }

}  // namespace linalg

}  // namespace eptlib

#endif  // EPTLIB_LINALG_MATRIX_H_
