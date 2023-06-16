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
             * Construct a matrix with zero rows and columns.
             */
            Matrix() : Matrix(0,0) {
                return;
            }
            
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
             * Copy constructor.
             * 
             * @param o matrix to be copied.
             */
            Matrix(const Matrix<Scalar> &o) :
            Matrix(o.GetNRow(), o.GetNCol()) {
                std::copy(o.data_.begin(), o.data_.end(), data_.begin());
                return;
            }

            /**
             * Move constructor.
             * 
             * @param o matrix to be moved.
             */
            Matrix(Matrix<Scalar> &&o) noexcept :
            Matrix() {
                swap(*this, o);
                return;
            }

            /**
             * Copy assignment.
             * 
             * @param o matrix to be copied.
             */
            Matrix<Scalar>& operator=(Matrix<Scalar> o) {
                swap(*this, o);
                return *this;
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
             * Get a copy of the idx-th entry of the matrix.
             * 
             * @param idx single index of the entry.
             * 
             * @return a copy of the idx-th entry.
             */
            inline Scalar operator()(const size_t idx) const {
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
             * Get a copy of the entry of the matrix at given row and column.
             * 
             * @param row row of the entry.
             * @param col column of the entry.
             * 
             * @return a copy of the entry.
             */
            inline Scalar operator()(const size_t row, const size_t col) const {
                return column_[col][row];
            }

            /**
             * @brief Get an iterator to the first element of the matrix.
             * 
             * @return an iterator to the first element of the matrix.
             */
            inline typename std::vector<Scalar>::iterator begin() {
                return data_.begin();
            }

            /**
             * @brief Get an iterator to the first element of the matrix.
             * 
             * @return an iterator to the first element of the matrix.
             */
            inline typename std::vector<Scalar>::const_iterator begin() const {
                return data_.begin();
            }

            /**
             * @brief Get an iterator following the last element of the matrix.
             * 
             * @return an iterator following the last element of the matrix.
             */
            inline typename std::vector<Scalar>::iterator end() {
                return data_.end();
            }

            /**
             * @brief Get an iterator following the last element of the matrix.
             * 
             * @return an iterator following the last element of the matrix.
             */
            inline typename std::vector<Scalar>::const_iterator end() const {
                return data_.end();
            }

            /**
             * @brief Get an iterator to the first element of a column of the matrix.
             * 
             * @param col column of interest.
             * 
             * @return an iterator to the first element of a column of the matrix.
             */
            inline typename std::vector<Scalar>::iterator begin(const size_t col) {
                return this->begin() + col*n_row_;
            }

            /**
             * @brief Get an iterator to the first element of a column of the matrix.
             * 
             * @param col column of interest.
             * 
             * @return an iterator to the first element of a column of the matrix.
             */
            inline typename std::vector<Scalar>::const_iterator begin(const size_t col) const {
                return this->begin() + col*n_row_;
            }

            /**
             * @brief Get an iterator following the last element of a column of the matrix.
             * 
             * @param col column of interest.
             * 
             * @return an iterator following the last element of a column of the matrix.
             */
            inline typename std::vector<Scalar>::iterator end(const size_t col) {
                return this->begin(col+1);
            }

            /**
             * @brief Get an iterator following the last element of a column of the matrix.
             * 
             * @param col column of interest.
             * 
             * @return an iterator following the last element of a column of the matrix.
             */
            inline typename std::vector<Scalar>::const_iterator end(const size_t col) const {
                return this->begin(col+1);
            }

            
            /**
             * Swap the content of two matrices.
             * 
             * @param first,second matrices to be swapped.
             */
            friend void swap(Matrix<Scalar> &first, Matrix<Scalar> &second) {
                using std::swap;
                swap(first.data_, second.data_);
                swap(first.column_, second.column_);
                swap(first.n_row_, second.n_row_);
            }
        private:
            /// Matrix entries ordered by columns.
            std::vector<Scalar> data_;
            /// Pointers to the beginning of each column.
            std::vector<Scalar*> column_;
            /// Number of rows.
            size_t n_row_;
    };

    /**
     * @brief Solve a square diagonal system.
     * 
     * @tparam ScalarMatrix numeric typename of the matrix entries.
     * @tparam ScalarVector numeric typename of the vector entries.
     * 
     * @param A matrix of the coefficients (non-diagonal elements assumed null).
     * @param b vector of the forcing terms.
     * 
     * @return vector of the solution. 
     */
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

    /**
     * @brief Solve an upper triangular square system.
     * 
     * @tparam ScalarMatrix numeric typename of the matrix entries.
     * @tparam ScalarVector numeric typename of the vector entries.
     * 
     * @param A matrix of the coefficients (non-upper triangular elements assumed null).
     * @param b vector of the forcing terms.
     * 
     * @return vector of the solution. 
     */
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

    /**
     * @brief Solve a lower triangular square system.
     * 
     * @tparam ScalarMatrix numeric typename of the matrix entries.
     * @tparam ScalarVector numeric typename of the vector entries.
     * 
     * @param A matrix of the coefficients (non-lower triangular elements assumed null).
     * @param b vector of the forcing terms.
     * 
     * @return vector of the solution. 
     */
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
