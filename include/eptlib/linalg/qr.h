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

#ifndef EPTLIB_LINALG_QR_H_
#define EPTLIB_LINALG_QR_H_

#include "eptlib/linalg/matrix.h"
#include "eptlib/linalg/vector.h"

#include <cstring>
#include <numeric>
#include <tuple>

namespace eptlib {

namespace linalg {

    /**
     * @brief Compute the elementary reflector.
     * 
     * The elementary reflector is stored in place of the range of elements, whose last
     * value is used to store the factor `pi'.
     * 
     * @tparam ForwardIt forward iterator typename.
     * 
     * @param first,last range of elements.
     */
    template <typename ForwardIt>
    void HouseholderReflector(ForwardIt first, ForwardIt last) {
        double sigma = eptlib::linalg::Norm2(first, last-1);
        sigma = std::copysign(sigma, *first);
        *first += sigma;
        *(last-1) = sigma * (*first);
        return;
    }

    /**
     * @brief Apply the elementary reflector to a vector.
     * 
     * @tparam InputIt1 input iterator typename for the vector.
     * @tparam InputIt2 input iterator typename for the elementary reflector.
     * 
     * @param x_first,xlast range of elements for the vector.
     * @param u_first beginning of the range of elements for the elementary reflector.
     * @param pi factor `pi' of the elementary reflector.
     */
    template <typename InputIt1, typename InputIt2>
    void HouseholderLeft(InputIt1 x_first, InputIt1 x_last, InputIt2 u_first, const double pi) {
        using tau_t = decltype(*x_first * (*u_first));
        tau_t tau = std::inner_product(x_first, x_last, u_first, static_cast<tau_t>(0.0)) / pi;
        while (x_first != x_last) {
            *x_first -= tau * (*u_first);
            ++x_first;
            ++u_first;
        }
        return;
    }

    /**
     * @brief Perform the QR decomposition of a vertical full-rank matrix.
     * 
     * The upper triangular part of the output stores the matrix R. The lower triangular
     * part of the output stores the elementary reflectors describing the transpose of Q.
     * 
     * @param A vertical full-rank matrix to be decomposed (column-major).
     * 
     * @return a matrix containing the QR decomposition.
     */
    Matrix<double> QRDecomposition(const Matrix<double> &A);

    /**
     * @brief Perform the QR decomposition with column pivoting.
     * 
     * The upper triangular part of the output stores the matrix R. The lower triangular
     * part of the output stores the elementary reflectors describing the transpose of Q.
     * 
     * @param A matrix to be decomposed (column-major).
     * 
     * @return a std::tuple with:
     *     1) a matrix containing the QR decomposition;
     *     2) a vector of indices defining the permutation;
     */
    std::tuple<Matrix<double>, std::vector<size_t>> QRDecompositionWithColumnPivoting(const Matrix<double> &A);

    /**
     * @brief Compute the numerical rank of a matrix given its rank-revealing QR decomposition.
     * 
     * For most cases, QR decomposition with column pivoting works as a rank-revealing QR decomposition.
     * 
     * @param RRQR a rank revealing QR decomposition of the matrix.
     * 
     * @return the numerical rank of the matrix.
     */
    size_t QRGetRank(const Matrix<double> &RRQR);

    /**
     * @brief Solve a linear system in the least square sense using the QR decomposition.
     * 
     * @tparam Scalar numerical type of the forcing term entries.
     * 
     * @param QR QR decomposition of the coefficient matrix of the linear system.
     * @param b forcing term of the linear system.
     * 
     * @return a std::tuple with:
     *     1) a vector solving the linear system in the least square sense;
     *     2) the Euclidean norm of the residual.
     */
    template <typename Scalar>
    std::tuple<std::vector<Scalar>, double> QRSolve(const Matrix<double> &QR, const std::vector<Scalar> &b) {
        const size_t n_row = b.size();
        const size_t n_col = QR.GetNCol();
        // compute c = Q^T * b
        std::vector<Scalar> c(n_row);
        std::memcpy(c.data(), b.data(), n_row*sizeof(Scalar));
        for (size_t col = 0; col < n_col; ++col) {
            HouseholderLeft(c.begin()+col, c.end(), QR.begin(col)+col+1, QR(n_row+1,col));
        }
        // solve R * x = c
        std::vector<Scalar> x = SolveTriU(QR, c);
        // compute the residual
        double chi = Norm2(c.begin()+n_col, c.end());
        return {x, chi};
    }

}  // namespace linalg

}  // namespace eptlib

#endif  // EPTLIB_LINALG_QR_H_
