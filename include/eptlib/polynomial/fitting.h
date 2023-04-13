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

#ifndef EPTLIB_POLYNOMIAL_FITTING_H_
#define EPTLIB_POLYNOMIAL_FITTING_H_

#include <tuple>
#include <vector>

#include "eptlib/shape.h"

#include "eptlib/linalg/matrix.h"
#include "eptlib/linalg/regression.h"
#include "eptlib/linalg/vector.h"

namespace eptlib {

namespace polynomial {

    /**
     * @brief Get the number of three-dimensional monomials of given maximum degree.
     * 
     * @param degree Maximum degree of the monomial.
     * 
     * @return the number of monomials. 
     */
    size_t GetNumberOfMonomials(const size_t degree);

    /**
     * @brief Evaluate all the monomials of given degree.
     * 
     * @param x x-coordinate of where the monomials are evaluted.
     * @param y y-coordinate of where the monomials are evaluted.
     * @param z z-coordinate of where the monomials are evaluted.
     * @param n degree of the monomials.
     * 
     * @return vector with the evaluated monomials.
     * 
     * The monomials are all those written as x^i * y^j * z^k, with
     * i+j+k = n. The elements are ordered as follow,
     * 1. (i=n  , j=0, k=0);
     * 2. (i=n-1, j=1, k=0);
     * 3. (i=n-2, j=2, k=0);
     * ...
     */
    std::vector<double> EvaluateMonomials(const double x, const double y, const double z, const size_t n);

    /**
     * @brief Fill the design matrix for polynomial fitting with a basis of monomials.
     * 
     * @param x x-coordinates where the polynomial is evaluated.
     * @param y y-coordinates where the polynomial is evaluated.
     * @param z z-coordinates where the polynomial is evaluated.
     * @param degree degree of the fitting polynomial.
     * 
     * @return the design matrix.
     * 
     * The columns of the matrix are ordered as follow,
     * (1, x, y, z, x^2, x*y, y^2, x*z, y*z, z^2, x^3, x^2*y, ...).
     */
    eptlib::linalg::Matrix<double> DesignMatrixWithMonomialsBasis(const std::vector<double> &x,
        const std::vector<double> &y, const std::vector<double> &z, const size_t degree);

    /**
     * @brief Compute the coefficients of the polynomial function fitting the data.
     * 
     * @param x x-coordinates where the data are provided.
     * @param y y-coordinates where the data are provided.
     * @param z z-coordinates where the data are provided.
     * @param degree degree of the fitting polynomial.
     * @param values values of the provided data.
     * 
     * @return a std::tuple with:
     *     1) a vector with the polynomial coefficients;
     *     2) the normalized chi-squared statistic of the fitting.
     */
    template <typename Scalar>
    std::tuple<std::vector<Scalar>, double> Fitting(const std::vector<double> &x,
        const std::vector<double> &y, const std::vector<double> &z, const size_t degree,
        const std::vector<Scalar> &values // , const std::vector<double> &weights
        ) {
        eptlib::linalg::Matrix<double> F = DesignMatrixWithMonomialsBasis(x, y, z, degree);
        return eptlib::linalg::LinearRegression(F, values);
    }

}  // namespace polynomial

}  // namespace eptlib

#endif  // EPTLIB_POLYNOMIAL_FITTING_H_
