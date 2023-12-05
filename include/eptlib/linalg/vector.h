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

#ifndef EPTLIB_LINALG_VECTOR_H_
#define EPTLIB_LINALG_VECTOR_H_

#include <iterator>
#include <limits>
#include <type_traits>

namespace eptlib {

namespace linalg {

    /**
     * @brief Look for the maximum absolute value in a vector.
     * 
     * @tparam ForwardIt forward iterator typename.
     * 
     * @param first,last range of elements.
     * 
     * @return the maximum absolute value in the vector.
     */
    template <typename ForwardIt>
    double MaxAbs(ForwardIt first, ForwardIt last) {
        double eta = 0;
        while (first != last) {
            double tmp = std::abs(*first++);
            if (std::isnan(tmp)) {
                return std::numeric_limits<double>::quiet_NaN();
            }
            eta = std::max(eta, tmp);
        }
        return eta;
    }

    /**
     * @brief Compute the Euclidean norm of a vector.
     * 
     * @tparam ForwardIt forward iterator typename.
     * 
     * @param first,last range of elements.
     * 
     * @return the Euclidean norm of the vector.
     */
    template <typename ForwardIt>
    double Norm2(ForwardIt first, ForwardIt last) {
        double eta = MaxAbs(first, last);
        if (eta == 0.0) {
            return eta;
        }
        double sigma = 0;
        while (first != last) {
            double tmp = std::abs(*first++);
            if (tmp >= std::sqrt(std::numeric_limits<double>::epsilon())*eta) {
                sigma += tmp*tmp/eta/eta;
            }
        }
        sigma = std::sqrt(sigma)*eta;
        return sigma;
    }

    /**
     * @brief Change the order of the elements according to a given permutation.
     * 
     * @tparam ForwardIt forward iterator typename.
     * 
     * @param first,last range of elements.
     * @param p permutation vector.
     */
    template <typename ForwardIt>
    void Permute(ForwardIt first, ForwardIt last, std::vector<size_t> p) {
        using std::swap;
        size_t i = 0;
        while (i < p.size()) {
            if (p[i] != i) {
                std::iter_swap(first, first+(p[i]-i));
                swap(p[i], p[p[i]]);
            } else {
                ++i;
                ++first;
            }
        }
        return;
    }

}  // namespace linalg

}  // namespace eptlib

#endif  // EPTLIB_LINALG_VECTOR_H_
