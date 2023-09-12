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

#ifndef EPTLIB_FILTER_WEIGHTS_H_
#define EPTLIB_FILTER_WEIGHTS_H_

#include <algorithm>
#include <functional>
#include <vector>

namespace eptlib {

namespace filter {

    /**
     * @brief Evaluate the hard threshold function.
     * 
     * @param x Abscissa at which the hard threshold function is evaluated.
     * @param threshold Threshold parameter.
     * 
     * @return Evaluated hard threshold function.
     */
    double HardThreshold(const double x, const double threshold);

    /**
     * @brief Evaluate the Gaussian function
     * 
     * @param x Abscissa at which the Gaussian function is evaluated.
     * @param sigma Standard deviation of the Gaussian.
     * 
     * @return Evaluated Gaussian function.
     */
    double Gaussian(const double x, const double sigma);

}  // namespace filter

}  // namespace eptlib

#endif  // EPTLIB_FILTER_WEIGHTS_H_
