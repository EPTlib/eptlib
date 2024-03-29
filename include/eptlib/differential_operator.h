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

#ifndef EPTLIB_DIFFERENTIAL_OPERATOR_H_
#define EPTLIB_DIFFERENTIAL_OPERATOR_H_

namespace eptlib {

    enum class DifferentialOperator {
        /// Zero order derivative (field approximation)
        Field,
        /// First order derivative along X
        GradientX,
        /// First order derivative along Y
        GradientY,
        /// First order derivative along Z
        GradientZ,
        /// Second order derivative along X
        GradientXX,
        /// Second order derivative along Y
        GradientYY,
        /// Second order derivative along Z
        GradientZZ,
        /// Laplacian
        Laplacian,
        /// Fictitious label to denote the end of differential operators
        END,
    };

}  // namespace eptlib

#endif  // EPTLIB_DIFFERENTIAL_OPERATOR_H_
