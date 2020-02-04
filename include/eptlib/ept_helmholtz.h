/*****************************************************************************
*
*     Program: EPTlib
*     Author: Alessandro Arduino <a.arduino@inrim.it>
*
*  MIT License
*
*  Copyright (c) 2020  Alessandro Arduino
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

#ifndef EPTLIB_EPT_HELMHOLTZ_H_
#define EPTLIB_EPT_HELMHOLTZ_H_

#include "eptlib/ept_interface.h"

#include "eptlib/config.h"
#include "eptlib/finite_difference.h"
#include "eptlib/shape.h"
#include "eptlib/util.h"

namespace eptlib {

/**
 * Implementation of the Helmholtz-based EPT method.
 */
class EPTHelmholtz : public EPTInterface {
    public:
        /**
         * Constructor.
         * 
         * @tparam T,U iterator typenames.
         * 
         * @param freq operative frequency of the MRI scanner.
         * @param nn number of voxels in each direction.
         * @param dd voxel sizes in each direction.
         * @param shape mask over which apply the finite difference scheme.
         * 
         * The number of Tx and Rx channels is fixed equal to one.
         */
        template <typename T, typename U>
        EPTHelmholtz(const real_t freq, const T &nn, const U &dd, const Shape &shape);
        /**
         * Virtual destructor.
         */
        virtual ~EPTHelmholtz();
        /**
         * Perform the Helmholtz-based EPT.
         * 
         * @return an error index about the state of the tomography.
         * 
         * Three variants of the EPT technique are implemented: when both the
         * tx sensitivity and the trx phase are known the complete method is
         * used, otherwise the magnitude-based or the phase-based
         * approximations are applied.
         */
        virtual EPTlibError_t Run() override;
    protected:
        /// Filter for the Laplacian computation.
        FDLaplacianKernel fd_lapl_;
};

// ---------------------------------------------------------------------------
// -------------------------  Implementation detail  -------------------------
// ---------------------------------------------------------------------------

// EPTHelmholtz constructor
template <typename T, typename U>
EPTHelmholtz::
EPTHelmholtz(const real_t freq, const T &nn, const U &dd, const Shape &shape) :
    EPTInterface(freq,nn,dd), fd_lapl_(shape) {
    return;
}

}  // namespace eptlib

#endif  // EPTLIB_EPT_HELMHOLTZ_H_
