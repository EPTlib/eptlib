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

#ifndef EPTLIB_EPT_CONVREACT_H_
#define EPTLIB_EPT_CONVREACT_H_

#include "eptlib/ept_interface.h"

#include <array>

#include "eptlib/finite_difference.h"
#include "eptlib/shape.h"
#include "eptlib/util.h"

namespace eptlib {

/**
 * Implementation of the convection-reaction EPT method.
 */
class EPTConvReact : public EPTInterface {
    public:
        /**
         * Constructor.
         * 
         * @param freq operative frequency of the MRI scanner.
         * @param nn number of voxels in each direction.
         * @param dd voxel sizes in each direction.
         * @param shape mask over which apply the finite difference scheme.
         * 
         * The number of Tx and Rx channels is fixed equal to one.
         */
        EPTConvReact(const double freq, const std::array<int,NDIM> &nn,
            const std::array<double,NDIM> &dd, const Shape &shape);
        /**
         * Virtual destructor.
         */
        virtual ~EPTConvReact();
        /**
         * Perform the convection-reaction EPT.
         * 
         * @return an error index about the state of the tomography.
         */
        virtual EPTlibError_t Run() override;
        /**
         * Set the dirichlet boundary conditions.
         * 
         * @return a Success or Unknown error.
         */
        EPTlibError_t SetDirichlet(const double dir_epsr, const double dir_sigma);
        /**
         * Set the selected plane index for plane tomography.
         * 
         * @return a Success or WrongDataFormat error.
         */
        EPTlibError_t SelectPlane(const int plane_idx);
        /**
         * Set the tomography on a volume domain.
         * 
         * @return a Success or Unknown error.
         */
        EPTlibError_t SetVolumeTomography();
        /**
         * Unset the tomography on a volume domain.
         * 
         * @return a Success or Unknown error.
         */
        EPTlibError_t UnsetVolumeTomography();
        /**
         * Set the artificial diffusion stabilisation.
         * 
         * @param diff_coeff the coefficient of the artificial diffusion stabilisation.
         * 
         * @return a Success or Unknown error.
         */
        EPTlibError_t SetArtificialDiffusion(const double diff_coeff);
        /**
         * Unset the artificial diffusion stabilisation.
         * 
         * @return a Success or Unknown error.
         */
        EPTlibError_t UnsetArtificialDiffusion();
    private:
        /// Dirichlet condition of relative permittivity.
        double dir_epsr_;
        /// Dirichlet condition of electric conductivity.
        double dir_sigma_;
        /// Selected plane index.
        int plane_idx_;
        /// Volume tomography flag.
        bool is_volume_;
        /// Artificial diffusion coefficient.
        double diff_coeff_;
        /// Artificial diffusion flag.
        bool thereis_diff_;
        /// Filter for the derivative computation.
        FDSavitzkyGolayFilter fd_filter_;
        /// Perform the complete convection-reaction EPT.
        EPTlibError_t CompleteEPTConvReact();
        /// Perform the phase approximated convection-reaction EPT.
        EPTlibError_t PhaseEPTConvReact();
};

}  // namespace eptlib

#endif  // EPTLIB_EPT_CONVREACT_H_
