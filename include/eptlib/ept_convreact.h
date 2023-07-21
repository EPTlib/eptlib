/*****************************************************************************
*
*     Program: EPTlib
*     Author: Alessandro Arduino <a.arduino@inrim.it>
*
*  MIT License
*
*  Copyright (c) 2020-2023  Alessandro Arduino
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

#include <optional>
#include <variant>

#include "eptlib/shape.h"
#include "eptlib/util.h"

#include "eptlib/filter/anatomical_savitzky_golay.h"
#include "eptlib/filter/savitzky_golay.h"

namespace eptlib {

/**
 * Implementation of the convection-reaction EPT method.
 */
class EPTConvReact : public EPTInterface {
    public:
        /**
         * Constructor.
         * 
         * @param n0 number of voxels along direction x.
         * @param n1 number of voxels along direction y.
         * @param n2 number of voxels along direction z.
         * @param d0 resolution in meter along direction x.
         * @param d1 resolution in meter along direction y.
         * @param d2 resolution in meter along direction z.
         * @param freq operative frequency of the MRI scanner.
         * @param window mask over which apply the finite difference scheme.
         * @param degree degree of the interpolating polynomial for the finite
         *     difference scheme (default: 2).
         * 
         * The number of Tx and Rx channels is fixed equal to one.
         */
        EPTConvReact(const size_t n0, const size_t n1, const size_t n2,
            const double d0, const double d1, const double d2,
            const double freq, const Shape &window, const int degree = 2);

        /**
         * Virtual destructor.
         */
        virtual ~EPTConvReact();

        /**
         * Perform the convection-reaction EPT.
         * 
         * @return an error index about the state of the tomography.
         */
        virtual EPTlibError Run() override;

        /**
         * @brief Get the volume tomography flag.
         * 
         * @return the volume tomography flag.
         */
        inline bool VolumeTomography() const {
            return !slice_index_.has_value();
        }

        /**
         * Set the selected plane index for plane tomography.
         * 
         * @return a Success or OutOfRange error.
         */
        inline EPTlibError SelectSlice(const size_t slice_index) {
            if (slice_index >= nn_[2]) {
                return EPTlibError::OutOfRange;
            }
            slice_index_ = slice_index;
            return EPTlibError::Success;
        }

        /**
         * Set the Dirichlet boundary conditions.
         */
        inline void SetDirichlet(const double dirichlet_epsr, const double dirichlet_sigma) {
            dirichlet_epsr_ = dirichlet_epsr;
            dirichlet_sigma_ = dirichlet_sigma;
        }

        /**
         * Set the artificial diffusion stabilisation.
         * 
         * @param artificial_diffusion the coefficient of the artificial diffusion stabilisation.
         */
        inline void SetArtificialDiffusion(const double artificial_diffusion) {
            artificial_diffusion_ = artificial_diffusion;
        }

        /**
         * @brief Check if the artificial diffusion stabilisation is set.
         * 
         * @return true if the artificial diffusion stabilisation is set.
         * @return false if the artificial diffusion stabilisation is not set.
         */
        inline bool ThereIsArtificialDiffusion() const {
            return artificial_diffusion_.has_value();
        }

        /**
         * Get the number of iterations to solve the linear system.
         * 
         * @return the number of iterations to solve the linear system.
         */
        inline std::ptrdiff_t GetSolverIterations() const {
            return solver_iterations_;
        }

        /**
         * Get the estimated error in the linear system solution.
         * 
         * @return the estimated error in the linear system solution.
         */
        inline double GetSolverResidual() const {
            return solver_residual_;
        }
    private:
        /// Selected slice index for plane tomography.
        std::optional<size_t> slice_index_;
        /// Artificial diffusion coefficient.
        std::optional<double> artificial_diffusion_;
        /// Dirichlet condition of relative permittivity.
        double dirichlet_epsr_;
        /// Dirichlet condition of electric conductivity.
        double dirichlet_sigma_;
        
        /// Savitzky-Golay filter for the derivative computation.
        std::variant<std::monostate, filter::SavitzkyGolay, filter::AnatomicalSavitzkyGolay> sg_filter_;
        /// Mask over which apply the Savitzky-Golay filter.
        Shape sg_window_;
        /// Degree of the interpolating polynomial for the Savitzky-Golay filter.
        int sg_degree_;

        /// The number of iterations to solve the linear system.
        std::ptrdiff_t solver_iterations_;
        /// The estimated error in the linear system solution.
        double solver_residual_;

        /// Perform the complete convection-reaction EPT.
        EPTlibError CompleteEPTConvReact();
        /// Perform the phase approximated convection-reaction EPT.
        EPTlibError PhaseEPTConvReact();
};

}  // namespace eptlib

#endif  // EPTLIB_EPT_CONVREACT_H_
