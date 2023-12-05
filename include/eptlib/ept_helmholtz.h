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

#ifndef EPTLIB_EPT_HELMHOLTZ_H_
#define EPTLIB_EPT_HELMHOLTZ_H_

#include "eptlib/ept_interface.h"

#include <variant>

#include "eptlib/shape.h"
#include "eptlib/util.h"

#include "eptlib/filter/anatomical_savitzky_golay.h"
#include "eptlib/filter/postprocessing.h"
#include "eptlib/filter/savitzky_golay.h"

namespace eptlib {

/**
 * Implementation of the Helmholtz-based EPT method.
 */
class EPTHelmholtz : public EPTInterface {
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
         * @param trx_phase_is_wrapped if true, the TRx phase is wrapped.
         * @param compute_variance if true, the variance map will be computed.
         * @param weight_param parameter of the weight function, used only with
         *     the anatomical Savitzky-Golay filter (default: 0.05).
         * 
         * The number of Tx and Rx channels is fixed equal to one.
         */
        EPTHelmholtz(const size_t n0, const size_t n1, const size_t n2,
            const double d0, const double d1, const double d2,
            const double freq, const Shape &window, const int degree = 2,
            const bool trx_phase_is_wrapped = false,
            const bool compute_variance = false,
            const double weight_param = 0.05);

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
         * Tx sensitivity and the TRx phase are known the complete method is
         * used, otherwise the magnitude-based or the phase-based
         * approximations are applied.
         */
        virtual EPTlibError Run() override;

        /**
         * @brief Set or unset the compute variance flag.
         * 
         * @return the updated compute variance flag.
         */
        inline bool ToggleComputeVariance() {
            compute_variance_ = !compute_variance_;
            return compute_variance_;
        }

        /**
         * @brief Get the compute variance flag.
         * 
         * @return the compute variance flag.
         */
        inline bool ComputeVariance() const {
            return compute_variance_;
        }

        /**
         * @brief Get the computed variance map of the electric conductivity.
         * 
         * @return reference of the pointer to the computed variance map of the electric conductivity.
         */
        inline std::unique_ptr<Image<double> >& GetElectricConductivityVariance() {
            return variance_sigma_;
        }

        /**
         * @brief Get the computed variance map of the relative permittivity.
         * 
         * @return reference of the pointer to the computed variance map of the relative permittivity.
         */
        inline std::unique_ptr<Image<double> >& GetRelativePermittivityVariance() {
            return variance_epsr_;
        }

        /**
         * @brief Check if the variance of the electric conductivity is set.
         * 
         * @return true if the variance of the electric conductivity is set.
         * @return false if the variance of the electric conductivity is not set.
         */
        inline bool ThereIsElectricConductivityVariance() const {
            return variance_sigma_!=nullptr;
        }

        /**
         * @brief Check if the variance of the relative permittivity is set.
         * 
         * @return true if the variance of the relative permittivity is set.
         * @return false if the variance of the relative permittivity is not set.
         */
        inline bool ThereIsRelativePermittivityVariance() const {
            return variance_epsr_!=nullptr;
        }
    private:
        /// Savitzky-Golay filter for the derivative computation.
        std::variant<std::monostate, filter::SavitzkyGolay, filter::AnatomicalSavitzkyGolay> sg_filter_;
        /// Mask over which apply the Savitzky-Golay filter.
        Shape sg_window_;
        /// Degree of the interpolating polynomial for the Savitzky-Golay filter.
        int sg_degree_;

        /// Variance map of the electric conductivity.
        std::unique_ptr<Image<double> > variance_sigma_;
        /// Variance map of the relative permittivity.
        std::unique_ptr<Image<double> > variance_epsr_;
        /// If true, compute the result variance.
        bool compute_variance_;

        /// Parameter of the weight function, used only with the anatomical Savitzky-Golay filter.
        double weight_param_;

        /// Perform the complete Helmholtz-based EPT.
        void CompleteEPTHelm();
        /// Perform the phase approximated Helmholtz-based EPT.
        void PhaseEPTHelm();
        /// Perform the magnitude approximated Helmholtz-based EPT.
        void MagnitudeEPTHelm();
};

}  // namespace eptlib

#endif  // EPTLIB_EPT_HELMHOLTZ_H_
