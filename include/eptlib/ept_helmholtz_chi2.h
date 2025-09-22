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

#ifndef EPTLIB_EPT_HELMHOLTZ_CHI2_H_
#define EPTLIB_EPT_HELMHOLTZ_CHI2_H_

#include "eptlib/ept_interface.h"

#include <vector>

#include "eptlib/shape.h"
#include "eptlib/util.h"

namespace eptlib {

/**
 * @brief Implementation of the Helmholtz-based EPT method with pixel-wise
 * optimised kernel shape.
 * 
 * The method is phase-based and provides only an estimation of the electric
 * conductivity. Because of the involved non-linearities, it can operate only
 * on unwrapped phase maps.
 */
class EPTHelmholtzChi2 : public EPTInterface {
    public:
        /**
         * @brief Constructor
         * 
         * @param n0 number of voxels along direction x.
         * @param n1 number of voxels along direction y.
         * @param n2 number of voxels along direction z.
         * @param d0 resolution in meter along direction x.
         * @param d1 resolution in meter along direction y.
         * @param d2 resolution in meter along direction z.
         * @param freq operative frequency of the MRI scanner.
         * @param windows list of masks over which apply the finite difference scheme.
         * @param degree degree of the interpolating polynomial for the finite
         *     difference scheme (default: 2).
         * @param admit_unphysical_values if true, unphysical results are admitted in the maps.
         * 
         * The number of Tx and Rx channels is fixed equal to one.
        */
        EPTHelmholtzChi2(const size_t n0, const size_t n1, const size_t n2,
            const double d0, const double d1, const double d2,
            const double freq, const std::vector<Shape> &windows,
            const int degree = 2, const bool admit_unphysical_values = false,
            const double weight_param = 0.05);

        /**
         * @brief Virtual destructor.
         */
        virtual ~EPTHelmholtzChi2();

        /**
         * @brief Perform the Helmholtz-based EPT with pixel-wise optimised kernel shape.
         * 
         * @return an error index about the state of the tomography.
         * 
         * Only the phase-based approximation of Helmholtz-based EPT is implemented.
         * The TRx phase must be provided.
         */
        virtual EPTlibError Run() override;

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

        /**
         * @brief Get the index map of the selected kernel shapes for electric conductivity.
         * 
         * @return reference of the pointer to the index map for electric conductivity.
         */
        inline std::unique_ptr<Image<int> >& GetElectricConductivityIndex() {
            return index_sigma_;
        }

        /**
         * @brief Check if the index map for electric conductivity is set.
         * 
         * @return true if the index map for electric conductivity is set.
         * @return false if the index map for electric conductivity is not set.
         */
        inline bool ThereIsElectricConductivityIndex() const {
            return index_sigma_!=nullptr;
        }

        /**
         * @brief Get the index map of the selected kernel shapes for relative permittivity.
         * 
         * @return reference of the pointer to the index map for the relative permittivity.
         */
        inline std::unique_ptr<Image<int> >& GetRelativePermittivityIndex() {
            return index_epsr_;
        }

        /**
         * @brief Check if the index map for relative permittivity is set.
         * 
         * @return true if the index map for relative permittivity is set.
         * @return false if the index map for relative permittivity is not set.
         */
        inline bool ThereIsRelativePermittivityIndex() const {
            return index_epsr_!=nullptr;
        }

        /**
         * @brief Set/unset unphysical values as admittable.
         * 
         * @return if unphysical values are admittable.
         */
        inline bool ToggleAdmitUnphysicalValues() {
            admit_unphysical_values_ = !admit_unphysical_values_;
            return admit_unphysical_values_;
        }

        /**
         * @brief Get the admit unphysical values flag.
         * 
         * @return the admit unphysical values flag.
         */
        inline bool AdmitUnphysicalValues() const {
            return admit_unphysical_values_;
        }
    private:
        /// List of masks over which apply the finite difference scheme.
        std::vector<Shape> windows_;
        /// Degree of the interpolating polynomial for the finite difference scheme.
        int degree_;

        /// Variance map of the electric conductivity.
        std::unique_ptr<Image<double> > variance_sigma_;
        /// Variance map of the relative permittivity.
        std::unique_ptr<Image<double> > variance_epsr_;

        /// Pixel-wise selected kernel shape for electric conductivity.
        std::unique_ptr<Image<int> > index_sigma_;
        /// Pixel-wise selected kernel shape for relative permittivity.
        std::unique_ptr<Image<int> > index_epsr_;

        /// Unphysical values flag.
        bool admit_unphysical_values_;

        /// Parameter of the weight function, used only with the anatomical Savitzky-Golay filter.
        double weight_param_;
};

}  // namespace eptlib

#endif  // EPTLIB_EPT_HELMHOLTZ_CHI2_H_
