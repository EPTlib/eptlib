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

#ifndef EPTLIB_EPT_GRADIENT_H_
#define EPTLIB_EPT_GRADIENT_H_

#include "eptlib/ept_interface.h"

#include <array>
#include <complex>
#include <optional>
#include <vector>

#include <boost/dynamic_bitset.hpp>

#include "eptlib/finite_difference.h"
#include "eptlib/shape.h"
#include "eptlib/util.h"

namespace eptlib {

/**
 * Structure of a seed point for the gradient inversion.
 */
struct SeedPoint {
    /// Integer coordinates of the seed point.
    std::array<size_t, N_DIM> ijk;
    /// Relative permittivity in the point.
    double epsr;
    /// Electric conductivity in the point.
    double sigma;
};

/**
 * Codes for the run modality of the gradient EPT method.
 */
enum class EPTGradientRun {
    /// Complete application of gradient EPT.
    FULL = 0,
    /// Solve the local system and return the first approximation.
    LOCAL,
    /// Solve the gradient inversion problem if a local solution is provided.
    GRADIENT,
};

/**
 * Implementation of the gradient EPT method.
 */
class EPTGradient : public EPTInterface {
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
         * @param shape mask over which apply the finite difference scheme.
         * @param n_tx_ch number of transmit channels.
         * @param degree degree of the interpolating polynomial for the finite
         *     difference scheme (default: 2).
         * 
         * The number of Rx channels is fixed equal to one.
         */
        EPTGradient(const size_t n0, const size_t n1, const size_t n2,
            const double d0, const double d1, const double d2,
            const double freq, const Shape &shape, const size_t n_tx_ch,
            const int degree = 2, const EPTGradientRun run_mode = EPTGradientRun::FULL);

        /**
         * Virtual destructor.
         */
        virtual ~EPTGradient();

        /**
         * Perform the gradient EPT.
         * 
         * @return an error index about the state of the tomography.
         */
        virtual EPTlibError Run() override;

        /**
         * Get the run mode of the next gradient EPT run.
         */
        inline EPTGradientRun GetRunMode() const {
            return run_mode_;
        }

        /**
         * Get the use seed points flag.
         * 
         * @return the seed points flag.
         */
        inline bool SeedPointsAreUsed() const {
            return seed_points_.has_value();
        }

        /**
         * Add a seed point to the list and set their use if not done before.
         * 
         * @param seed_point seed point to the added to the list.
         */
        inline void AddSeedPoint(const SeedPoint seed_point) {
            if (!SeedPointsAreUsed()) {
                seed_points_.emplace(0);
            }
            seed_points_->push_back(seed_point);
        }

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
         * Set the regularization coefficient.
         * 
         * @param lambda regularization coefficient.
         * 
         * @return a Success or WrongDataFormat error.
         */
        inline EPTlibError SetRegularizationCoefficient(const double lambda) {
            if (lambda < 0.0) {
                return EPTlibError::WrongDataFormat;
            }
            lambda_ = lambda;
            return EPTlibError::Success;
        }

        /**
         * Set the gradient tolerance for estimating the mask of homogeneous regions.
         * 
         * @return a Success or WrongDataFormat error.
         */
        inline EPTlibError SetGradientTolerance(const double gradient_tolerance) {
            if (gradient_tolerance < 0.0 || gradient_tolerance > 1.0) {
                return EPTlibError::WrongDataFormat;
            }
            gradient_tolerance_ = gradient_tolerance;
            return EPTlibError::Success;
        }

        /**
         * Get the homogeneous regions mask.
         * 
         * @return a constant reference to the mask.
         */
        inline const boost::dynamic_bitset<>& GetMask() const {
            return mask_;
        }

        /**
         * Get the cost functional.
         * 
         * @return the cost functional.
         */
        inline double GetCostFunctional() const {
            return cost_functional_;
        }

        /**
         * Get the cost regularization.
         * 
         * @return the cost regularization.
         */
        inline double GetCostRegularization() const {
            return cost_regularization_;
        }

        /**
         * Get the complex permittivity.
         * 
         * @return a constant reference to the complex permittivity.
         */
        inline const std::vector<std::complex<double> >& GetEpsC() const {
            return epsc_;
        }

        /**
         * Get the positive derivative of the complex permittivity logarithm.
         * 
         * @return a constant reference to the positive derivative.
         */
        inline const std::vector<std::complex<double> >& GetGPlus() const {
            return g_plus_;
        }

        /**
         * Get the longitudinal derivative of the complex permittivity logarithm.
         * 
         * @return a constant reference to the longitudinal derivative.
         */
        inline const std::vector<std::complex<double> >& GetGZ() const {
            return g_z_;
        }
    private:
        /// Selected slice index for plane tomography.
        std::optional<size_t> slice_index_;
        /// List of seed points.
        std::optional<std::vector<SeedPoint> > seed_points_;
        /// Weight of local estimate in global step (used if seed points are unset).
        double lambda_;
        /// Gradient tolerance w.r.t maximum (used if seed points are unset).
        double gradient_tolerance_;
        /// Mask of the homogeneous regions (used if seed points are unset).
        boost::dynamic_bitset<> mask_;
        /// Filter for the derivative computation.
        FDSavitzkyGolayFilter fd_filter_;
        /// Gradient EPT run flag.
        EPTGradientRun run_mode_;

        /// First estimate of the complex permittivity.
        std::vector<std::complex<double> > epsc_;
        /// Positive derivative of the complex permittivity logarithm.
        std::vector<std::complex<double> > g_plus_;
        /// Longitudinal derivative of the complex permittivity logarithm.
        std::vector<std::complex<double> > g_z_;

        /// Cost functional main term.
        double cost_functional_;
        /// Cost functional regularization term.
        double cost_regularization_;
        /// First estimate flag
        bool thereis_epsc_;

        /// Perform the pixel-by-pixel recovery.
        void LocalRecovery(
            std::array<std::vector<double>,N_DIM> *grad_phi0,
            std::vector<std::complex<double> > *g_plus,
            std::vector<std::complex<double> > *g_z,
            std::vector<std::complex<double> > *theta,
            const int iref);
        /// Perform pixel-by-pixel recovery in a slice.
        void LocalRecoverySlice(
            std::array<std::vector<double>,N_DIM> *grad_phi0,
            std::vector<std::complex<double> > *g_plus,
            std::vector<std::complex<double> > *g_z,
            std::vector<std::complex<double> > *theta,
            const int iref, const int i2);
        /// Estimate the complex permittivity from theta local recovery.
        EPTlibError Theta2Epsc(std::vector<std::complex<double> > *theta,
            const std::array<std::vector<double>,N_DIM> &grad_phi0,
            const std::vector<std::complex<double> > &g_plus,
            const std::vector<std::complex<double> > &g_z);
        /// Select the degrees of freedom.
        void SelectDoF(std::vector<int> *dof, std::vector<int> *ele,
            int *n_dof, const std::vector<int> &locstep);
        /// Solve the global minimisation problem.
        void GlobalMinimisation();
        /// Extract the electric properties from the complex permittivity.
        void ExtractElectricProperties();
};

}  // namespace eptlib

#endif  // EPTLIB_EPT_GRADIENT_H_
