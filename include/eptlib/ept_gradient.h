/*****************************************************************************
*
*     Program: EPTlib
*     Author: Alessandro Arduino <a.arduino@inrim.it>
*
*  MIT License
*
*  Copyright (c) 2020-2021  Alessandro Arduino
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
	std::array<int,NDIM> ijk;
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
		 * @param freq operative frequency of the MRI scanner.
		 * @param nn number of voxels in each direction.
		 * @param dd voxel sizes in each direction.
		 * @param tx_ch number of transmit channels.
		 * @param shape mask over which apply the finite difference scheme.
		 * @param is_2d 2D assumption flag.
		 * 
		 * The number of Rx channels is fixed equal to one.
		 */
		EPTGradient(const double freq, const std::array<int,NDIM> &nn,
			const std::array<double,NDIM> &dd, const int tx_ch,
			const Shape &shape, const bool is_2d);
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
		 * Set the run mode of the next gradient EPT run.
		 * 
		 * @param run_mode run mode to be set.
		 */
		void SetRunMode(EPTGradientRun run_mode);
		/**
		 * Get the run mode of the next gradient EPT run.
		 */
		EPTGradientRun GetRunMode();
		/**
		 * Set/unset the use of seed points to invert the gradient.
		 * 
		 * @return if it has been set.
		 */
		bool ToggleSeedPoints();
		/**
		 * Get the use seed points flag.
		 * 
		 * @return the seed points flag.
		 */
		bool SeedPointsAreUsed();
		/**
		 * Add a seed point to the list and set their use if not done before.
		 * 
		 * @param seed_point seed point to the added to the list.
		 */
		void AddSeedPoint(const SeedPoint seed_point);
		/**
         * Set the selected plane index for plane tomography.
         * 
         * @return a Success or WrongDataFormat error.
         */
        EPTlibError SelectSlice(const int slice_idx);
		/**
		 * Set the regularization coefficient.
		 * 
		 * @param lambda regularization coefficient.
		 * 
		 * @return a Success or WrongDataFormat error.
		 */
		EPTlibError SetRegularizationCoefficient(const double lambda);
		/**
		 * Set the gradient tolerance for estimating the mask of homogeneous regions.
		 * 
		 * @return a Success or WrongDataFormat error.
		 */
		EPTlibError SetGradientTolerance(const double gradient_tolerance);
		/**
		 * Get the homogeneous regions mask.
		 * 
		 * @return a constant reference to the mask.
		 */
		const boost::dynamic_bitset<>& GetMask();
		/**
		 * Get the cost functional.
		 * 
		 * @return the cost functional.
		 */
		double GetCostFunctional();
		/**
		 * Get the cost regularization.
		 * 
		 * @return the cost regularization.
		 */
		double GetCostRegularization();
		/**
		 * Get the complex permittivity.
		 * 
		 * @param[out] epsc pointer to the complex permittivity destination.
		 * 
		 * @return a Success or MissingData error.
		 */
		EPTlibError GetEpsC(std::vector<std::complex<double> > *epsc);
		/**
		 * Get the positive derivative of the complex permittivity logarithm.
		 * 
		 * @param[out] g_plus pointer to the positive derivative destination.
		 * 
		 * @return a Success or MissingData error.
		 */
		EPTlibError GetGPlus(std::vector<std::complex<double> > *g_plus);
		/**
		 * Get the longitudinal derivative of the complex permittivity logarithm.
		 * 
		 * @param[out] g_z pointer to the longitudinal derivative destination.
		 * 
		 * @return a Success or MissingData error.
		 */
		EPTlibError GetGZ(std::vector<std::complex<double> > *g_z);
	private:
		/// 2D assumption flag.
		const bool is_2d_;
		/// Seed points flag.
		bool use_seed_points_;
        /// Selected plane index (used if 2D is set).
        int plane_idx_;
		/// First estimate weight (used if seed points are unset).
		double lambda_;
		/// Gradient tolerance w.r.t maximum (used if seed points are unset).
		double gradient_tolerance_;
		/// Mask of the homogeneous regions (used if seed points are unset).
		boost::dynamic_bitset<> mask_;
		/// List of seed points (used if seed points is set).
		std::vector<SeedPoint> seed_points_;
		/// Filter for the derivative computation.
		FDSavitzkyGolayFilter fd_filter_;
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
		/// Gradient EPT run flag.
		EPTGradientRun run_mode_;

		// Auxiliary methods
		/// Perform the pixel-by-pixel recovery.
		void LocalRecovery(
			std::array<std::vector<double>,NDIM> *grad_phi0,
			std::vector<std::complex<double> > *g_plus,
			std::vector<std::complex<double> > *g_z,
			std::vector<std::complex<double> > *theta,
			const int iref);
		/// Perform pixel-by-pixel recovery in a slice.
		void LocalRecoverySlice(
			std::array<std::vector<double>,NDIM> *grad_phi0,
			std::vector<std::complex<double> > *g_plus,
			std::vector<std::complex<double> > *g_z,
			std::vector<std::complex<double> > *theta,
			const int iref, const int i2);
		/// Estimate the complex permittivity from theta local recovery.
		EPTlibError Theta2Epsc(std::vector<std::complex<double> > *theta,
			const std::array<std::vector<double>,NDIM> &grad_phi0,
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
