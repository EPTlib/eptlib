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

#ifndef EPTLIB_EPT_GRADIENT_H_
#define EPTLIB_EPT_GRADIENT_H_

#include "eptlib/ept_interface.h"

#include <array>
#include <complex>
#include <vector>

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
		 * 
		 * The number of Rx channels is fixed equal to one.
		 */
		EPTGradient(const double freq, const std::array<int,NDIM> &nn,
			const std::array<double,NDIM> &dd, const int tx_ch,
			const Shape &shape);
		/**
		 * Virtual destructor.
		 */
		virtual ~EPTGradient();
		/**
		 * Perform the gradient EPT.
		 * 
		 * @return an error index about the state of the tomography.
		 */
		virtual EPTlibError_t Run() override;
		/**
         * Set the selected plane index for plane tomography.
         * 
         * @return a Success or WrongDataFormat error.
         */
        EPTlibError_t SelectPlane(const int plane_idx);
		/**
		 * Set the two-dimensional assumption.
		 * 
		 * @return a Success or Unknown error.
		 */
		EPTlibError_t Set2D();
		/**
		 * Unset the two-dimensional assumption.
		 * 
		 * @return a Success or Unknown error.
		 */
		EPTlibError_t Unset2D();
		/**
		 * Set the first estimate weight.
		 * 
		 * @return a Success or Unknown error.
		 */
		EPTlibError_t SetLambda(const double lambda);
		/**
		 * Set the L-curve method to determine the optimal lambda.
		 * 
		 * @return a Success or Unknown error.
		 */
		EPTlibError_t SetLCurve();
		/**
		 * Unset the L-curve method to determine the optimal lambda.
		 * 
		 * @return a Success or Unknown error.
		 */
		EPTlibError_t UnsetLCurve();
		/*
		 * Set/unset the median to average the local system results.
		 * 
		 * @return a Success or Unknown error.
		 */
		EPTlibError_t ToggleAverageWithMedian();
		/*
		 * Add a seed point to the list and set their use if not done before.
		 * 
		 * @return a Success or Unknown error.
		 */
		EPTlibError_t AddSeedPoint(const SeedPoint seed_point);
	private:
        /// Selected plane index (for 2D assumption only).
        int plane_idx_;
		/// 2D assumption flag.
		bool is_2d_;
		/// First estimate weight.
		double lambda_;
		/// Use l-curve method.
		bool use_lcurve_;
		/// Use the median to average the local system results.
		bool average_with_median_;
		/// List of seed points.
		std::vector<SeedPoint> seed_points_;
		/// Use the seed points.
		bool use_seed_points_;
		/// Filter for the derivative computation.
		FDSavitzkyGolayFilter fd_filter_;
		/// Gradient of the reference phase.
		std::vector<std::vector<double> > gradx_phi0_;
		std::vector<std::vector<double> > grady_phi0_;
		std::vector<std::vector<double> > gradz_phi0_;
		/// Positive derivative of the complex permittivity logarithm.
		std::vector<std::vector<std::complex<double> > > g_plus_;
		/// Longitudinal derivative of the complex permittivity logarithm.
		std::vector<std::vector<std::complex<double> > > g_z_;
		/// First estimate of the complex permittivity.
		std::vector<std::vector<std::complex<double> > > theta_;
		// Auxiliary methods
		// Perform the pixel-by-pixel recovery.
		EPTlibError_t LocalRecovery(std::vector<double> *gradx_phi0,
			std::vector<double> *grady_phi0, std::vector<double> *gradz_phi0,
			std::vector<std::complex<double> > *g_plus,
			std::vector<std::complex<double> > *g_z,
			std::vector<std::complex<double> > *theta,
			const int iref);
		/// Perform pixel-by-pixel recovery in a slice.
		EPTlibError_t LocalRecoverySlice(std::vector<double> *gradx_phi0,
			std::vector<double> *grady_phi0, std::vector<double> *gradz_phi0,
			std::vector<std::complex<double> > *g_plus,
			std::vector<std::complex<double> > *g_z,
			std::vector<std::complex<double> > *theta,
			const int iref, const int i2);
		/// Estimate the complex permittivity from theta local recovery.
		EPTlibError_t Theta2Epsc(std::vector<std::complex<double> > *epsc,
			const std::array<std::vector<double>*,NDIM> &grad_phi0,
			const std::vector<std::complex<double> > &g_plus,
			const std::vector<std::complex<double> > &g_z,
			const std::vector<std::complex<double> > &theta);
		/// Average a quantity to improve its quality.
		EPTlibError_t AverageQuantity(std::vector<std::complex<double> > *avg,
			const std::vector<std::vector<std::complex<double> > > &src);
		/// Select the degrees of freedom.
		EPTlibError_t SelectDoF(std::vector<int> *dof, std::vector<int> *ele,
			int *n_dof, const std::vector<std::complex<double> > &epsc,
			const std::vector<int> &locstep);
		/// Solve the global minimisation problem.
		EPTlibError_t GlobalMinimisation(
			std::vector<std::complex<double> > *epsc,
			const std::vector<std::complex<double> > &g_plus,
			const std::vector<std::complex<double> > &g_z);
		/// Extract the electric properties of a slice.
		EPTlibError_t ExtractElectricPropertiesSlice(const int i2,
			const int out_i2, const std::vector<std::complex<double> > &epsc);
};

}  // namespace eptlib

#endif  // EPTLIB_EPT_GRADIENT_H_
