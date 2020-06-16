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
		 * Set the first estimate weight.
		 * 
		 * @return a Success or Unknown error.
		 */
		EPTlibError_t SetLambda(const double lambda);
		/**
		 * Get a reference to the gradient of the reference phase.
		 * 
		 * @return a reference to the gradient of the reference phase.
		 */
		std::array<std::vector<double>,NDIM>& GetGradPhi0();
		/**
		 * Get a const reference to the gradient of the reference phase.
		 * 
		 * @return a const reference to the gradient of the reference phase.
		 */
		const std::array<std::vector<double>,NDIM>& GetGradPhi0() const;
		/**
		 * Get a reference to the deegrees of freedom.
		 * 
		 * @return a reference to the deegrees of freedom.
		 */
		std::vector<int>& GetDof();
		/**
		 * Get a const reference to the deegrees of freedom.
		 * 
		 * @return a const reference to the deegrees of freedom.
		 */
		const std::vector<int>& GetDof() const;
		/**
		 * Get a reference to the positive derivative of the complex permittivity logarithm.
		 * 
		 * @return a reference to the positive derivative of the complex permittivity logarithm.
		 */
		std::vector<std::complex<double> >& GetGPlus();
		/**
		 * Get a const reference to the positive derivative of the complex permittivity logarithm.
		 * 
		 * @return a const reference to the positive derivative of the complex permittivity logarithm.
		 */
		const std::vector<std::complex<double> >& GetGPlus() const;
		/**
		 * Get a reference to the longitudinal derivative of the complex permittivity logarithm.
		 * 
		 * @return a reference to the longitudinal derivative of the complex permittivity logarithm.
		 */
		std::vector<std::complex<double> >& GetGZ();
		/**
		 * Get a const reference to the longitudinal derivative of the complex permittivity logarithm.
		 * 
		 * @return a const reference to the longitudinal derivative of the complex permittivity logarithm.
		 */
		const std::vector<std::complex<double> >& GetGZ() const;
		/**
		 * Get a reference to the first estimate of the complex permittivity.
		 * 
		 * @return a reference to the first estimate of the complex permittivity.
		 */
		std::vector<std::complex<double> >& GetTheta();
		/**
		 * Get a const reference to the first estimate of the complex permittivity.
		 * 
		 * @return a const reference to the first estimate of the complex permittivity.
		 */
		const std::vector<std::complex<double> >& GetTheta() const;
	private:
		/// First estimate weight.
		double lambda_;
		/// Filter for the derivative computation.
		FDSavitzkyGolayFilter fd_filter_;
		/// Gradient of the reference phase.
		std::array<std::vector<double>,NDIM> grad_phi0_;
		/// Deegrees of freedom of `g_plus_', `g_z_' and `theta_'.
		std::vector<int> dof_;
		/// Positive derivative of the complex permittivity logarithm.
		std::vector<std::complex<double> > g_plus_;
		/// Longitudinal derivative of the complex permittivity logarithm.
		std::vector<std::complex<double> > g_z_;
		/// First estimate of the complex permittivity.
		std::vector<std::complex<double> > theta_;
};

}  // namespace eptlib

#endif  // EPTLIB_EPT_GRADIENT_H_
