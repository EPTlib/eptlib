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

#include "eptlib/ept_gradient.h"

#include <cmath>
#include <limits>

#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace eptlib;

// EPTGradient constructor
EPTGradient::
EPTGradient(const double freq, const std::array<int,NDIM> &nn,
	const std::array<double,NDIM> &dd, const int tx_ch,
	const Shape &shape) :
	EPTInterface(freq,nn,dd,tx_ch, 1), lambda_(0.0), fd_filter_(shape), dof_(n_vox_) {
	if (tx_ch_<5) {
		throw std::runtime_error("Impossible to solve the linear system: at least 5 transmit channels are needed.");
	}
	return;
}

// EPTGradient destructor
EPTGradient::
~EPTGradient() {
	return;
}

// EPTGradient run
EPTlibError_t EPTGradient::
Run() {
	if (thereis_tx_sens_.all() && thereis_trx_phase_.all()) {
		thereis_epsr_ = true;
		thereis_sigma_ = true;
		epsr_ = Image<double>(nn_[0],nn_[1],nn_[2]);
		sigma_ = Image<double>(nn_[0],nn_[1],nn_[2]);
	} else {
		return EPTlibError::MissingData;
	}
	// compute the relative sensitivities
	std::vector<std::vector<std::complex<double> > > tx_sens_c(tx_ch_);
	tx_sens_c[0].resize(n_vox_);
	for (int idx = 0; idx<n_vox_; ++idx) {
		tx_sens_c[0][idx] = (*tx_sens_[0])[idx];
	}
	for (int tx = 1; tx<tx_ch_; ++tx) {
		tx_sens_c[tx].resize(n_vox_);
		for (int idx = 0; idx<n_vox_; ++idx) {
			tx_sens_c[tx][idx] = (*tx_sens_[tx])[idx]*std::exp(std::complex<double>(0.0,(*trx_phase_[tx])[idx]-(*trx_phase_[0])[idx]));
		}
	}
	// compute the derivatives
	std::vector<std::array<std::vector<std::complex<double> >,NDIM> > grad(tx_ch_);
	std::vector<std::vector<std::complex<double> > > lapl(tx_ch_);
	for (int tx = 0; tx<tx_ch_; ++tx) {
		for (int d = 0; d<NDIM; ++d) {
			grad[tx][d].resize(n_vox_,std::numeric_limits<double>::quiet_NaN());
			DifferentialOperator_t diff_op = static_cast<DifferentialOperator_t>(d);
			EPTlibError_t error = fd_filter_.Apply(diff_op,grad[tx][d].data(),tx_sens_c[tx].data(),nn_,dd_);
			if (error!=EPTlibError::Success) {
				return error;
			}
		}
		lapl[tx].resize(n_vox_,std::numeric_limits<double>::quiet_NaN());
		DifferentialOperator_t diff_op = DifferentialOperator::Laplacian;
		EPTlibError_t error = fd_filter_.Apply(diff_op,lapl[tx].data(),tx_sens_c[tx].data(),nn_,dd_);
		if (error!=EPTlibError::Success) {
			return error;
		}
	}
	// Select the "degrees of freedom"
	int n_dof = 0;
	int n_dop = 0;
	for (int idx = 0; idx<n_vox_; ++idx) {
		if (lapl[0][idx]==lapl[0][idx]) {
			dof_[idx] = ++n_dof;
		} else {
			dof_[idx] = --n_dop;
		}
	}
	// Write and solve the system voxel-by-voxel
	for (int d = 0; d<NDIM; ++d) {
		grad_phi0_[d].resize(n_vox_,0.0);
	}
	g_plus_.resize(n_dof,0.0);
	g_z_.resize(n_dof,0.0);
	theta_.resize(n_dof,0.0);
	for (int idx = 0; idx<n_vox_; ++idx) {
		if (dof_[idx]<0) {
			continue;
		}
		// fill matrix...
		Eigen::MatrixXd A(2*tx_ch_,9);
		Eigen::VectorXd b(2*tx_ch_);
		for (int tx = 0; tx<tx_ch_; ++tx) {
			int row = 2*tx;
			b[row] = lapl[tx][idx].real();
			b[row+1] = lapl[tx][idx].imag();
			// ...grad phi0
			for (int d = 0; d<NDIM; ++d) {
				A(row,d) = 2.0*grad[tx][d][idx].imag();
				A(row+1,d) = -2.0*grad[tx][d][idx].real();
			}
			// ...xi plus
			A(row,3) = grad[tx][0][idx].real()+grad[tx][1][idx].imag();
			A(row+1,3) = grad[tx][0][idx].imag()-grad[tx][1][idx].real();
			// ...zeta plus
			A(row,4) = grad[tx][1][idx].real()-grad[tx][0][idx].imag();
			A(row+1,4) = grad[tx][1][idx].imag()+grad[tx][0][idx].real();
			// ...xi z
			A(row,5) = grad[tx][2][idx].real();
			A(row+1,5) = grad[tx][2][idx].imag();
			// ...zeta z
			A(row,6) = -grad[tx][2][idx].imag();
			A(row+1,6) = grad[tx][2][idx].real();
			// ...alpha plus
			A(row,7) = tx_sens_c[tx][idx].real();
			A(row+1,7) = tx_sens_c[tx][idx].imag();
			// ...beta plus
			A(row,8) = -tx_sens_c[tx][idx].imag();
			A(row+1,8) = tx_sens_c[tx][idx].real();
		}
		// solve the system
		Eigen::VectorXd x = A.bdcSvd(Eigen::ComputeThinU|Eigen::ComputeThinV).solve(b);
		int idx_dof = dof_[idx]-1;
		for (int d = 0; d<NDIM; ++d) {
			grad_phi0_[d][idx] = x[d];
		}
		g_plus_[idx_dof] = std::complex<double>(x[3],x[4]);
		g_z_[idx_dof] = std::complex<double>(x[5],x[6]);
		theta_[idx_dof] = std::complex<double>(x[7],x[8]);
	}
	// Estimate the complex permittivity from theta
	for (int d = 0; d<NDIM; ++d) {
		std::vector<double> lapl_phi0(n_vox_);
		DifferentialOperator_t diff_op = static_cast<DifferentialOperator_t>(d);
		EPTlibError_t error = fd_filter_.Apply(diff_op,lapl_phi0.data(),grad_phi0_[d].data(),nn_,dd_);
		if (error!=EPTlibError::Success) {
			return error;
		}
		for (int idx = 0; idx<n_vox_; ++idx) {
			if (dof_[idx]<0) {
				continue;
			}
			int idx_dof = dof_[idx]-1;
			theta_[idx_dof] += std::complex<double>(0.0,1.0)*lapl_phi0[idx];
		}
	}
	for (int idx = 0; idx<n_vox_; ++idx) {
		if (dof_[idx]<0) {
			continue;
		}
		int idx_dof = dof_[idx]-1;
		for (int d = 0; d<NDIM; ++d) {
			theta_[idx_dof] -= grad_phi0_[d][idx]*grad_phi0_[d][idx];
		}
		theta_[idx_dof] -= std::complex<double>(0.0,1.0)*grad_phi0_[0][idx]*g_plus_[idx_dof];
		theta_[idx_dof] -= grad_phi0_[1][idx]*g_plus_[idx_dof];
		theta_[idx_dof] -= std::complex<double>(0.0,1.0)*grad_phi0_[2][idx]*g_z_[idx_dof];
		theta_[idx_dof] /= -omega_*omega_*MU0;
	}
	// Invert the gradient through minimisation
	// ...set-up variables
	int vox_nodes = 8;
	std::array<int,NDIM> step{1,nn_[0],nn_[0]*nn_[1]};
	std::vector<int> locstep(vox_nodes,0);
	for (int locid = 0; locid<vox_nodes; ++locid) {
		for (int d = 0; d<NDIM; ++d) {
			locstep[locid] += locid&(1<<d) ? step[d] : 0;
		}
	}
	// ...update the degrees of freedom
	std::vector<int> dof_g(n_vox_);
	int n_dof_g = 0;
	int n_dop_g = 0;
	for (int idx = 0; idx<n_vox_; ++idx) {
		if (dof_[idx]<0) {
			dof_g[idx] = --n_dop_g;
		} else {
			int idx_dof = dof_[idx]-1;
			if (theta_[idx_dof]==theta_[idx_dof]) {
				dof_g[idx] = ++n_dof_g;
			} else {
				dof_g[idx] = --n_dop_g;
			}
		}
	}
	// ...select the nan-free voxels
	int n_max_ele = 1;
	for (int d = 0; d<NDIM; ++d) {
		n_max_ele *= nn_[d]-1;
	}
	std::vector<int> ele(n_max_ele);
	int n_ele = 0;
	{
		int idx_ele = 0;
		for (int id2 = 0; id2<nn_[2]-1; ++id2) {
			for (int id1 = 0; id1<nn_[1]-1; ++id1) {
				for (int id0 = 0; id0<nn_[0]-1; ++id0) {
					int idx = id0 + nn_[0]*(id1 + nn_[1]*id2);
					bool nan_flag = false;
					for (int locid = 0; locid<vox_nodes; ++locid) {
						if (dof_g[idx+locstep[locid]] < 0) {
							nan_flag = true;
							continue;
						}
					}
					if (nan_flag) {
						ele[idx_ele] = -1;
					} else {
						ele[idx_ele] = idx;
						++n_ele;
					}
					++idx_ele;
				}
			}
		}
	}
	// ...build coefficient matrix and forcing term
	Eigen::SparseMatrix<double> A(n_dof_g,n_dof_g);
	std::vector<Eigen::Triplet<double> > A_trip(0);
	Eigen::VectorXcd b(n_dof);
	for (int idof_g = 0; idof_g<n_dof_g; ++idof_g) {
		b[idof_g] = std::complex<double>(0.0,0.0);
	}
	for (int idx_ele = 0; idx_ele<n_max_ele; ++idx_ele) {
		int idx = ele[idx_ele];
		if (idx<0) {
			continue;
		}
		// loop for the equation index...
		for (int i = 0; i<vox_nodes; ++i) {
			int iidx = idx+locstep[i];
			int idof = dof_[iidx]-1;
			int ii = dof_g[iidx]-1;
			// loop for the unknown index...
			for (int j = 0; j<vox_nodes; ++j) {
				int jidx = idx+locstep[j];
				int jdof = dof_[jidx]-1;
				int jj = dof_g[jidx]-1;
				std::array<bool,NDIM> same_side;
				for (int d = 0; d<NDIM; ++d) {
					same_side[d] = (i&(1<<d))==(j&(1<<d));
				}
				// first estimate term
				double Aij = 0.0;
				if (lambda_>0) {
					double tmp = lambda_;
					for (int d = 0; d<NDIM; ++d) {
						tmp *= same_side[d] ? dd_[d]/3.0 : dd_[d]/6.0;
					}
					Aij += tmp;
					// forcing term
					b[ii] += tmp*std::log(theta_[jdof]);
				}
				// derivative terms
				for (int d1 = 0; d1<NDIM; ++d1) {
					double tmp = 1.0;
					for (int d = 0; d<NDIM; ++d) {
						if (d==d1) {
							tmp *= same_side[d] ? 1.0/dd_[d] : -1.0/dd_[d];
						} else {
							tmp *= same_side[d] ? dd_[d]/3.0 : dd_[d]/6.0;
						}
					}
					Aij += tmp;
				}
				// push back the triplet
				A_trip.push_back(Eigen::Triplet<double>(ii,jj,Aij));
				// forcing term
				for (int d1 = 0; d1<NDIM; ++d1) {
					double tmp = 1.0;
					for (int d = 0; d<NDIM; ++d) {
						if (d==d1) {
							tmp *= i&(1<<d) ? 0.5 : -0.5;
						} else {
							tmp *= same_side[d] ? dd_[d]/3.0 : dd_[d]/6.0;
						}
					}
					if (d1 == 0) {
						b[ii] += tmp*g_plus_[jdof];
					} else if (d1 == 1) {
						b[ii] += -std::complex<double>(0.0,tmp)*g_plus_[jdof];
					} else {
						b[ii] += tmp*g_z_[jdof];
					}
				}
			}
		}
	}
	A.setFromTriplets(A_trip.begin(),A_trip.end());
	A.makeCompressed();
	// Solve the linear system
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,Eigen::Upper|Eigen::Lower,Eigen::IncompleteCholesky<double> > solver;
	solver.compute(A);
	Eigen::VectorXcd x = solver.solve(b);
//	std::cout << "\n  #iterations:     " << solver.iterations() << std::endl;
//	std::cout << "  estimated error: " << solver.error()      << std::endl;
	// Extract the electric properties
	for (int idx = 0; idx<n_vox_; ++idx) {
		int idof = dof_g[idx];
		if (idof<0) {
			epsr_[idx] = std::numeric_limits<double>::quiet_NaN();
			sigma_[idx] = std::numeric_limits<double>::quiet_NaN();
		} else {
			std::complex<double> epsc = std::exp(x[idof-1]);
			epsr_[idx] = epsc.real()/EPS0;
			sigma_[idx] = -epsc.imag()*omega_;
		}
	}
	return EPTlibError::Success;
}

EPTlibError_t EPTGradient::
SetLambda(const double lambda) {
	if (lambda < 0.0) {
		return EPTlibError::Unknown;
	}
	lambda_ = lambda;
	return EPTlibError::Success;
}

std::array<std::vector<double>,NDIM>& EPTGradient::
GetGradPhi0() {
	return grad_phi0_;
}
const std::array<std::vector<double>,NDIM>& EPTGradient::
GetGradPhi0() const {
	return grad_phi0_;
}
std::vector<int>& EPTGradient::
GetDof() {
	return dof_;
}
const std::vector<int>& EPTGradient::
GetDof() const {
	return dof_;
}
std::vector<std::complex<double> >& EPTGradient::
GetGPlus() {
	return g_plus_;
}
const std::vector<std::complex<double> >& EPTGradient::
GetGPlus() const {
	return g_plus_;
}
std::vector<std::complex<double> >& EPTGradient::
GetGZ() {
	return g_z_;
}
const std::vector<std::complex<double> >& EPTGradient::
GetGZ() const {
	return g_z_;
}
std::vector<std::complex<double> >& EPTGradient::
GetTheta() {
	return theta_;
}
const std::vector<std::complex<double> >& EPTGradient::
GetTheta() const {
	return theta_;
}
