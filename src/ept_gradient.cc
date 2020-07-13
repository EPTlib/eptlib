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
#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace eptlib;

// EPTGradient constructor
EPTGradient::
EPTGradient(const double freq, const std::array<int,NDIM> &nn,
	const std::array<double,NDIM> &dd, const int tx_ch,
	const Shape &shape) :
	EPTInterface(freq,nn,dd,tx_ch, 1),
	is_2d_(false), lambda_(0.0), use_lcurve_(false), average_with_median_(false),
	use_seed_points_(false), seed_points_(0),
	fd_filter_(shape), gradx_phi0_(tx_ch_), grady_phi0_(tx_ch_),
	gradz_phi0_(tx_ch_), g_plus_(tx_ch_), g_z_(tx_ch_), theta_(tx_ch_) {
	if (tx_ch_<5) {
		throw std::runtime_error("Impossible to solve the linear system: at least 5 transmit channels are needed.");
	}
	plane_idx_ = nn_[2]/2;
	for (int itx = 0; itx<tx_ch_; ++itx) {
		gradx_phi0_[itx].resize(n_vox_,std::numeric_limits<double>::quiet_NaN());
		grady_phi0_[itx].resize(n_vox_,std::numeric_limits<double>::quiet_NaN());
		gradz_phi0_[itx].resize(n_vox_,std::numeric_limits<double>::quiet_NaN());
		g_plus_[itx].resize(n_vox_,std::numeric_limits<double>::quiet_NaN());
		g_z_[itx].resize(n_vox_,std::numeric_limits<double>::quiet_NaN());
		theta_[itx].resize(n_vox_,std::numeric_limits<double>::quiet_NaN());
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
		if (is_2d_) {
			epsr_ = Image<double>(nn_[0],nn_[1]);
			sigma_ = Image<double>(nn_[0],nn_[1]);
		} else {
			epsr_ = Image<double>(nn_[0],nn_[1],nn_[2]);
			sigma_ = Image<double>(nn_[0],nn_[1],nn_[2]);
		}
	} else {
		return EPTlibError::MissingData;
	}
	// Perform the local recovery using each transmit channel as reference
	for (int iref = 0; iref<tx_ch_; ++iref) {
		LocalRecovery(&(gradx_phi0_[iref]),&(grady_phi0_[iref]),
		&(gradz_phi0_[iref]),&(g_plus_[iref]),&(g_z_[iref]),&(theta_[iref]),
		iref);
	}
	// Estimate the complex permittivity from theta for each reference
	std::vector<std::vector<std::complex<double> > > epsc(tx_ch_);
	for (int iref = 0; iref<tx_ch_; ++iref) {
		epsc[iref].resize(n_vox_);
		std::array<std::vector<double>*,NDIM> grad_phi0;
		grad_phi0[0] = &(gradx_phi0_[iref]);
		grad_phi0[1] = &(grady_phi0_[iref]);
		grad_phi0[2] = &(gradz_phi0_[iref]);
		Theta2Epsc(&(epsc[iref]),grad_phi0,g_plus_[iref],g_z_[iref],theta_[iref]);
	}
	// Average the epsc, g_plus and g_z quantities
	std::vector<std::complex<double> > epsc_avg(n_vox_);
	std::vector<std::complex<double> > g_plus_avg(n_vox_);
	std::vector<std::complex<double> > g_z_avg(0);
	AverageQuantity(&epsc_avg,epsc);
	AverageQuantity(&g_plus_avg,g_plus_);
	if (!is_2d_) {
		g_z_avg.resize(n_vox_);
		AverageQuantity(&g_z_avg,g_z_);
	}
	// Optimise the estimate using the gradient information
	GlobalMinimisation(&epsc_avg, g_plus_avg,g_z_avg);
	// Extract the electric properties
	if (is_2d_) {
		ExtractElectricPropertiesSlice(plane_idx_,0,epsc_avg);
	} else {
		for (int i2 = 0; i2<nn_[2]; ++i2) {
			ExtractElectricPropertiesSlice(i2,i2,epsc_avg);
		}
	}
	return EPTlibError::Success;
}

/////////////////////////////////////////////////////////////////////////////
// EPTGradient LOCAL RECOVERY ///////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
namespace { //details
	static constexpr int col_theta = 0;
	static constexpr int col_gplus = 2;
	static constexpr int col_phi = 4;
	static constexpr int col_gz = 7;

	// Fill the local matrix for pixel-by-pixel recovery.
	EPTlibError_t FillLocalMatrix(Eigen::MatrixXd *A, Eigen::VectorXd *b,
		const std::vector<std::vector<std::complex<double> > > &field_crop,
		const FDSavitzkyGolayFilter &fd_filter,
		const std::array<double,NDIM> &dd, const int tx_ch, const bool is_2d) {
		int mid_vox = fd_filter.GetShape().GetBoxVolume()/2;
		int n_col = is_2d?6:9;
		*A = Eigen::MatrixXd::Zero(2*tx_ch,n_col);
		*b = Eigen::VectorXd::Zero(2*tx_ch);
		// for each transmit channel...
		for (int itx = 0; itx<tx_ch; ++itx) {
			int row = 2*itx;
			// ...laplacian
			std::complex<double> tmp = fd_filter.Laplacian(field_crop[itx],dd);
			(*b)[row] = tmp.real();
			(*b)[row+1] = tmp.imag();
			// ...grad_x
			tmp = fd_filter.FirstOrder(0,field_crop[itx],dd);
			(*A)(row,col_phi) = 2.0*tmp.imag();
			(*A)(row+1,col_phi) = -2.0*tmp.real();
			(*A)(row,col_gplus) = tmp.real();
			(*A)(row+1,col_gplus) = tmp.imag();
			(*A)(row,col_gplus+1) = -tmp.imag();
			(*A)(row+1,col_gplus+1) = tmp.real();
			// ...grad_y
			tmp = fd_filter.FirstOrder(1,field_crop[itx],dd);
			(*A)(row,col_phi+1) = 2.0*tmp.imag();
			(*A)(row+1,col_phi+1) = -2.0*tmp.real();
			(*A)(row,col_gplus) += tmp.imag();
			(*A)(row+1,col_gplus) += -tmp.real();
			(*A)(row,col_gplus+1) += tmp.real();
			(*A)(row+1,col_gplus+1) += tmp.imag();
			// ...grad_z
			if (!is_2d) {
				tmp = fd_filter.FirstOrder(2,field_crop[itx],dd);
				(*A)(row,col_phi+2) = 2.0*tmp.imag();
				(*A)(row+1,col_phi+2) = -2.0*tmp.real();
				(*A)(row,col_gz) = tmp.real();
				(*A)(row+1,col_gz) = tmp.imag();
				(*A)(row,col_gz+1) = -tmp.imag();
				(*A)(row+1,col_gz+1) = tmp.real();
			}
			// ...value
			tmp = field_crop[itx][mid_vox];
			(*A)(row,col_theta) = tmp.real();
			(*A)(row+1,col_theta) = tmp.imag();
			(*A)(row,col_theta+1) = -tmp.imag();
			(*A)(row+1,col_theta+1) = tmp.real();
		}
		return EPTlibError::Success;
	}
	// Solve the local system for pixel-by-pixel recovery.
	EPTlibError_t SolveLocalSystem(double *gradx_phi0, double *grady_phi0,
		double *gradz_phi0, std::complex<double> *g_plus,
		std::complex<double> *g_z, std::complex<double> *theta,
		Eigen::MatrixXd A, Eigen::VectorXd b, const bool is_2d) {
		// solve the system
		b = A.transpose()*b;
		A = A.transpose()*A;
		Eigen::VectorXd x = A.ldlt().solve(b);
		// distribute the result
		*theta = std::complex<double>(x[col_theta],x[col_theta+1]);
		*g_plus = std::complex<double>(x[col_gplus],x[col_gplus+1]);
		*gradx_phi0 = x[col_phi];
		*grady_phi0 = x[col_phi+1];
		if (!is_2d) {
			*gradz_phi0 = x[col_phi+2];
			*g_z = std::complex<double>(x[col_gz],x[col_gz+1]);
		}
		return EPTlibError::Success;
	}
}  // namespace
// Perform the pixel-by-pixel recovery
EPTlibError_t EPTGradient::
LocalRecovery(std::vector<double> *gradx_phi0,
	std::vector<double> *grady_phi0,
	std::vector<double> *gradz_phi0,
	std::vector<std::complex<double> > *g_plus,
	std::vector<std::complex<double> > *g_z,
	std::vector<std::complex<double> > *theta,
	const int iref) {
	int r2 = fd_filter_.GetShape().GetSize()[2]/2;
	if (is_2d_) {
		for (int i2 = plane_idx_-r2; i2<=plane_idx_+r2; ++i2) {
			LocalRecoverySlice(gradx_phi0,grady_phi0,gradz_phi0,g_plus,g_z,
				theta,iref,i2);
		}
	} else {
		for (int i2 = r2; i2<nn_[2]-r2; ++i2) {
			LocalRecoverySlice(gradx_phi0,grady_phi0,gradz_phi0,g_plus,g_z,
				theta,iref,i2);
		}
	}
	return EPTlibError::Success;
}
// Perform pixel-by-pixel recovery in a slice.
EPTlibError_t EPTGradient::
LocalRecoverySlice(std::vector<double> *gradx_phi0,
	std::vector<double> *grady_phi0,
	std::vector<double> *gradz_phi0,
	std::vector<std::complex<double> > *g_plus,
	std::vector<std::complex<double> > *g_z,
	std::vector<std::complex<double> > *theta,
	const int iref, const int i2) {
	std::array<int,NDIM> rr = fd_filter_.GetShape().GetSize();
	for (int d = 0; d<NDIM; ++d) {
		rr[d] /= 2;
	}
	std::array<int,NDIM> ii;
	std::array<int,NDIM> m_inc;
	int m_vox = fd_filter_.GetShape().GetBoxVolume();
	m_inc[0] = 1;
	m_inc[1] = nn_[0]-fd_filter_.GetShape().GetSize()[0];
	m_inc[2] = nn_[0]*(nn_[1]-fd_filter_.GetShape().GetSize()[1]);
	ii[2] = i2;
	for (ii[1] = rr[1]; ii[1]<nn_[1]-rr[1]; ++ii[1]) {
		for (ii[0] = rr[0]; ii[0]<nn_[0]-rr[0]; ++ii[0]) {
			std::array<int,NDIM> ii_l;
			std::vector<std::vector<std::complex<double> > > field_crop(tx_ch_);
			for (int itx = 0; itx<tx_ch_; ++itx) {
				field_crop[itx].resize(m_vox,0.0);
			}
			// inner loop over the kernel voxels
			int idx_l = 0;
			int idx_g = ii[0]-rr[0] + nn_[0]*(ii[1]-rr[1] + nn_[1]*(ii[2]-rr[2]));
			for (ii_l[2] = -rr[2]; ii_l[2]<=rr[2]; ++ii_l[2]) {
				for (ii_l[1] = -rr[1]; ii_l[1]<=rr[1]; ++ii_l[1]) {
					for (ii_l[0] = -rr[0]; ii_l[0]<=rr[0]; ++ii_l[0]) {
						if (fd_filter_.GetShape()[idx_l]) {
							// set the field value
							for (int itx = 0; itx<tx_ch_; ++itx) {
								field_crop[itx][idx_l] = (*tx_sens_[itx])[idx_g]*
									std::exp(std::complex<double>(0.0,(*trx_phase_[itx])[idx_g]-
									(*trx_phase_[iref])[idx_g]));
							}
						}
						++idx_l;
						idx_g += m_inc[0];
					}
					idx_g += m_inc[1];
				}
				idx_g += m_inc[2];
			}
			// fill the matrix
			Eigen::MatrixXd A;
			Eigen::VectorXd b;
			::FillLocalMatrix(&A,&b, field_crop,fd_filter_,dd_,tx_ch_,is_2d_);
			// solve the system
			int idx = ii[0] + nn_[0]*(ii[1] + nn_[1]*ii[2]);
			::SolveLocalSystem(&((*gradx_phi0)[idx]),&((*grady_phi0)[idx]),
				&((*gradz_phi0)[idx]),&((*g_plus)[idx]),&((*g_z)[idx]),
				&((*theta)[idx]), A,b,is_2d_);
		}
	}
	return EPTlibError::Success;
}

/////////////////////////////////////////////////////////////////////////////
// EPTGradient THETA TO COMPLEX PERMITTIVITY ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
// Estimate the complex permittivity from theta local recovery.
EPTlibError_t EPTGradient::
Theta2Epsc(std::vector<std::complex<double> > *epsc,
	const std::array<std::vector<double>*,NDIM> &grad_phi0,
	const std::vector<std::complex<double> > &g_plus,
	const std::vector<std::complex<double> > &g_z,
	const std::vector<std::complex<double> > &theta) {
	int n_dim = is_2d_?2:3;
	std::copy(theta.begin(),theta.end(),epsc->begin());
	// grad_phi0 laplacian...
	for (int d = 0; d<n_dim; ++d) {
		std::vector<double> lapl_phi0(n_vox_);
		DifferentialOperator_t diff_op = static_cast<DifferentialOperator_t>(d);
		EPTlibError_t error = fd_filter_.Apply(diff_op,lapl_phi0.data(),
			grad_phi0[d]->data(),nn_,dd_);
		if (error!=EPTlibError::Success) {
			return error;
		}
		std::transform(epsc->begin(),epsc->end(),lapl_phi0.begin(),epsc->begin(),
			[](const std::complex<double> &a,const double &b)->std::complex<double>{
				return a + std::complex<double>(0.0,b);
			});
	}
	// grad_phi0 magnitude...
	for (int d = 0; d<n_dim; ++d) {
		std::vector<double> tmp(n_vox_);
		std::transform(grad_phi0[d]->begin(),grad_phi0[d]->end(),tmp.begin(),
			[](const double &a)->double{
				return a*a;
			});
		std::transform(epsc->begin(),epsc->end(),tmp.begin(),epsc->begin(),
		[](const std::complex<double> &a,const double &b)->std::complex<double>{
			return a-b;
		});
	}
	// grad_phi0, g inner product...
	std::vector<std::complex<double> > tmp(n_vox_);
	// ...component 0
	std::transform(grad_phi0[0]->begin(),grad_phi0[0]->end(),g_plus.begin(),tmp.begin(),
		[](const double &a,const std::complex<double> &b)->std::complex<double>{
			return std::complex<double>(0.0,a)*b;
		});
	std::transform(epsc->begin(),epsc->end(),tmp.begin(),epsc->begin(),
		std::minus<std::complex<double> >());
	// ...component 1
	std::transform(grad_phi0[1]->begin(),grad_phi0[1]->end(),g_plus.begin(),tmp.begin(),
		[](const double &a,const std::complex<double> &b)->std::complex<double>{
			return a*b;
		});
	std::transform(epsc->begin(),epsc->end(),tmp.begin(),epsc->begin(),
		std::minus<std::complex<double> >());
	if (!is_2d_) {
		// ...component 2
		std::transform(grad_phi0[2]->begin(),grad_phi0[2]->end(),g_z.begin(),tmp.begin(),
		[](const double &a,const std::complex<double> &b)->std::complex<double>{
			return std::complex<double>(0.0,a)*b;
		});
		std::transform(epsc->begin(),epsc->end(),tmp.begin(),epsc->begin(),
			std::minus<std::complex<double> >());
	}
	// multiplicative factor...
	std::transform(epsc->begin(),epsc->end(),epsc->begin(),
		[this](const std::complex<double> &a)->std::complex<double>{
			return -a/omega_/omega_/MU0;
		});
	return EPTlibError::Success;
}
// Average a quantity to improve its quality
EPTlibError_t EPTGradient::
AverageQuantity(std::vector<std::complex<double> > *avg,
	const std::vector<std::vector<std::complex<double> > > &src) {
	int med = tx_ch_/2;
	bool is_odd = tx_ch_%2;
	for (int idx = 0; idx<n_vox_; ++idx) {
		if (average_with_median_) {
			// select the elements
			std::vector<double> real(tx_ch_);
			std::vector<double> imag(tx_ch_);
			for (int iref = 0; iref<tx_ch_; ++iref) {
				real[iref] = src[iref][idx].real();
				imag[iref] = src[iref][idx].imag();
			}
			// find the median
			std::nth_element(real.begin(),real.begin()+med,real.end());
			std::nth_element(imag.begin(),imag.begin()+med,imag.end());
			(*avg)[idx] = std::complex<double>(real[med],imag[med]);
			if (!is_odd) {
				auto it_real = std::max_element(real.begin(),real.begin()+med);
				auto it_imag = std::max_element(imag.begin(),imag.begin()+med);
				(*avg)[idx] += std::complex<double>(*it_real,*it_imag);
				(*avg)[idx] /= 2.0;
			}
		} else {
			(*avg)[idx] = 0.0;
			double total = 0.0;
			for (int iref = 0; iref<tx_ch_; ++iref) {
				(*avg)[idx] += src[iref][idx]*(*tx_sens_[iref])[idx];
				total += (*tx_sens_[iref])[idx];
			}
			(*avg)[idx] /= total;
		}
	}
	return EPTlibError::Success;
}

/////////////////////////////////////////////////////////////////////////////
// EPTGradient GLOBAL OPTIMISATION //////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
namespace { //details
	static constexpr int vox_nodes = 8;
	static constexpr int quad_nodes = 4;

	// Select the degrees of freedom.
	EPTlibError_t FillDoF(std::vector<int> *dof, int *idx_dof, int *n_dof,
		int *n_dop, const int i2, const std::array<int,NDIM> &nn,
		const std::vector<std::complex<double> > &epsc) {
		int nxy = nn[0]*nn[1];
		int idx = nxy*i2;
		for (int out_idx = 0; out_idx<nxy; ++out_idx) {
			if (epsc[idx]==epsc[idx]) {
				(*dof)[(*idx_dof)++] = ++(*n_dof);
			} else {
				(*dof)[(*idx_dof)++] = --(*n_dop);
			}
			++idx;
		}
		return EPTlibError::Success;
	}
	// Select the NaN-free elements.
	EPTlibError_t NaNFreeEle(std::vector<int> *ele, int *idx_ele, const int i2,
		const std::array<int,NDIM> &nn, const std::vector<int> &locstep,
		const std::vector<int> &dof, const bool is_2d) {
		int ele_nodes = is_2d?quad_nodes:vox_nodes;
		int nxy = nn[0]*nn[1];
		int idx_out = is_2d?0:nxy*i2;
		for (int i1 = 0; i1<nn[1]-1; ++i1) {
			for (int i0 = 0; i0<nn[0]-1; ++i0) {
				bool thereis_nan = false;
				for (int locid = 0; locid<ele_nodes; ++locid) {
					if (dof[idx_out+locstep[locid]]<0) {
						thereis_nan = true;
						break;
					}
				}
				if (thereis_nan) {
					(*ele)[(*idx_ele)++] = -1;
				} else {
					(*ele)[(*idx_ele)++] = idx_out;
				}
				++idx_out;
			}
			++idx_out;
		}
		return EPTlibError::Success;
	}
	// Clean dof by counting to how many elements they belong.
	EPTlibError_t CleanDoF(std::vector<int> *dof, int *n_dof, int *n_dop,
		const std::vector<int> &ele, const std::vector<int> &locstep,
		const bool is_2d) {
		int ele_nodes = is_2d?quad_nodes:vox_nodes;
		// count the elements
		std::vector<int> howmany_ele(dof->size(),0);
		for (int idx_ele = 0; idx_ele<ele.size(); ++idx_ele) {
			int idx = ele[idx_ele];
			if (idx>=0) {
				for (int locid = 0; locid<ele_nodes; ++locid) {
					++howmany_ele[idx+locstep[locid]];
				}
			}
		}
		// clean the dof
		*n_dof = 0;
		*n_dop = 0;
		for (int idx = 0; idx<dof->size(); ++idx) {
			if ((*dof)[idx]>0 && howmany_ele[idx]>0) {
				(*dof)[idx] = ++(*n_dof);
			} else {
				(*dof)[idx] = --(*n_dop);
			}
		}
		return EPTlibError::Success;
	}

	// Associate each seed point to its dof.
	EPTlibError_t SeedPointsToDoF(std::vector<std::complex<double> > *dir,
		std::vector<int> *dof, int *n_dof,
		const std::vector<SeedPoint> &seed_points, const std::array<int,NDIM> &nn,
		const int i2, const double omega, const bool is_2d) {
		// for each seed find the corresponding idx
		auto n_out = dof->size();
		std::vector<int> idx2sp(n_out,-1);
		for (int idx_sp = 0; idx_sp<seed_points.size(); ++idx_sp) {
			SeedPoint sp = seed_points[idx_sp];
			if (is_2d && sp.ijk[2]!=i2) {
				continue;
			}
			int idx = sp.ijk[0] + nn[0]*sp.ijk[1];
			if (!is_2d) {
				idx += nn[0]*nn[1]*sp.ijk[2];
			}
			std::cout << idx << ": " << std::flush;
			idx2sp[idx] = idx_sp;
			(*dir)[idx_sp] = std::log(std::complex<double>(sp.epsr*EPS0,-sp.sigma/omega));

			std::cout << idx_sp << ", " << (*dir)[idx_sp] << std::endl;
		}
		// correct the dof
		*n_dof = 0;
		for (int idx = 0; idx<n_out; ++idx) {
			if ((*dof)[idx]>0) {
				if (idx2sp[idx]<0) {
					(*dof)[idx] = ++(*n_dof);
				} else {
					(*dof)[idx] = -idx2sp[idx]-1;
				}
			} else {
				(*dof)[idx] = 0;
			}
		}
		return EPTlibError::Success;
	}
	
	// Fill the global matrix for minimisation.
	EPTlibError_t FillGlobalMatrix(Eigen::SparseMatrix<std::complex<double> > *A,
		Eigen::SparseMatrix<std::complex<double> > *B, Eigen::VectorXcd *b,
		Eigen::VectorXcd *l, const std::vector<int> &dof, const int n_dof,
		const std::vector<int> &ele, const std::vector<int> &locstep,
		const std::vector<std::complex<double> > &epsc,
		const std::vector<std::complex<double> > &g_plus,
		const std::vector<std::complex<double> > &g_z,
		const std::array<double,NDIM> &dd, const int idx0,
		const bool is_2d, const bool use_seed_points_,
		const std::vector<std::complex<double> > &dir) {
		int ele_nodes = is_2d?quad_nodes:vox_nodes;
		int n_dim = is_2d?2:3;
		std::vector<Eigen::Triplet<std::complex<double> > > A_trip(0);
		std::vector<Eigen::Triplet<std::complex<double> > > B_trip(0);
		*b = Eigen::VectorXcd::Zero(n_dof);
		*l = Eigen::VectorXcd::Zero(n_dof);
		// loop the elements
		for (int idx_ele = 0; idx_ele<ele.size(); ++idx_ele) {
			int idx_out = ele[idx_ele];
			if (idx_out<0) {
				continue;
			}
			// loop the equation index
			for (int i = 0; i<ele_nodes; ++i) {
				int iidx_out = idx_out+locstep[i];
				int iidx = idx0+iidx_out;
				int ii = dof[iidx_out];
				if (ii>0) {
					--ii;
					// loop the unknown index
					for (int j = 0; j<ele_nodes; ++j) {
						int jidx_out = idx_out+locstep[j];
						int jidx = idx0+jidx_out;
						int jj = dof[jidx_out];
						std::array<bool,NDIM> same_side;
						for (int d = 0; d<NDIM; ++d) {
							same_side[d] = (i&(1<<d))==(j&(1<<d));
						}
						// matrix B and vector l
						if (!use_seed_points_) {
							std::complex<double> Bij = 1.0;
							for (int d = 0; d<n_dim; ++d) {
								Bij *= same_side[d] ? dd[d]/3.0 : dd[d]/6.0;
							}
							B_trip.push_back(Eigen::Triplet<std::complex<double> >(ii,jj-1,Bij));
							(*l)[ii] += Bij*std::log(epsc[jidx]);
						}
						// matrix A
						std::complex<double> Aij = 0.0;
						for (int d1 = 0; d1<n_dim; ++d1) {
							double tmp = 1.0;
							for (int d = 0; d<n_dim; ++d) {
								if (d==d1) {
									tmp *= same_side[d] ? 1.0/dd[d] : -1.0/dd[d];
								} else {
									tmp *= same_side[d] ? dd[d]/3.0 : dd[d]/6.0;
								}
							}
							Aij += tmp;
						}
						{
							double tmp = 1.0;
							tmp *= i&(1<<0) ? 0.5 : -0.5;
							tmp *= j&(1<<1) ? 0.5 : -0.5;
							tmp *= same_side[2] ? dd[2]/3.0 : dd[2]/6.0;
							Aij += std::complex<double>(0.0,tmp);
						}
						{
							double tmp = -1.0;
							tmp *= j&(1<<0) ? 0.5 : -0.5;
							tmp *= i&(1<<1) ? 0.5 : -0.5;
							tmp *= same_side[2] ? dd[2]/3.0 : dd[2]/6.0;
							Aij += std::complex<double>(0.0,tmp);
						}
						if (jj>0) {
							A_trip.push_back(Eigen::Triplet<std::complex<double> >(ii,jj-1,Aij));
						} else {
							(*b)[ii] += -Aij*dir[-jj-1];
						}
						// vector b
						for (int d1 = 0; d1<n_dim; ++d1) {
							double tmp = 1.0;
							for (int d = 0; d<n_dim; ++d) {
								if (d==d1) {
									tmp *= i&(1<<d) ? 0.5 : -0.5;
								} else {
									tmp *= same_side[d] ? dd[d]/3.0 : dd[d]/6.0;
								}
							}
							if (d1==0) {
								(*b)[ii] += tmp*g_plus[jidx];
							} else if (d1==1) {
								(*b)[ii] += -std::complex<double>(0.0,tmp)*g_plus[jidx];
							} else {
								(*b)[ii] += tmp*g_z[jidx];
							}
						}
					}
				}
			}
		}
		(*A).setFromTriplets(A_trip.begin(),A_trip.end());
		(*A).makeCompressed();
		(*B).setFromTriplets(B_trip.begin(),B_trip.end());
		(*B).makeCompressed();
		return EPTlibError::Success;
	}
	// Solve the global system for minimisation.
//	EPTlibError_t SolveGlobalSystem(Eigen::VectorXcd *x,
//		const Eigen::SparseMatrix<std::complex<double> > &A,
//		const Eigen::SparseMatrix<double> &B,
//		const Eigen::VectorXcd &b, const Eigen::VectorXcd &l,
//		const double lambda) {
	EPTlibError_t SolveGlobalSystem(Eigen::VectorXcd *x,
		const Eigen::SparseMatrix<std::complex<double> > A,
		const Eigen::VectorXcd b) {
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<std::complex<double> > > solver;
		solver.compute(A);
		*x = solver.solve(b);
		return EPTlibError::Success;
	}
}  // namespace
// Select the degree of freedom.
EPTlibError_t EPTGradient::
SelectDoF(std::vector<int> *dof, std::vector<int> *ele,
	int *n_dof, const std::vector<std::complex<double> > &epsc,
	const std::vector<int> &locstep) {
	// select the dof
	int idx_dof = 0;
	int n_dop = 0;
	*n_dof = 0;
	if (is_2d_) {
		::FillDoF(dof,&idx_dof,n_dof,&n_dop, plane_idx_,nn_,epsc);
	} else {
		for (int i2 = 0; i2<nn_[2]; ++i2) {
			::FillDoF(dof,&idx_dof,n_dof,&n_dop, i2,nn_,epsc);
		}
	}
	// determine the nan-free elements
	int idx_ele = 0;
	if (is_2d_) {
		::NaNFreeEle(ele,&idx_ele, plane_idx_,nn_,locstep,*dof,is_2d_);
	} else {
		for (int i2 = 0; i2<nn_[2]-1; ++i2) {
			::NaNFreeEle(ele,&idx_ele, i2,nn_,locstep,*dof,is_2d_);
		}
	}
	// clean the dof by counting to how many elements they belong
	::CleanDoF(dof,n_dof,&n_dop, *ele,locstep,is_2d_);
	return EPTlibError::Success;
}
// Solve the global minimisation problem.
EPTlibError_t EPTGradient::
GlobalMinimisation(std::vector<std::complex<double> > *epsc,
	const std::vector<std::complex<double> > &g_plus,
	const std::vector<std::complex<double> > &g_z) {
	int ele_nodes = is_2d_?quad_nodes:vox_nodes;
	int n_dim = is_2d_?2:3;
	std::array<int,NDIM> step{1,nn_[0],nn_[0]*nn_[1]};
	int n_out = is_2d_?step[2]:n_vox_;
	// define the computational elements
	std::vector<int> locstep(ele_nodes,0);
	for (int locid = 0; locid<ele_nodes; ++locid) {
		for (int d = 0; d<n_dim; ++d) {
			locstep[locid] += locid&(1<<d) ? step[d] : 0;
		}
	}
	// select the degrees of freedom
	int n_ele = (nn_[0]-1)*(nn_[1]-1);
	if (!is_2d_) {
		n_ele *= (nn_[2]-1);
	}
	std::vector<int> dof(n_out);
	std::vector<int> ele(n_ele);
	int ndof;
	SelectDoF(&dof,&ele,&ndof, *epsc,locstep);
	// set the dirichlet values in case of seed points
	std::vector<std::complex<double> > dir(seed_points_.size());
	if (use_seed_points_) {
		::SeedPointsToDoF(&dir,&dof,&ndof, seed_points_,nn_,plane_idx_,omega_,is_2d_);
	}
	// create the matrix of coefficients and the forcing term
	int idx0 = is_2d_?plane_idx_*step[2]:0;
	Eigen::SparseMatrix<std::complex<double> > A(ndof,ndof);
	Eigen::SparseMatrix<std::complex<double> > B(ndof,ndof);
	Eigen::VectorXcd b;
	Eigen::VectorXcd l;
	::FillGlobalMatrix(&A,&B,&b,&l, dof,ndof,ele,locstep,*epsc,g_plus,g_z,dd_,
		idx0,is_2d_, use_seed_points_,dir);
	// solve the linear system
	Eigen::VectorXcd x(ndof);
	{
		int idx = idx0;
		for (int idx_out = 0; idx_out<n_out; ++idx_out) {
			if (dof[idx_out]>0) {
				x[dof[idx_out]-1] = std::log((*epsc)[idx]);
			}
			++idx;
		}
	}
	if (use_lcurve_) {
		Eigen::VectorXcd g(ndof);
		Eigen::VectorXcd x0(ndof);
		Eigen::VectorXcd z(ndof);
		for (int idx = 0; idx<n_out; ++idx) {
			if (dof[idx]>0) {
				g[dof[idx]-1] = g_plus[idx0+idx];
				if (!is_2d_) {
					g[dof[idx]-1] += g_z[idx0+idx];
				}
				x0[dof[idx]-1] = std::log((*epsc)[idx0+idx]);
			}
		}
		for (int idx = 0; idx<ndof; ++idx) {
			z[idx] = x[idx];
		}
		double lambda_k = lambda_;
		for (int k = 0; k<20; ++k) {
			::SolveGlobalSystem(&x, A+lambda_k*B,b+lambda_k*l);
			::SolveGlobalSystem(&z, A+lambda_k*B,A*x-b);

			std::complex<double> Fg = x.adjoint()*(A*x-2.0*b);
			Fg += g.adjoint()*(B*g);
			double rho = std::real(Fg);
			double eta = (x-x0).norm();
			Fg = -4.0/lambda_k*(x-x0).adjoint()*z;
			double eta_d = std::real(Fg);

			double kappa = 2.0*eta*rho/eta_d*(lambda_k*lambda_k*eta_d*rho + 2.0*lambda_k*eta*rho + lambda_k*lambda_k*lambda_k*lambda_k*eta*eta_d)/std::sqrt((lambda_k*lambda_k*eta*eta+rho*rho)*(lambda_k*lambda_k*eta*eta+rho*rho)*(lambda_k*lambda_k*eta*eta+rho*rho));
			std::cout<<lambda_k<<
				","<<rho<<
				","<<eta<<
				","<<Fg<<
				","<<kappa<<
				","<<std::endl;

			lambda_k *= 10.0;
		}
	} else {
		if (use_seed_points_) {
			::SolveGlobalSystem(&x, A+lambda_*B,b+lambda_*B*x);
		} else {
			::SolveGlobalSystem(&x, A,b);
		}
	}
	// extract electric properties
	{
		int idx = idx0;
		for (int idx_out = 0; idx_out<n_out; ++idx_out) {
			if (dof[idx_out]>0) {
				(*epsc)[idx] = std::exp(x[dof[idx_out]-1]);
			} else {
				if (use_seed_points_ && dof[idx_out]<0) {
					(*epsc)[idx] = std::exp(dir[-dof[idx_out]-1]);
				} else {
					(*epsc)[idx] = std::numeric_limits<double>::quiet_NaN();
				}
			}
			++idx;
		}
	}
	return EPTlibError::Success;
}

// Extract the electric properties of a slice.
EPTlibError_t EPTGradient::
ExtractElectricPropertiesSlice(const int i2, const int out_i2,
	const std::vector<std::complex<double> > &epsc) {
	int nxy = nn_[0]*nn_[1];
	int idx = nn_[0]*nn_[1]*i2;
	int out_idx = nn_[0]*nn_[1]*out_i2;
	for (int s = 0; s<nxy; ++s) {
		epsr_[out_idx] = epsc[idx].real()/EPS0;
		sigma_[out_idx] = -epsc[idx].imag()*omega_;
		++idx;
		++out_idx;
	}
	return EPTlibError::Success;
}

EPTlibError_t EPTGradient::
SelectPlane(const int plane_idx) {
    plane_idx_ = plane_idx;
    return EPTlibError::Success;
}
EPTlibError_t EPTGradient::
Set2D() {
	is_2d_ = true;
	return EPTlibError::Success;
}
EPTlibError_t EPTGradient::
Unset2D() {
	is_2d_ = false;
	return EPTlibError::Success;
}
EPTlibError_t EPTGradient::
SetLCurve() {
	if (use_seed_points_) {
		return EPTlibError::Unknown;
	}
	use_lcurve_ = true;
	return EPTlibError::Success;
}
EPTlibError_t EPTGradient::
UnsetLCurve() {
	use_lcurve_ = false;
	return EPTlibError::Success;
}
EPTlibError_t EPTGradient::
ToggleAverageWithMedian() {
	average_with_median_ = !average_with_median_;
	return EPTlibError::Success;
};
EPTlibError_t EPTGradient::
SetLambda(const double lambda) {
	if (lambda<0.0) {
		return EPTlibError::Unknown;
	}
	if (use_seed_points_) {
		return EPTlibError::Unknown;
	}
	lambda_ = lambda;
	return EPTlibError::Success;
}
EPTlibError_t EPTGradient::
AddSeedPoint(const SeedPoint seed_point) {
	use_seed_points_ = true;
	seed_points_.push_back(seed_point);
	use_lcurve_ = false;
	lambda_ = 0.0;
	return EPTlibError::Success;
}
