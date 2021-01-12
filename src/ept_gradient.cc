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

namespace {
	static double nand = std::numeric_limits<double>::quiet_NaN();
	static std::complex<double> nancd = std::complex<double>(nand,nand);
}

// EPTGradient constructor
EPTGradient::
EPTGradient(const double freq, const std::array<int,NDIM> &nn,
	const std::array<double,NDIM> &dd, const int tx_ch,
	const Shape &shape, const bool is_2d) :
	EPTInterface(freq,nn,dd,tx_ch, 1),
	is_2d_(is_2d), use_seed_points_(false),
	plane_idx_(nn[2]/2), lambda_(0.0), seed_points_(0),
	fd_filter_(shape),
	epsc_(0), g_plus_(0), g_z_(0), thereis_epsc_(false),
	run_mode_(EPTGradientRun::FULL) {
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
EPTlibError EPTGradient::
Run() {
	// Check the run mode conditions
	if (run_mode_==EPTGradientRun::GRADIENT && !thereis_epsc_) {
		return EPTlibError::MissingData;
	}
	// Set-up the output
	int n_out = is_2d_?nn_[0]*nn_[1]:n_vox_;
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
	// Solution of the local systems
	if (run_mode_==EPTGradientRun::FULL || run_mode_==EPTGradientRun::LOCAL) {
		thereis_epsc_ = true;
		// Clear and re-allocate the intermediate results
		epsc_.clear();
		epsc_.resize(n_out,0.0);
		g_plus_.clear();
		g_plus_.resize(n_out,0.0);
		g_z_.clear();
		if (!is_2d_) {
			g_z_.resize(n_out,0.0);
		}
		// For each transmit channel as reference...
		int r2 = fd_filter_.GetShape().GetSize()[2];
		int n_int = is_2d_?nn_[0]*nn_[1]*r2:n_vox_;
		std::vector<double> totweight(n_out,0.0);
		for (int iref = 0; iref<tx_ch_; ++iref) {
			// ...allocate memory for auxiliary quantities
			std::array<std::vector<double>,NDIM> grad_phi0;
			for (int d = 0; d<NDIM; ++d) {
				grad_phi0[d].resize(n_int,::nand);
			}
			std::vector<std::complex<double> > g_plus(n_int,::nancd);
			std::vector<std::complex<double> > g_z(n_int,::nancd);
			std::vector<std::complex<double> > theta(n_int,::nancd);
			// ...perform the local recovery
			LocalRecovery(&grad_phi0,&g_plus,&g_z,&theta,iref);
			// ...estimate the complex permittivity from theta
			EPTlibError error = Theta2Epsc(&theta,grad_phi0,g_plus,g_z);
			if (error!=EPTlibError::Success) {
				return error;
			}
			// ...contribute in averaging
			for (int idx = 0; idx<n_out; ++idx) {
				int idx_tx = is_2d_?(idx+nn_[0]*nn_[1]*plane_idx_):idx;
				int idx_int = is_2d_?(idx+nn_[0]*nn_[1]*(r2/2)):idx;
				double weight = (*tx_sens_[iref])[idx_tx];
				epsc_[idx] += theta[idx_int]*weight;
				g_plus_[idx] += g_plus[idx_int]*weight;
				if (!is_2d_) {
					g_z_[idx] += g_z[idx_int]*weight;
				}
				totweight[idx] += weight;
			}
		}
		// Average the epsc, g_plus and g_z quantities
		for (int idx = 0; idx<n_out; ++idx) {
			epsc_[idx] /= totweight[idx];
			g_plus_[idx] /= totweight[idx];
			if (!is_2d_) {
				g_z_[idx] /= totweight[idx];
			}
		}
	}
	// Optimise the estimate using the gradient information
	if (run_mode_==EPTGradientRun::FULL || run_mode_==EPTGradientRun::GRADIENT) {
		GlobalMinimisation();
	}
	// Extract the electric properties
	ExtractElectricProperties();
	return EPTlibError::Success;
}

void EPTGradient::
SetRunMode(EPTGradientRun run_mode) {
	run_mode_ = run_mode;
	return;
}
EPTGradientRun EPTGradient::
GetRunMode() {
	return run_mode_;
}
bool EPTGradient::
ToggleSeedPoints() {
	use_seed_points_ = !use_seed_points_;
	return use_seed_points_;
}
void EPTGradient::
AddSeedPoint(const SeedPoint seed_point) {
	seed_points_.push_back(seed_point);
	return;
}
EPTlibError EPTGradient::
SelectSlice(const int slice_idx) {
	int r2 = fd_filter_.GetShape().GetSize()[2]/2;
	if (slice_idx<2*r2||slice_idx>nn_[2]-1-2*r2) {
		return EPTlibError::WrongDataFormat;
	}
    plane_idx_ = slice_idx;
    return EPTlibError::Success;
}
EPTlibError EPTGradient::
SetLambda(const double lambda) {
	if (lambda<0.0) {
		return EPTlibError::WrongDataFormat;
	}
	lambda_ = lambda;
	return EPTlibError::Success;
}

EPTlibError EPTGradient::
GetEpsC(std::vector<std::complex<double> > *epsc) {
	if (!thereis_epsc_) {
		return EPTlibError::MissingData;
	}
	*epsc = epsc_;
	return EPTlibError::Success;
}
EPTlibError EPTGradient::
GetGPlus(std::vector<std::complex<double> > *g_plus) {
	if (!thereis_epsc_) {
		return EPTlibError::MissingData;
	}
	*g_plus = g_plus_;
	return EPTlibError::Success;
}
EPTlibError EPTGradient::
GetGZ(std::vector<std::complex<double> > *g_z) {
	if (!thereis_epsc_ || is_2d_) {
		return EPTlibError::MissingData;
	}
	*g_z = g_z_;
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
	void FillLocalMatrix(Eigen::MatrixXd *A, Eigen::VectorXd *b,
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
		return;
	}
	// Solve the local system for pixel-by-pixel recovery.
	void SolveLocalSystem(double *gradx_phi0, double *grady_phi0,
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
		return;
	}
}  // namespace
// Perform the pixel-by-pixel recovery
void EPTGradient::
LocalRecovery(std::array<std::vector<double>,NDIM> *grad_phi0,
	std::vector<std::complex<double> > *g_plus,
	std::vector<std::complex<double> > *g_z,
	std::vector<std::complex<double> > *theta,
	const int iref) {
	int r2 = fd_filter_.GetShape().GetSize()[2]/2;
	if (is_2d_) {
		for (int i2 = plane_idx_-r2; i2<=plane_idx_+r2; ++i2) {
			LocalRecoverySlice(grad_phi0,g_plus,g_z,theta,iref,i2);
		}
	} else {
		for (int i2 = r2; i2<nn_[2]-r2; ++i2) {
			LocalRecoverySlice(grad_phi0,g_plus,g_z,theta,iref,i2);
		}
	}
	return;
}
// Perform pixel-by-pixel recovery in a slice.
void EPTGradient::
LocalRecoverySlice(std::array<std::vector<double>,NDIM> *grad_phi0,
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
			int idx = ii[0] + nn_[0]*ii[1];
			if (is_2d_) {
				idx += nn_[0]*nn_[1]*(ii[2]-plane_idx_+rr[2]);
			} else {
				idx += nn_[0]*nn_[1]*ii[2];
			}
			::SolveLocalSystem(&((*grad_phi0)[0][idx]),&((*grad_phi0)[1][idx]),
				&((*grad_phi0)[2][idx]),&((*g_plus)[idx]),&((*g_z)[idx]),
				&((*theta)[idx]), A,b,is_2d_);
		}
	}
	return;
}

/////////////////////////////////////////////////////////////////////////////
// EPTGradient THETA TO COMPLEX PERMITTIVITY ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
// Estimate the complex permittivity from theta local recovery.
EPTlibError EPTGradient::
Theta2Epsc(std::vector<std::complex<double> > *theta,
	const std::array<std::vector<double>,NDIM> &grad_phi0,
	const std::vector<std::complex<double> > &g_plus,
	const std::vector<std::complex<double> > &g_z) {
	int n_dim = is_2d_?2:3;
	int r2 = fd_filter_.GetShape().GetSize()[2];
	std::array<int,NDIM> nn{nn_[0],nn_[1],is_2d_?r2:nn_[2]};
	int n_vox = nn[0]*nn[1]*nn[2];
	int offset = is_2d_?nn[0]*nn[1]*(r2/2):0;
	// grad_phi0 laplacian...
	for (int d = 0; d<n_dim; ++d) {
		std::vector<double> lapl_phi0(n_vox);
		DifferentialOperator diff_op = static_cast<DifferentialOperator>(d);
		EPTlibError error = fd_filter_.Apply(diff_op,lapl_phi0.data(),
			grad_phi0[d].data(),nn,dd_);
		if (error!=EPTlibError::Success) {
			return error;
		}
		std::transform(theta->begin()+offset,
			theta->end()-offset,
			lapl_phi0.begin()+offset,
			theta->begin()+offset,
			[](const std::complex<double> &a,const double &b)->std::complex<double>{
				return a + std::complex<double>(0.0,b);
			});
	}
	// grad_phi0 magnitude...
	for (int d = 0; d<n_dim; ++d) {
		std::vector<double> tmp(n_vox);
		std::transform(grad_phi0[d].begin()+offset,
			grad_phi0[d].end()-offset,
			tmp.begin()+offset,
			[](const double &a)->double{
				return a*a;
			});
		std::transform(theta->begin()+offset,
			theta->end()-offset,
			tmp.begin()+offset,
			theta->begin()+offset,
			[](const std::complex<double> &a,const double &b)->std::complex<double>{
				return a-b;
			});
	}
	// grad_phi0, g inner product...
	std::vector<std::complex<double> > tmp(n_vox);
	// ...component 0
	std::transform(grad_phi0[0].begin()+offset,
		grad_phi0[0].end()-offset,
		g_plus.begin()+offset,
		tmp.begin()+offset,
		[](const double &a,const std::complex<double> &b)->std::complex<double>{
			return std::complex<double>(0.0,a)*b;
		});
	std::transform(theta->begin()+offset,
		theta->end()-offset,
		tmp.begin()+offset,
		theta->begin()+offset,
		std::minus<std::complex<double> >());
	// ...component 1
	std::transform(grad_phi0[1].begin()+offset,
		grad_phi0[1].end()-offset,
		g_plus.begin()+offset,
		tmp.begin()+offset,
		[](const double &a,const std::complex<double> &b)->std::complex<double>{
			return a*b;
		});
	std::transform(theta->begin()+offset,
		theta->end()-offset,
		tmp.begin()+offset,
		theta->begin()+offset,
		std::minus<std::complex<double> >());
	if (!is_2d_) {
		// ...component 2
		std::transform(grad_phi0[2].begin()+offset,
			grad_phi0[2].end()-offset,
			g_z.begin()+offset,
			tmp.begin()+offset,
			[](const double &a,const std::complex<double> &b)->std::complex<double>{
				return std::complex<double>(0.0,a)*b;
			});
		std::transform(theta->begin()+offset,
			theta->end()-offset,
			tmp.begin()+offset,
			theta->begin()+offset,
			std::minus<std::complex<double> >());
	}
	// multiplicative factor...
	std::transform(theta->begin()+offset,
		theta->end()-offset,
		theta->begin()+offset,
		[this](const std::complex<double> &a)->std::complex<double>{
			return -a/omega_/omega_/MU0;
		});
	return EPTlibError::Success;
}

/////////////////////////////////////////////////////////////////////////////
// EPTGradient GLOBAL OPTIMISATION //////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
namespace { //details
	static constexpr int vox_nodes = 8;
	static constexpr int quad_nodes = 4;

	// Select the NaN-free elements.
	void NaNFreeEle(std::vector<int> *ele, int *idx_ele, const int i2,
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
		return;
	}
	// Clean dof by counting to how many elements they belong.
	void CleanDoF(std::vector<int> *dof, int *n_dof, int *n_dop,
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
		return;
	}

	// Associate each seed point to its dof.
	void SeedPointsToDoF(std::vector<std::complex<double> > *dir,
		std::vector<int> *dof, int *n_dof,
		const std::vector<SeedPoint> &seed_points,
		const std::array<int,NDIM> &nn, const int i2, const double omega,
		const bool is_2d) {
		// for each seed find the corresponding idx
		auto n_out = dof->size();
		std::vector<int> idx2sp(n_out,-1);
		for (int idx_sp = 0; idx_sp<seed_points.size(); ++idx_sp) {
			SeedPoint sp = seed_points[idx_sp];
			int idx = sp.ijk[0] + nn[0]*sp.ijk[1];
			if (!is_2d) {
				idx += nn[0]*nn[1]*sp.ijk[2];
			}
			idx2sp[idx] = idx_sp;
			(*dir)[idx_sp] = std::log(std::complex<double>(sp.epsr*EPS0,
				-sp.sigma/omega));
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
		return;
	}
	
	// Fill the global matrix for minimisation.
	void FillGlobalMatrix(Eigen::SparseMatrix<std::complex<double> > *A,
		Eigen::SparseMatrix<std::complex<double> > *B, Eigen::VectorXcd *b,
		Eigen::VectorXcd *l, const std::vector<int> &dof, const int n_dof,
		const std::vector<int> &ele, const std::vector<int> &locstep,
		const std::vector<std::complex<double> > &epsc,
		const std::vector<std::complex<double> > &g_plus,
		const std::vector<std::complex<double> > &g_z,
		const std::array<double,NDIM> &dd, const bool is_2d,
		const bool use_seed_points_,
		const std::vector<std::complex<double> > &dir) {
		const int ele_nodes = is_2d?quad_nodes:vox_nodes;
		const int n_dim = is_2d?2:3;
		std::vector<Eigen::Triplet<std::complex<double> > > A_trip(0);
		std::vector<Eigen::Triplet<std::complex<double> > > B_trip(0);
		*b = Eigen::VectorXcd::Zero(n_dof);
		*l = Eigen::VectorXcd::Zero(n_dof);
		// loop the elements
		for (int idx_ele = 0; idx_ele<ele.size(); ++idx_ele) {
			int idx = ele[idx_ele];
			if (idx<0) {
				continue;
			}
			// loop the equation index
			for (int i = 0; i<ele_nodes; ++i) {
				int iidx = idx+locstep[i];
				int ii = dof[iidx];
				if (ii>0) {
					--ii;
					// loop the unknown index
					for (int j = 0; j<ele_nodes; ++j) {
						int jidx = idx+locstep[j];
						int jj = dof[jidx];
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
		return;
	}
	// Solve the global system for minimisation.
	void SolveGlobalSystem(Eigen::VectorXcd *x,
		const Eigen::SparseMatrix<std::complex<double> > &A,
		const Eigen::VectorXcd &b) {
		Eigen::ConjugateGradient<Eigen::SparseMatrix<std::complex<double> > > solver;
		solver.compute(A);
		*x = solver.solve(b);
		return;
	}

	// Implements the l-curve method. !NOT IMPLEMENTED YET!
	void LCurve(Eigen::VectorXcd *x,
		const Eigen::SparseMatrix<std::complex<double> > &A,
		const Eigen::SparseMatrix<std::complex<double> > &B,
		const Eigen::VectorXcd &b,
		const Eigen::VectorXcd &l,
		const double lambda, const std::vector<int> &dof, const int n_dof,
		const std::vector<std::complex<double> > &epsc,
		const std::vector<std::complex<double> > &g_plus,
		const std::vector<std::complex<double> > &g_z,
		const bool is_2d) {
		const auto n_vox = dof.size();
		Eigen::VectorXcd g(n_dof);
		Eigen::VectorXcd x0(n_dof);
//		Eigen::VectorXcd z(n_dof);
		for (int idx = 0; idx<n_vox; ++idx) {
			if (dof[idx]>0) {
				g[dof[idx]-1] = g_plus[idx];
				if (!is_2d) {
					g[dof[idx]-1] += g_z[idx];
				}
				x0[dof[idx]-1] = std::log(epsc[idx]);
			}
		}
//		for (int idx = 0; idx<n_dof; ++idx) {
//			z[idx] = (*x)[idx];
//		}
		double lambda_k = lambda;
		for (int k = 0; k<20; ++k) {
			::SolveGlobalSystem(x, A+lambda_k*B,b+lambda_k*l);
//			::SolveGlobalSystem(&z, A+lambda_k*B,A*(*x)-b);
			std::complex<double> Fg = x->adjoint()*(A*(*x)-2.0*b);
			std::complex<double> tmp = g.adjoint()*(B*g);
			Fg += tmp;
			double rho = std::real(Fg);
			double eta = (*x-x0).norm();
//			Fg = -4.0/lambda_k*(*x-x0).adjoint()*z;
//			double eta_d = std::real(Fg);
//			double kappa = 2.0*eta*rho/eta_d*(lambda_k*lambda_k*eta_d*rho + 2.0*lambda_k*eta*rho + lambda_k*lambda_k*lambda_k*lambda_k*eta*eta_d)/std::sqrt((lambda_k*lambda_k*eta*eta+rho*rho)*(lambda_k*lambda_k*eta*eta+rho*rho)*(lambda_k*lambda_k*eta*eta+rho*rho));
			lambda_k *= 10.0;
		}
		return;
	}
}  // namespace
// Select the degree of freedom.
void EPTGradient::
SelectDoF(std::vector<int> *dof, std::vector<int> *ele, int *n_dof,
	const std::vector<int> &locstep) {
	// select the dof
	int n_dop = 0;
	*n_dof = 0;
	for (int idx = 0; idx<dof->size(); ++idx) {
		if (epsc_[idx]==epsc_[idx]) {
			(*dof)[idx] = ++(*n_dof);
		} else {
			(*dof)[idx] = --n_dop;
		}
	}
	// determine the nan-free elements
	int idx_ele = 0;
	if (is_2d_) {
		::NaNFreeEle(ele,&idx_ele, 0,nn_,locstep,*dof,is_2d_);
	} else {
		for (int i2 = 0; i2<nn_[2]-1; ++i2) {
			::NaNFreeEle(ele,&idx_ele, i2,nn_,locstep,*dof,is_2d_);
		}
	}
	// clean the dof by counting to how many elements they belong
	::CleanDoF(dof,n_dof,&n_dop, *ele,locstep,is_2d_);
	return;
}
// Solve the global minimisation problem.
void EPTGradient::
GlobalMinimisation() {
	const int ele_nodes = is_2d_?quad_nodes:vox_nodes;
	const int n_dim = is_2d_?2:3;
	const std::array<int,NDIM> step{1,nn_[0],nn_[0]*nn_[1]};
	const int n_out = is_2d_?step[2]:n_vox_;
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
	SelectDoF(&dof,&ele,&ndof,locstep);
	// set the dirichlet values in case of seed points
	std::vector<std::complex<double> > dir(seed_points_.size());
	if (use_seed_points_) {
		::SeedPointsToDoF(&dir,&dof,&ndof, seed_points_,nn_,plane_idx_,omega_,is_2d_);
	}
	// create the matrix of coefficients and the forcing term
	Eigen::SparseMatrix<std::complex<double> > A(ndof,ndof);
	Eigen::SparseMatrix<std::complex<double> > B(ndof,ndof);
	Eigen::VectorXcd b;
	Eigen::VectorXcd l;
	::FillGlobalMatrix(&A,&B,&b,&l, dof,ndof,ele,locstep,epsc_,g_plus_,g_z_,
		dd_,is_2d_, use_seed_points_,dir);
	// solve the linear system
	Eigen::VectorXcd x(ndof);
	if (!use_seed_points_) {
		/* !NOT IMPLEMENTED YET! */
		for (int idx = 0; idx<n_out; ++idx) {
			if (dof[idx]>0) {
				x[dof[idx]-1] = std::log(epsc_[idx]);
			}
			++idx;
		}
		::LCurve(&x, A,B,b,l,lambda_,dof,ndof,epsc_,g_plus_,g_z_,is_2d_);
	} else {
		::SolveGlobalSystem(&x, A,b);
	}
	// extract electric properties
	for (int idx = 0; idx<n_out; ++idx) {
		if (dof[idx]>0) {
			epsc_[idx] = std::exp(x[dof[idx]-1]);
		} else {
			if (use_seed_points_ && dof[idx]<0) {
				epsc_[idx] = std::exp(dir[-dof[idx]-1]);
			} else {
				epsc_[idx] = ::nancd;
			}
		}
	}
	return;
}

// Extract the electric properties of a slice.
void EPTGradient::
ExtractElectricProperties() {
	for (int idx = 0; idx<epsc_.size(); ++idx) {
		epsr_[idx] = epsc_[idx].real()/EPS0;
		sigma_[idx] = -epsc_[idx].imag()*omega_;
	}
	return;
}
