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

#include "eptlib/ept_gradient.h"

#include <cmath>
#include <limits>

#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace eptlib;
using namespace std::literals::complex_literals;

namespace {
    static double nand = std::numeric_limits<double>::quiet_NaN();
    static std::complex<double> nancd = std::complex<double>(nand,nand);
}

// EPTGradient constructor
EPTGradient::
EPTGradient(const size_t n0, const size_t n1, const size_t n2,
    const double d0, const double d1, const double d2,
    const double freq, const Shape &window, const size_t n_tx_ch,
    const int degree, const EPTGradientRun run_mode) :
    EPTInterface(n0,n1,n2, d0,d1,d2, freq, n_tx_ch,1, false),
    slice_index_(),
    seed_points_(),
    lambda_(0.0),
    gradient_tolerance_(0.0),
    mask_(0),
    sg_filter_(d0,d1,d2, window, degree),
    run_mode_(run_mode),
    epsc_(0),
    g_plus_(0),
    g_z_(0),
    cost_functional_(0.0),
    cost_regularization_(0.0),
    thereis_epsc_(false) {
    if (n_tx_ch_<5) {
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
    if (run_mode_==EPTGradientRun::GRADIENT && !thereis_epsc_) {
        return EPTlibError::MissingData;
    }
    if (!ThereAreAllTxSens() || !ThereAreAllTRxPhase()) {
        return EPTlibError::MissingData;
    }
    // Set-up the output
    epsr_  = std::make_unique<Image<double> >(nn_[0], nn_[1], VolumeTomography() ? nn_[2] : 1);
    sigma_ = std::make_unique<Image<double> >(nn_[0], nn_[1], VolumeTomography() ? nn_[2] : 1);
    size_t n_vox = GetTRxPhase(0,0)->GetNVox();
    size_t n_out = sigma_->GetNVox();
    // Solution of the local systems
    if (run_mode_==EPTGradientRun::FULL || run_mode_==EPTGradientRun::LOCAL) {
        thereis_epsc_ = true;
        // Clear and re-allocate the intermediate results
        epsc_  .assign(n_out, 0.0);
        g_plus_.assign(n_out, 0.0);
        if (VolumeTomography()) {
            g_z_.assign(n_out, 0.0);
        }
        // For each transmit channel as reference...
        size_t m2 = sg_filter_.GetWindow().GetSize(2);
        size_t r2 = m2/2;
        std::array<size_t, N_DIM> nn{nn_[0], nn_[1], VolumeTomography() ? nn_[2] : m2};
        size_t n_int = nn[0]*nn[1]*nn[2];
        std::vector<double> totweight(n_out, 0.0);
        for (size_t tx_ch_ref = 0; tx_ch_ref<n_tx_ch_; ++tx_ch_ref) {
            // ...allocate memory for auxiliary quantities
            std::array<Image<double>, N_DIM> grad_phi0;
            for (size_t d = 0; d<N_DIM; ++d) {
                grad_phi0[d] = Image<double>(nn[0], nn[1], nn[2]);
                grad_phi0[d].GetData().assign(n_int, ::nand);
            }
            std::vector<std::complex<double> > g_plus(n_int, ::nancd);
            std::vector<std::complex<double> > g_z   (n_int, ::nancd);
            std::vector<std::complex<double> > theta (n_int, ::nancd);
            // ...perform the local recovery
            LocalRecovery(&grad_phi0,&g_plus,&g_z,&theta, tx_ch_ref);
            // ...estimate the complex permittivity from theta
            EPTlibError error = Theta2Epsc(&theta, grad_phi0,g_plus,g_z);
            if (error!=EPTlibError::Success) {
                return error;
            }
            // ...contribute in averaging
            for (size_t idx = 0; idx<n_out; ++idx) {
                size_t idx_tx  = VolumeTomography() ? idx : idx+n_out*slice_index_.value();
                size_t idx_int = VolumeTomography() ? idx : idx+n_out*r2;
                double weight = GetTxSens(tx_ch_ref)->At(idx_tx);
                epsc_[idx] += theta[idx_int]*weight;
                g_plus_[idx] += g_plus[idx_int]*weight;
                if (VolumeTomography()) {
                    g_z_[idx] += g_z[idx_int]*weight;
                }
                totweight[idx] += weight;
            }
        }
        // Average the epsc, g_plus and g_z quantities
        for (int idx = 0; idx<n_out; ++idx) {
            epsc_[idx] /= totweight[idx];
            g_plus_[idx] /= totweight[idx];
            if (VolumeTomography()) {
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
        const std::vector<std::complex<double> > &field_mid_values,
        const filter::SavitzkyGolay &sg_filter,
        const std::array<double,N_DIM> &dd, const int tx_ch, const bool is_2d) {
        int n_col = is_2d?6:9;
        *A = Eigen::MatrixXd::Zero(2*tx_ch,n_col);
        *b = Eigen::VectorXd::Zero(2*tx_ch);
        // for each transmit channel...
        for (int itx = 0; itx<tx_ch; ++itx) {
            int row = 2*itx;
            // ...laplacian
            std::complex<double> tmp = sg_filter.Laplacian(field_crop[itx]);
            (*b)[row  ] = tmp.real();
            (*b)[row+1] = tmp.imag();
            // ...grad_x
            tmp = sg_filter.FirstOrderDerivative(0, field_crop[itx]);
            (*A)(row  , col_phi    ) =  tmp.imag()*2.0;
            (*A)(row+1, col_phi    ) = -tmp.real()*2.0;
            (*A)(row  , col_gplus  ) =  tmp.real();
            (*A)(row+1, col_gplus  ) =  tmp.imag();
            (*A)(row  , col_gplus+1) = -tmp.imag();
            (*A)(row+1, col_gplus+1) =  tmp.real();
            // ...grad_y
            tmp = sg_filter.FirstOrderDerivative(1, field_crop[itx]);
            (*A)(row  , col_phi+1  )  =  tmp.imag()*2.0;
            (*A)(row+1, col_phi+1  )  = -tmp.real()*2.0;
            (*A)(row  , col_gplus  ) +=  tmp.imag();
            (*A)(row+1, col_gplus  ) += -tmp.real();
            (*A)(row  , col_gplus+1) +=  tmp.real();
            (*A)(row+1, col_gplus+1) +=  tmp.imag();
            // ...grad_z
            if (!is_2d) {
                tmp = sg_filter.FirstOrderDerivative(2, field_crop[itx]);
                (*A)(row  , col_phi+2) =  tmp.imag()*2.0;
                (*A)(row+1, col_phi+2) = -tmp.real()*2.0;
                (*A)(row  , col_gz   ) =  tmp.real();
                (*A)(row+1, col_gz   ) =  tmp.imag();
                (*A)(row  , col_gz+1 ) = -tmp.imag();
                (*A)(row+1, col_gz+1 ) =  tmp.real();
            }
            // ...value
            tmp = field_mid_values[itx];
            (*A)(row  , col_theta  ) =  tmp.real();
            (*A)(row+1, col_theta  ) =  tmp.imag();
            (*A)(row  , col_theta+1) = -tmp.imag();
            (*A)(row+1, col_theta+1) =  tmp.real();
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
        *theta  = std::complex<double>(x[col_theta], x[col_theta+1]);
        *g_plus = std::complex<double>(x[col_gplus], x[col_gplus+1]);
        *gradx_phi0 = x[col_phi];
        *grady_phi0 = x[col_phi+1];
        if (!is_2d) {
            *gradz_phi0 = x[col_phi+2];
            *g_z = std::complex<double>(x[col_gz], x[col_gz+1]);
        }
        return;
    }

}  // namespace

// Perform the pixel-by-pixel recovery
void EPTGradient::
LocalRecovery(std::array<Image<double>, N_DIM> *grad_phi0,
    std::vector<std::complex<double> > *g_plus,
    std::vector<std::complex<double> > *g_z,
    std::vector<std::complex<double> > *theta,
    const int iref) {
    size_t r2 = sg_filter_.GetWindow().GetSize(2)/2;
    if (!VolumeTomography()) {
        for (size_t i2 = slice_index_.value()-r2; i2<=slice_index_.value()+r2; ++i2) {
            LocalRecoverySlice(grad_phi0,g_plus,g_z,theta,iref,i2);
        }
    } else {
        for (size_t i2 = r2; i2<nn_[2]-r2; ++i2) {
            LocalRecoverySlice(grad_phi0,g_plus,g_z,theta,iref,i2);
        }
    }
    return;
}

// Perform pixel-by-pixel recovery in a slice.
void EPTGradient::
LocalRecoverySlice(std::array<Image<double>, N_DIM> *grad_phi0,
    std::vector<std::complex<double> > *g_plus,
    std::vector<std::complex<double> > *g_z,
    std::vector<std::complex<double> > *theta,
    const int iref, const int i2) {
    const Shape& window = sg_filter_.GetWindow();
    std::array<int, N_DIM> rr;
    for (size_t d = 0; d<N_DIM; ++d) {
        rr[d] = window.GetSize(d)/2;
    }
    std::array<size_t, N_DIM> ii;
    std::array<size_t, N_DIM> m_inc;
    size_t m_vox = window.GetNVox();
    m_inc[0] = 1;
    m_inc[1] = nn_[0] - window.GetSize(0);
    m_inc[2] = nn_[0] * (nn_[1] - window.GetSize(1));
    ii[2] = i2;
    for (ii[1] = rr[1]; ii[1]<nn_[1]-rr[1]; ++ii[1]) {
        for (ii[0] = rr[0]; ii[0]<nn_[0]-rr[0]; ++ii[0]) {
            std::array<int, N_DIM> ii_l;
            std::vector<std::vector<std::complex<double> > > field_crop(n_tx_ch_);
            std::vector<std::complex<double> > field_mid_values(n_tx_ch_, 0.0);
            for (int tx_ch = 0; tx_ch<n_tx_ch_; ++tx_ch) {
                field_crop[tx_ch].assign(window.GetVolume(), 0.0);
            }
            // inner loop over the kernel voxels
            int idx_f = 0;
            int idx_l = 0;
            int idx_g = ii[0]-rr[0] + nn_[0]*(ii[1]-rr[1] + nn_[1]*(ii[2]-rr[2]));
            for (ii_l[2] = -rr[2]; ii_l[2]<=rr[2]; ++ii_l[2]) {
                for (ii_l[1] = -rr[1]; ii_l[1]<=rr[1]; ++ii_l[1]) {
                    for (ii_l[0] = -rr[0]; ii_l[0]<=rr[0]; ++ii_l[0]) {
                        if (ii_l[0]==0 && ii_l[1]==0 && ii_l[2]==0) {
                            for (int tx_ch = 0; tx_ch<n_tx_ch_; ++tx_ch) {
                                std::complex<double> exponent = std::complex<double>(0.0, GetTRxPhase(tx_ch,0)->At(idx_g) - GetTRxPhase(iref,0)->At(idx_g));
                                field_mid_values[tx_ch] = GetTxSens(tx_ch)->At(idx_g) * std::exp(exponent);
                            }
                        }
                        if (window(idx_l)) {
                            // set the field value
                            for (int tx_ch = 0; tx_ch<n_tx_ch_; ++tx_ch) {
                                std::complex<double> exponent = std::complex<double>(0.0, GetTRxPhase(tx_ch,0)->At(idx_g) - GetTRxPhase(iref,0)->At(idx_g));
                                field_crop[tx_ch][idx_f] = GetTxSens(tx_ch)->At(idx_g) * std::exp(exponent);
                            }
                            ++idx_f;
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
            ::FillLocalMatrix(&A,&b, field_crop,field_mid_values, sg_filter_,dd_,n_tx_ch_,!VolumeTomography());
            // solve the system
            int idx = ii[0] + nn_[0]*ii[1];
            if (!VolumeTomography()) {
                idx += nn_[0]*nn_[1]*(ii[2]-slice_index_.value()+rr[2]);
            } else {
                idx += nn_[0]*nn_[1]*ii[2];
            }
            ::SolveLocalSystem(&((*grad_phi0)[0](idx)), &((*grad_phi0)[1](idx)),
                &((*grad_phi0)[2](idx)), &((*g_plus)[idx]), &((*g_z)[idx]),
                &((*theta)[idx]), A, b, !VolumeTomography());
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
    const std::array<Image<double>, N_DIM> &grad_phi0,
    const std::vector<std::complex<double> > &g_plus,
    const std::vector<std::complex<double> > &g_z) {
    int n_dim = VolumeTomography() ? 3 : 2;
    size_t m2 = sg_filter_.GetWindow().GetSize(2);
    std::array<size_t, N_DIM> nn{nn_[0], nn_[1], VolumeTomography() ? nn_[2] : m2};
    int n_vox = nn[0]*nn[1]*nn[2];
    int offset = VolumeTomography() ? 0 : nn[0]*nn[1]*(m2/2);
    // grad_phi0 laplacian...
    for (int d = 0; d<n_dim; ++d) {
        Image<double> lapl_phi0(nn[0], nn[1], nn[2]);
        filter::DifferentialOperator diff_op = static_cast<filter::DifferentialOperator>(d+1);
        EPTlibError error = sg_filter_.Apply(diff_op, &lapl_phi0, grad_phi0[d]);
        if (error!=EPTlibError::Success) {
            return error;
        }
        std::transform(theta->begin()+offset,
            theta->end()-offset,
            lapl_phi0.GetData().begin()+offset,
            theta->begin()+offset,
            [](const std::complex<double> &a,const double &b)->std::complex<double>{
                return a + 1.0i*b;
            });
    }
    // grad_phi0 magnitude...
    for (int d = 0; d<n_dim; ++d) {
        std::vector<double> tmp(n_vox);
        std::transform(grad_phi0[d].GetData().begin()+offset,
            grad_phi0[d].GetData().end()-offset,
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
    std::transform(grad_phi0[0].GetData().begin()+offset,
        grad_phi0[0].GetData().end()-offset,
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
    std::transform(grad_phi0[1].GetData().begin()+offset,
        grad_phi0[1].GetData().end()-offset,
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
    if (VolumeTomography()) {
        // ...component 2
        std::transform(grad_phi0[2].GetData().begin()+offset,
            grad_phi0[2].GetData().end()-offset,
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
        const std::array<int,N_DIM> &nn, const std::vector<int> &locstep,
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
        const std::array<int,N_DIM> &nn, const double omega,
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
        Eigen::SparseMatrix<std::complex<double> > *B, 
        Eigen::SparseMatrix<std::complex<double> > *Bout, Eigen::VectorXcd *b,
        Eigen::VectorXcd *b_dir, Eigen::VectorXcd *l, const std::vector<int> &dof,
        const int n_dof, const std::vector<int> &ele,
        const std::vector<int> &locstep,
        const std::vector<std::complex<double> > &epsc,
        const std::vector<std::complex<double> > &g_plus,
        const std::vector<std::complex<double> > &g_z,
        const std::array<double,N_DIM> &dd, const bool is_2d,
        const bool use_seed_points_,
        const std::vector<std::complex<double> > &dir,
        const boost::dynamic_bitset<> &mask) {
        const int ele_nodes = is_2d?quad_nodes:vox_nodes;
        const int n_dim = is_2d?2:3;
        std::vector<Eigen::Triplet<std::complex<double> > > A_trip(0);
        std::vector<Eigen::Triplet<std::complex<double> > > B_trip(0);
        std::vector<Eigen::Triplet<std::complex<double> > > Bout_trip(0);
        *b = Eigen::VectorXcd::Zero(n_dof);
        *b_dir = Eigen::VectorXcd::Zero(n_dof);
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
                        std::array<bool,N_DIM> same_side;
                        for (int d = 0; d<N_DIM; ++d) {
                            same_side[d] = (i&(1<<d))==(j&(1<<d));
                        }
                        // matrix B and vector l
                        std::complex<double> Bij = 1.0;
                        for (int d = 0; d<n_dim; ++d) {
                            Bij *= same_side[d] ? dd[d]/3.0 : dd[d]/6.0;
                        }
                        if (use_seed_points_ || !mask[idx_ele]) {
                            if (jj>0) {
                                Bout_trip.push_back(Eigen::Triplet<std::complex<double> >(ii,jj-1,Bij));
                            }
                        } else {
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
                            if (!is_2d) {
                                tmp *= same_side[2] ? dd[2]/3.0 : dd[2]/6.0;
                            }
                            Aij += std::complex<double>(0.0,tmp);
                        }
                        {
                            double tmp = -1.0;
                            tmp *= j&(1<<0) ? 0.5 : -0.5;
                            tmp *= i&(1<<1) ? 0.5 : -0.5;
                            if (!is_2d) {
                                tmp *= same_side[2] ? dd[2]/3.0 : dd[2]/6.0;
                            }
                            Aij += std::complex<double>(0.0,tmp);
                        }
                        if (jj>0) {
                            A_trip.push_back(Eigen::Triplet<std::complex<double> >(ii,jj-1,Aij));
                        } else {
                            (*b_dir)[ii] += -Aij*dir[-jj-1];
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
        (*Bout).setFromTriplets(Bout_trip.begin(),Bout_trip.end());
        (*Bout).makeCompressed();
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

    // Set a gradient threshold.
    double SetGradientThreshold(const double relative_threshold,
        const std::vector<std::complex<double> > &g_plus,
        const std::vector<std::complex<double> > &g_z, const bool is_2d) {
        double Gmax = 0;
        // loop the g_plus values
        for (int idx = 0; idx<g_plus.size(); ++idx) {
            if (g_plus[idx]==g_plus[idx]) {
                double G = std::norm(g_plus[idx]);
                if (!is_2d) {
                    G += std::norm(g_z[idx]);
                }
                if (G>Gmax) {
                    Gmax = G;
                }
            }
        }
        return relative_threshold*Gmax;
    }

    // Fill the gradient mask.
    void FillGradientMask(boost::dynamic_bitset<> *mask,
        const double threshold, const std::vector<int> &ele,
        const std::vector<int> &locstep,
        const std::vector<std::complex<double> > &g_plus,
        const std::vector<std::complex<double> > &g_z, const bool is_2d) {
        const int ele_nodes = is_2d?quad_nodes:vox_nodes;
        const int n_dim = is_2d?2:3;
        // loop the elements
        for (int idx_ele = 0; idx_ele<ele.size(); ++idx_ele) {
            int idx = ele[idx_ele];
            if (idx<0) {
                continue;
            }
            mask->set(idx_ele,true);
            // loop the equation index
            for (int i = 0; i<ele_nodes; ++i) {
                int iidx = idx+locstep[i];
                double G = std::norm(g_plus[iidx]);
                if (!is_2d) {
                    G += std::norm(g_z[iidx]);
                }
                if (G>threshold) {
                    mask->set(idx_ele,false);
                    break;
                }
            }
        }
        return;
    }

    // Evaluate the costs functional terms.
    void EvaluateCostFunctional(double *cost_functional,
        double *cost_regularization, const Eigen::VectorXcd &x,
        const Eigen::SparseMatrix<std::complex<double> > &A,
        const Eigen::SparseMatrix<std::complex<double> > &B,
        const Eigen::SparseMatrix<std::complex<double> > &Bout,
        const Eigen::VectorXcd &b, const std::vector<int> &dof,
        const int n_dof, const std::vector<std::complex<double> > &epsc,
        const std::vector<std::complex<double> > &g_plus,
        const std::vector<std::complex<double> > &g_z, const bool is_2d) {
        const auto n_vox = dof.size();
        Eigen::VectorXcd gp(n_dof);
        Eigen::VectorXcd gz(is_2d?0:n_dof);
        Eigen::VectorXcd x0(n_dof);
        for (int idx = 0; idx<n_vox; ++idx) {
            if (dof[idx]>0) {
                gp[dof[idx]-1] = g_plus[idx];
                if (!is_2d) {
                    gz[dof[idx]-1] = g_z[idx];
                }
                x0[dof[idx]-1] = std::log(epsc[idx]);
            }
        }
        double norm2_Dpu = std::real(std::complex<double>(x.adjoint()*(A*x)));
        double norm2_gp = std::real(std::complex<double>(gp.adjoint()*((B+Bout)*gp)));
        *cost_functional = norm2_Dpu + norm2_gp - 2.0*std::real(std::complex<double>(x.adjoint()*b));
        if (!is_2d) {
            *cost_functional += std::real(std::complex<double>(gz.adjoint()*((B+Bout)*gz)));
        }
        *cost_regularization = std::real(std::complex<double>(x.adjoint()*(B*(x-2.0*x0)))+std::complex<double>(x0.adjoint()*(B*x0)));
        return;
    }

}  // namespace

// Select the degree of freedom.
void EPTGradient::
SelectDoF(std::vector<int> *dof, std::vector<int> *ele, int *n_dof,
    const std::vector<int> &locstep) {
    const std::array<int,N_DIM> nn{nn_[0],nn_[1],nn_[2]};
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
    if (!VolumeTomography()) {
        ::NaNFreeEle(ele,&idx_ele, 0,nn,locstep,*dof, !VolumeTomography());
    } else {
        for (int i2 = 0; i2<nn_[2]-1; ++i2) {
            ::NaNFreeEle(ele,&idx_ele, i2,nn,locstep,*dof, !VolumeTomography());
        }
    }
    // clean the dof by counting to how many elements they belong
    ::CleanDoF(dof,n_dof,&n_dop, *ele,locstep, !VolumeTomography());
    return;
}

// Solve the global minimisation problem.
void EPTGradient::
GlobalMinimisation() {
    const std::array<int,N_DIM> nn{nn_[0],nn_[1],nn_[2]};
    const int ele_nodes = VolumeTomography() ? vox_nodes : quad_nodes;
    const int n_dim = VolumeTomography() ? 3 : 2;
    const std::array<int,N_DIM> step{1,nn_[0],nn_[0]*nn_[1]};
    const int n_vox = GetTxSens(0)->GetNVox();
    const int n_out = sigma_->GetNVox();
    // define the computational elements
    std::vector<int> locstep(ele_nodes,0);
    for (int locid = 0; locid<ele_nodes; ++locid) {
        for (int d = 0; d<n_dim; ++d) {
            locstep[locid] += locid&(1<<d) ? step[d] : 0;
        }
    }
    // select the degrees of freedom
    int n_ele = (nn_[0]-1)*(nn_[1]-1);
    if (VolumeTomography()) {
        n_ele *= (nn_[2]-1);
    }
    std::vector<int> dof(n_out);
    std::vector<int> ele(n_ele);
    int ndof;
    SelectDoF(&dof,&ele,&ndof,locstep);
    // set the dirichlet values in case of seed points...
    std::vector<std::complex<double> > dir(0);
    if (SeedPointsAreUsed()) {
        dir.resize(seed_points_->size());
        ::SeedPointsToDoF(&dir,&dof,&ndof, seed_points_.value(), nn, omega_, !VolumeTomography());
    } else {
        // ... otherwise set the mask
        mask_.resize(ele.size(),false);
        double threshold = ::SetGradientThreshold(gradient_tolerance_,g_plus_,g_z_, !VolumeTomography());
        ::FillGradientMask(&mask_, threshold,ele,locstep,g_plus_,g_z_, !VolumeTomography());
    }
    // create the matrix of coefficients and the forcing term
    Eigen::SparseMatrix<std::complex<double> > A(ndof,ndof);
    Eigen::SparseMatrix<std::complex<double> > B(ndof,ndof);
    Eigen::SparseMatrix<std::complex<double> > Bout(ndof,ndof); // Bout accounts for the masked out mass contributions
    Eigen::VectorXcd b;
    Eigen::VectorXcd b_dir;
    Eigen::VectorXcd l;
    ::FillGlobalMatrix(&A,&B,&Bout,&b,&b_dir,&l, dof,ndof,ele,locstep,epsc_,g_plus_,g_z_,
        dd_, !VolumeTomography(), SeedPointsAreUsed(),dir,mask_);
    // solve the linear system
    Eigen::VectorXcd x(ndof);
    if (!SeedPointsAreUsed()) {
        ::SolveGlobalSystem(&x, A+lambda_*B,b+lambda_*l);
    } else {
        ::SolveGlobalSystem(&x, A,b+b_dir);
    }
    // evaluate the cost functional terms
    ::EvaluateCostFunctional(&cost_functional_,&cost_regularization_, x,A,B,Bout,b,
        dof,ndof,epsc_,g_plus_,g_z_, !VolumeTomography());
    // extract electric properties
    for (int idx = 0; idx<n_out; ++idx) {
        if (dof[idx]>0) {
            epsc_[idx] = std::exp(x[dof[idx]-1]);
        } else {
            if (SeedPointsAreUsed() && dof[idx]<0) {
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
        epsr_->At(idx) = epsc_[idx].real()/EPS0;
        sigma_->At(idx) = -epsc_[idx].imag()*omega_;
    }
    return;
}
