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

#include "eptlib/ept_convreact.h"

#include <complex>
#include <limits>

#include <iostream>
#include <fstream>

#include <Eigen/Sparse>

using namespace eptlib;

// EPTConvReact constructor
EPTConvReact::
EPTConvReact(const double freq, const std::array<int,NDIM> &nn,
    const std::array<double,NDIM> &dd, const Shape &shape) :
    EPTInterface(freq,nn,dd, 1,1), dir_epsr_(1.0), dir_sigma_(0.0),
    plane_idx_(nn[2]/2), is_volume_(false), diff_coeff_(0.0),
    thereis_diff_(false), fd_filter_(shape),
    solver_iterations_(0), solver_residual_(0.0) {
    return;
}

// EPTConvReact destructor
EPTConvReact::
~EPTConvReact() {
    return;
}

// EPTConvReact run
EPTlibError EPTConvReact::
Run() {
    if (thereis_tx_sens_.all()) {
        thereis_epsr_ = true;
        if (is_volume_) {
            epsr_ = Image<double>(nn_[0],nn_[1],nn_[2]);
        } else {
            epsr_ = Image<double>(nn_[0],nn_[1]);
        }
    }
    if (thereis_trx_phase_.all()) {
        thereis_sigma_ = true;
        if (is_volume_) {
            sigma_ = Image<double>(nn_[0],nn_[1],nn_[2]);
        } else {
            sigma_ = Image<double>(nn_[0],nn_[1]);
        }
    } else {
        return EPTlibError::MissingData;
    }
    if (thereis_epsr_) {
        // ...complete convection-reaction EPT
        return CompleteEPTConvReact();
    } else {
        // ...phase-based convection-reaction EPT
        return PhaseEPTConvReact();
    }
    return EPTlibError::Success;
}

// EPTConvReact set volume tomography
bool EPTConvReact::
Toggle3D() {
    is_volume_ = !is_volume_;
    return is_volume_;
}

// EPTConvReact set the selected plane index
EPTlibError EPTConvReact::
SelectSlice(const int slice_idx) {
    int r2 = fd_filter_.GetShape().GetSize()[2]/2;
    if (slice_idx<r2||slice_idx>nn_[2]-1-r2) {
        return EPTlibError::OutOfRange;
    }
    plane_idx_ = slice_idx;
    return EPTlibError::Success;
}

// EPTConvReact set the dirichlet boundary condition
EPTlibError EPTConvReact::
SetDirichlet(const double dir_epsr, const double dir_sigma) {
    if (dir_epsr<1.0||dir_sigma<0.0) {
        return EPTlibError::WrongDataFormat;
    }
    dir_epsr_ = dir_epsr;
    dir_sigma_ = dir_sigma;
    return EPTlibError::Success;
}

// EPTConvReact set artificial diffusion
void EPTConvReact::
SetArtificialDiffusion(const double diff_coeff) {
    diff_coeff_ = diff_coeff;
    thereis_diff_ = true;
    return;
}
// EPTConvReact unset artificial diffusion
void EPTConvReact::
UnsetArtificialDiffusion() {
    thereis_diff_ = false;
    return;
}

// EPTConvReact get the number of iterations
int EPTConvReact::
GetSolverIterations() {
    return solver_iterations_;
}
// EPTConvReact get estimated error
double EPTConvReact::
GetSolverResidual() {
    return solver_residual_;
}

namespace { // details
    template <typename T>
    void FillDoF(std::vector<int> *dof, int *idx_dof, int *n_dof, int *n_dop,
        const int n_dim, const int i2, const std::array<int,NDIM> &step,
        const std::vector<T> &beta) {
        int idx = step[2]*i2;
        for (int idx_out = 0; idx_out<step[2]; ++idx_out) {
            if (beta[idx]==beta[idx]) {
                // if beta is not a NaN, it could be a DoF...
                (*dof)[(*idx_dof)++] = ++(*n_dof);
                for (int d = 0; d<n_dim; ++d) {
                    if (beta[idx+step[d]]!=beta[idx+step[d]] || beta[idx-step[d]]!=beta[idx-step[d]]) {
                        // ...unless it is near to a NaN
                        (*dof)[*idx_dof-1] = --(*n_dop);
                        --(*n_dof);
                        break;
                    }
                }
            } else {
                (*dof)[(*idx_dof)++] = --(*n_dop);
            }
            ++idx;
        }
        return;
    }
}  // namespace

// EPTConvReact complete EPT
EPTlibError EPTConvReact::
CompleteEPTConvReact() {
    const std::array<int,NDIM> step{1,nn_[0],nn_[0]*nn_[1]};
    int n_dim = is_volume_?NDIM:NDIM-1;
    int n_out = is_volume_?n_vox_:step[2];
    std::complex<double> dir_x = 1.0/std::complex<double>(dir_epsr_*EPS0,-dir_sigma_/omega_);
    // compute the gradient
    std::vector<std::complex<double> > tx_sens_c(n_vox_);
    std::array<std::vector<std::complex<double> >,NDIM> beta;
    for (int idx = 0; idx<n_vox_; ++idx) {
        tx_sens_c[idx] = (*tx_sens_[0])[idx]*std::exp(std::complex<double>(0.0,0.5*(*trx_phase_[0])[idx]));
    }
    for (int d = 0; d<n_dim; ++d) {
        beta[d].resize(n_vox_,std::numeric_limits<double>::quiet_NaN());
        DifferentialOperator diff_op = static_cast<DifferentialOperator>(d);
        EPTlibError error = fd_filter_.Apply(diff_op,beta[d].data(),tx_sens_c.data(),nn_,dd_);
        if (error!=EPTlibError::Success) {
            return error;
        }
    }
    for (int idx = 0; idx<n_vox_; ++idx) {
        beta[0][idx] = beta[0][idx] - std::complex<double>(0.0,1.0)*beta[1][idx];
        beta[1][idx] = std::complex<double>(0.0,1.0)*beta[0][idx];
    }
    if (!is_volume_) {
        beta[2].resize(n_vox_,std::numeric_limits<double>::quiet_NaN());
        DifferentialOperator diff_op = DifferentialOperator::GradientZZ;
        EPTlibError error = fd_filter_.Apply(diff_op,beta[2].data(),tx_sens_c.data(),nn_,dd_);
        if (error!=EPTlibError::Success) {
            return error;
        }
    }
    // select the degrees of freedom
    std::vector<int> dof(n_out);
    int idx_dof = 0;
    int n_dof = 0;
    int n_dop = 0;
    if (is_volume_) {
        for (int i2 = 0; i2<nn_[2]; ++i2) {
            ::FillDoF(&dof,&idx_dof,&n_dof,&n_dop, n_dim,i2,step,beta[0]);
        }
    } else {
        ::FillDoF(&dof,&idx_dof,&n_dof,&n_dop, n_dim,plane_idx_,step,beta[0]);
    }
    // build coefficient matrix and forcing term
    Eigen::SparseMatrix<std::complex<double> > A(n_dof,n_dof);
    std::vector<Eigen::Triplet<std::complex<double> > > A_trip(0);
    Eigen::VectorXcd b(n_dof);
    for (int idx_out = 0; idx_out<n_out; ++idx_out) {
        int idx = is_volume_?idx_out:idx_out+step[2]*plane_idx_;
        int idof = dof[idx_out];
        if (idof > 0) {
            b[idof-1] = 0.0;
            // coefficient matrix
            for (int d = 0; d<n_dim; ++d) {
                if (thereis_diff_) {
                    A_trip.push_back(Eigen::Triplet<std::complex<double> >(idof-1,idof-1,2.0*diff_coeff_/dd_[d]/dd_[d]));
                }
                int jdx_out = idx_out+step[d];
                int jdx = idx+step[d];
                int jdof = dof[jdx_out];
                if (jdof > 0) {
                    A_trip.push_back(Eigen::Triplet<std::complex<double> >(idof-1,jdof-1,beta[d][jdx]/2.0/dd_[d]));
                    if (thereis_diff_) {
                        A_trip.push_back(Eigen::Triplet<std::complex<double> >(idof-1,jdof-1,-diff_coeff_/dd_[d]/dd_[d]));
                    }
                } else {
                    b[idof-1] += -dir_x*beta[d][jdx]/2.0/dd_[d];
                    if (thereis_diff_) {
                        b[idof-1] += dir_x*diff_coeff_/dd_[d]/dd_[d];
                    }
                }
                jdx_out = idx_out-step[d];
                jdx = idx-step[d];
                jdof = dof[jdx_out];
                if (jdof > 0) {
                    A_trip.push_back(Eigen::Triplet<std::complex<double> >(idof-1,jdof-1,-beta[d][jdx]/2.0/dd_[d]));
                    if (thereis_diff_) {
                        A_trip.push_back(Eigen::Triplet<std::complex<double> >(idof-1,jdof-1,-diff_coeff_/dd_[d]/dd_[d]));
                    }
                } else {
                    b[idof-1] += dir_x*beta[d][jdx]/2.0/dd_[d];
                    if (thereis_diff_) {
                        b[idof-1] += dir_x*diff_coeff_/dd_[d]/dd_[d];
                    }
                }
            }
            if (!is_volume_) {
                A_trip.push_back(Eigen::Triplet<std::complex<double> >(idof-1,idof-1,beta[2][idx]));
            }
            // forcing term
            b[idof-1] += -omega_*omega_*MU0*tx_sens_c[idx];
        }
    }
    A.setFromTriplets(A_trip.begin(),A_trip.end());
    A.makeCompressed();
    // Solve the linear system
    Eigen::BiCGSTAB<Eigen::SparseMatrix<std::complex<double> > > solver;
    solver.compute(A);
    Eigen::VectorXcd x = solver.solve(b);
    solver_iterations_ = solver.iterations();
    solver_residual_ = solver.error();
    // Extract the electric properties from the result
    for (int idx = 0; idx<n_out; ++idx) {
        int idof = dof[idx];
        if (idof > 0) {
            std::complex<double> epsc = 1.0/x[idof-1];
            epsr_[idx] = epsc.real()/EPS0;
            sigma_[idx] = -epsc.imag()*omega_;
        } else {
            epsr_[idx] = dir_epsr_;
            sigma_[idx] = dir_sigma_;
        }
    }
    return EPTlibError::Success;
}

// EPTConvReact phase-based EPT
EPTlibError EPTConvReact::
PhaseEPTConvReact() {
    const std::array<int,NDIM> step{1,nn_[0],nn_[0]*nn_[1]};
    int n_dim = is_volume_?NDIM:NDIM-1;
    int n_out = is_volume_?n_vox_:step[2];
    double dir_x = 1.0/dir_sigma_;
    // compute the gradient
    std::array<std::vector<double>,NDIM> beta;
    for (int d = 0; d<n_dim; ++d) {
        beta[d].resize(n_vox_,std::numeric_limits<double>::quiet_NaN());
        DifferentialOperator diff_op = static_cast<DifferentialOperator>(d);
        EPTlibError error = fd_filter_.Apply(diff_op,beta[d].data(),trx_phase_[0]->GetData().data(),nn_,dd_);
        if (error!=EPTlibError::Success) {
            return error;
        }
    }
    if (!is_volume_) {
        beta[2].resize(n_vox_,std::numeric_limits<double>::quiet_NaN());
        DifferentialOperator diff_op = DifferentialOperator::GradientZZ;
        EPTlibError error = fd_filter_.Apply(diff_op,beta[2].data(),trx_phase_[0]->GetData().data(),nn_,dd_);
        if (error!=EPTlibError::Success) {
            return error;
        }
    }
    // select the degrees of freedom
    std::vector<int> dof(n_out);
    int idx_dof = 0;
    int n_dof = 0;
    int n_dop = 0;
    if (is_volume_) {
        for (int i2 = 0; i2<nn_[2]; ++i2) {
            ::FillDoF(&dof,&idx_dof,&n_dof,&n_dop, n_dim,i2,step,beta[0]);
        }
    } else {
        ::FillDoF(&dof,&idx_dof,&n_dof,&n_dop, n_dim,plane_idx_,step,beta[0]);
    }
    // build coefficient matrix and forcing term
    Eigen::SparseMatrix<double> A(n_dof,n_dof);
    std::vector<Eigen::Triplet<double> > A_trip(0);
    Eigen::VectorXd b(n_dof);
    for (int idx_out = 0; idx_out<n_out; ++idx_out) {
        int idx = is_volume_?idx_out:idx_out+step[2]*plane_idx_;
        int idof = dof[idx_out];
        if (idof > 0) {
            b[idof-1] = 0.0;
            // coefficient matrix
            for (int d = 0; d<n_dim; ++d) {
                // diffusion (central term)
                if (thereis_diff_) {
                    A_trip.push_back(Eigen::Triplet<double>(idof-1,idof-1,2.0*diff_coeff_/dd_[d]/dd_[d]));
                }
                std::array<int,2> sides{+1,-1};
                for (int s = 0; s<2; ++s) {
                    int jdx_out = idx_out+sides[s]*step[d];
                    int jdx = idx+sides[s]*step[d];
                    int jdof = dof[jdx_out];
                    // convection
                    if (beta[d][idx]*beta[d][jdx]>0) {
                        if (sides[s]*beta[d][idx]>0) {
                            double A_tmp = sides[s]*beta[d][idx]/dd_[d];
                            A_trip.push_back(Eigen::Triplet<double>(idof-1,idof-1,A_tmp));
                        } else {
                            double A_tmp = sides[s]*beta[d][jdx]/dd_[d];
                            if (jdof > 0) {
                                A_trip.push_back(Eigen::Triplet<double>(idof-1,jdof-1,A_tmp));
                            } else {
                                b[idof-1] += -dir_x*A_tmp;
                            }
                        }
                    } else {
                        double A_tmp = sides[s]*beta[d][idx]/dd_[d]/2.0;
                        A_trip.push_back(Eigen::Triplet<double>(idof-1,idof-1,A_tmp));
                        A_tmp = sides[s]*beta[d][jdx]/dd_[d]/2.0;
                        if (jdof > 0) {
                            A_trip.push_back(Eigen::Triplet<double>(idof-1,jdof-1,A_tmp));
                        } else {
                            b[idof-1] += -dir_x*A_tmp;
                        }
                    }
                    // diffusion
                    if (thereis_diff_) {
                        if (jdof > 0) {
                            A_trip.push_back(Eigen::Triplet<double>(idof-1,jdof-1,-diff_coeff_/dd_[d]/dd_[d]));
                        } else {
                            b[idof-1] += dir_x*diff_coeff_/dd_[d]/dd_[d];
                        }
                    }
                }
            }
            if (!is_volume_) {
                A_trip.push_back(Eigen::Triplet<double>(idof-1,idof-1,beta[2][idx]));
            }
            // forcing term
            b[idof-1] += 2*omega_*MU0;
        }
    }
    A.setFromTriplets(A_trip.begin(),A_trip.end());
    A.makeCompressed();
    // Solve the linear system
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
    solver.compute(A);
    Eigen::VectorXd x = solver.solve(b);
    solver_iterations_ = solver.iterations();
    solver_residual_ = solver.error();
    // Extract the electric properties from the result
    for (int idx = 0; idx<n_out; ++idx) {
        int idof = dof[idx];
        if (idof > 0) {
            sigma_[idx] = 1.0/x[idof-1];
        } else {
            sigma_[idx] = dir_sigma_;
        }
    }
    return EPTlibError::Success;
}
