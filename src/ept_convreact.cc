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

#include "eptlib/ept_convreact.h"

#include <complex>
#include <limits>

#include <Eigen/Sparse>

using namespace eptlib;

// EPTConvReact constructor
EPTConvReact::
EPTConvReact(const size_t n0, const size_t n1, const size_t n2,
    const double d0, const double d1, const double d2,
    const double freq, const Shape &window, const int degree,
    const size_t max_iterations, const double tolerance,
    const double weight_param) :
    EPTInterface(n0,n1,n2, d0,d1,d2, freq, 1,1, false),
    slice_index_(),
    artificial_diffusion_(),
    dirichlet_epsr_(1.0),
    dirichlet_sigma_(0.0),
    sg_filter_(),
    sg_window_(window),
    sg_degree_(degree),
    weight_param_(weight_param),
    solver_max_iterations_(max_iterations),
    solver_tolerance_(tolerance),
    solver_iterations_(0),
    solver_residual_(0.0) {
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
    // check that the sg filter is correctly set up
    if (ThereIsReferenceImage() && !std::holds_alternative<filter::AnatomicalSavitzkyGolay>(sg_filter_)) {
        sg_filter_.emplace<filter::AnatomicalSavitzkyGolay>(dd_[0],dd_[1],dd_[2], sg_window_, sg_degree_, weight_param_);
    } else
    if (!ThereIsReferenceImage() && !std::holds_alternative<filter::SavitzkyGolay>(sg_filter_)) {
        sg_filter_.emplace<filter::SavitzkyGolay>(dd_[0],dd_[1],dd_[2], sg_window_, sg_degree_);
    }
    // check that the essential input is available
    if (!ThereIsTRxPhase(0,0)) {
        return EPTlibError::MissingData;
    }
    // setup the output variables and perform ept
    sigma_ = std::make_unique<Image<double> >(nn_[0], nn_[1], VolumeTomography() ? nn_[2] : 1);
    if (ThereIsTxSens(0)) {
        epsr_ = std::make_unique<Image<double> >(nn_[0], nn_[1], VolumeTomography() ? nn_[2] : 1);
        return CompleteEPTConvReact();
    } else {
        return PhaseEPTConvReact();
    }
}

namespace { // details

    template <typename Scalar>
    void FillDoF(std::vector<int> *dof, int *idx_dof, int *n_dof, int *n_dop,
        const int n_dim, const size_t i2, const std::array<std::ptrdiff_t, N_DIM> &step,
        const std::array<std::size_t, N_DIM> &nn, const Image<Scalar> &beta) {
        std::array<std::size_t, N_DIM> ii{0, 0, i2};
        int idx = step[2]*ii[2];
        for (ii[0] = 0; ii[0]<nn[0]; ++ii[0]) {
            for (ii[1] = 0; ii[1]<nn[1]; ++ii[1]) {
                if (beta(idx)==beta(idx)) {
                    // if beta is not a NaN, it could be a DoF...
                    (*dof)[(*idx_dof)++] = ++(*n_dof);
                    for (size_t d = 0; d<n_dim; ++d) {
                        if (ii[d]==0 || ii[d]==nn[d]-1 || beta(idx+step[d])!=beta(idx+step[d]) || beta(idx-step[d])!=beta(idx-step[d])) {
                            // ...unless it is near to a boundary or a NaN
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
        }
        return;
    }

}  // namespace

// EPTConvReact complete EPT
EPTlibError EPTConvReact::
CompleteEPTConvReact() {
    const std::array<std::ptrdiff_t, N_DIM> step{1, nn_[0], nn_[0]*nn_[1]};
    size_t n_dim = VolumeTomography() ? N_DIM : N_DIM-1;
    size_t n_vox = GetTRxPhase(0,0)->GetNVox();
    size_t n_out = sigma_->GetNVox();
    std::complex<double> dirichlet_x = 1.0/std::complex<double>(dirichlet_epsr_*EPS0, -dirichlet_sigma_/omega_);
    // compute the gradient
    Image<std::complex<double> > tx_sens_c(nn_[0], nn_[1], nn_[2]);
    std::array<Image<std::complex<double> >, N_DIM> beta;
    for (size_t idx = 0; idx<n_vox; ++idx) {
        std::complex<double> exponent = std::complex<double>(0.0, 0.5 * GetTRxPhase(0,0)->At(idx));
        tx_sens_c(idx) = GetTxSens(0)->At(idx) * std::exp(exponent);
    }
    for (size_t d = 0; d<n_dim; ++d) {
        beta[d] = Image<std::complex<double> >(nn_[0], nn_[1], nn_[2]);
        beta[d].GetData().assign(n_vox, nand);
        DifferentialOperator diff_op = static_cast<DifferentialOperator>(d+1);
        EPTlibError error;
        if (ThereIsReferenceImage()) {
            error = std::get<filter::AnatomicalSavitzkyGolay>(sg_filter_).Apply(diff_op, &beta[d], tx_sens_c, *reference_image_);
        } else {
            error = std::get<filter::SavitzkyGolay>(sg_filter_).Apply(diff_op, &beta[d], tx_sens_c);
        }
        if (error!=EPTlibError::Success) {
            return error;
        }
    }
    for (size_t idx = 0; idx<n_vox; ++idx) {
        beta[0](idx) = beta[0](idx) - std::complex<double>(0.0, 1.0) * beta[1](idx);
        beta[1](idx) = std::complex<double>(0.0, 1.0) * beta[0](idx);
    }
    if (!VolumeTomography()) {
        beta[2] = Image<std::complex<double> >(nn_[0], nn_[1], nn_[2]);
        beta[2].GetData().assign(n_vox, nand);
        DifferentialOperator diff_op = DifferentialOperator::GradientZZ;
        EPTlibError error;
        if (ThereIsReferenceImage()) {
            error = std::get<filter::AnatomicalSavitzkyGolay>(sg_filter_).Apply(diff_op, &beta[2], tx_sens_c, *reference_image_);
        } else {
            error = std::get<filter::SavitzkyGolay>(sg_filter_).Apply(diff_op, &beta[2], tx_sens_c);
        }
        if (error!=EPTlibError::Success) {
            return error;
        }
    }
    // select the degrees of freedom
    std::vector<int> dof(n_out);
    int idx_dof = 0;
    int n_dof = 0;
    int n_dop = 0;
    if (VolumeTomography()) {
        for (size_t i2 = 0; i2<nn_[2]; ++i2) {
            ::FillDoF(&dof,&idx_dof,&n_dof,&n_dop, n_dim,i2,step,nn_,beta[0]);
        }
    } else {
        ::FillDoF(&dof,&idx_dof,&n_dof,&n_dop, n_dim,slice_index_.value(),step,nn_,beta[0]);
    }
    // build coefficient matrix and forcing term
    Eigen::SparseMatrix<std::complex<double> > A(n_dof,n_dof);
    std::vector<Eigen::Triplet<std::complex<double> > > A_trip(0);
    Eigen::VectorXcd b(n_dof);
    for (size_t idx_out = 0; idx_out<n_out; ++idx_out) {
        int idx = VolumeTomography() ? idx_out : idx_out+step[2]*slice_index_.value();
        int idof = dof[idx_out];
        if (idof > 0) {
            b[idof-1] = 0.0;
            // coefficient matrix
            for (size_t d = 0; d<n_dim; ++d) {
                if (ThereIsArtificialDiffusion()) {
                    A_trip.push_back(Eigen::Triplet<std::complex<double> >(idof-1, idof-1, 2.0*artificial_diffusion_.value()/dd_[d]/dd_[d]));
                }
                int jdx_out = idx_out+step[d];
                int jdx = idx+step[d];
                int jdof = dof[jdx_out];
                if (jdof > 0) {
                    A_trip.push_back(Eigen::Triplet<std::complex<double> >(idof-1,jdof-1,beta[d](jdx)/2.0/dd_[d]));
                    if (ThereIsArtificialDiffusion()) {
                        A_trip.push_back(Eigen::Triplet<std::complex<double> >(idof-1, jdof-1, -artificial_diffusion_.value()/dd_[d]/dd_[d]));
                    }
                } else {
                    b[idof-1] += -dirichlet_x*beta[d](jdx)/2.0/dd_[d];
                    if (ThereIsArtificialDiffusion()) {
                        b[idof-1] += dirichlet_x*artificial_diffusion_.value()/dd_[d]/dd_[d];
                    }
                }
                jdx_out = idx_out-step[d];
                jdx = idx-step[d];
                jdof = dof[jdx_out];
                if (jdof > 0) {
                    A_trip.push_back(Eigen::Triplet<std::complex<double> >(idof-1,jdof-1,-beta[d](jdx)/2.0/dd_[d]));
                    if (ThereIsArtificialDiffusion()) {
                        A_trip.push_back(Eigen::Triplet<std::complex<double> >(idof-1, jdof-1, -artificial_diffusion_.value()/dd_[d]/dd_[d]));
                    }
                } else {
                    b[idof-1] += dirichlet_x*beta[d](jdx)/2.0/dd_[d];
                    if (ThereIsArtificialDiffusion()) {
                        b[idof-1] += dirichlet_x*artificial_diffusion_.value()/dd_[d]/dd_[d];
                    }
                }
            }
            if (!VolumeTomography()) {
                A_trip.push_back(Eigen::Triplet<std::complex<double> >(idof-1,idof-1,beta[2](idx)));
            }
            // forcing term
            b[idof-1] += -omega_*omega_*MU0*tx_sens_c(idx);
        }
    }
    A.setFromTriplets(A_trip.begin(),A_trip.end());
    A.makeCompressed();
    // Solve the linear system
    Eigen::BiCGSTAB<Eigen::SparseMatrix<std::complex<double> > > solver;
    solver.compute(A);
    solver.setMaxIterations(solver_max_iterations_);
    solver.setTolerance(solver_tolerance_);
    Eigen::VectorXcd x = solver.solve(b);
    solver_iterations_ = solver.iterations();
    solver_residual_ = solver.error();
    // Extract the electric properties from the result
    for (size_t idx = 0; idx<n_out; ++idx) {
        int idof = dof[idx];
        if (idof > 0) {
            std::complex<double> epsc = 1.0/x[idof-1];
            epsr_ ->At(idx) =  epsc.real()/EPS0;
            sigma_->At(idx) = -epsc.imag()*omega_;
        } else {
            epsr_ ->At(idx) = nand;
            sigma_->At(idx) = nand;
        }
    }
    return EPTlibError::Success;
}

// EPTConvReact phase-based EPT
EPTlibError EPTConvReact::
PhaseEPTConvReact() {
    const std::array<std::ptrdiff_t, N_DIM> step{1, nn_[0], nn_[0]*nn_[1]};
    size_t n_dim = VolumeTomography() ? N_DIM : N_DIM-1;
    size_t n_vox = GetTRxPhase(0,0)->GetNVox();
    size_t n_out = sigma_->GetNVox();
    double dir_x = 1.0/dirichlet_sigma_;
    // compute the gradient
    std::array<Image<double>, N_DIM> beta;
    for (size_t d = 0; d<n_dim; ++d) {
        beta[d] = Image<double>(nn_[0], nn_[1], nn_[2]);
        beta[d].GetData().assign(n_vox, nand);
        DifferentialOperator diff_op = static_cast<DifferentialOperator>(d+1);
        EPTlibError error;
        if (PhaseIsWrapped()) {
            if (ThereIsReferenceImage()) {
                throw std::runtime_error("Feature not implemented!");
            } else {
                error = std::get<filter::SavitzkyGolay>(sg_filter_).ApplyWrappedPhase(diff_op, &beta[d], *GetTRxPhase(0,0));
            }
        } else {
            if (ThereIsReferenceImage()) {
                error = std::get<filter::AnatomicalSavitzkyGolay>(sg_filter_).Apply(diff_op, &beta[d], *GetTRxPhase(0,0), *reference_image_);
            } else {
                error = std::get<filter::SavitzkyGolay>(sg_filter_).Apply(diff_op, &beta[d], *GetTRxPhase(0,0));
            }
        }
        if (error!=EPTlibError::Success) {
            return error;
        }
    }
    if (!VolumeTomography()) {
        beta[2] = Image<double>(nn_[0], nn_[1], nn_[2]);
        beta[2].GetData().assign(n_vox, nand);
        DifferentialOperator diff_op = DifferentialOperator::GradientZZ;
        EPTlibError error;
        if (PhaseIsWrapped()) {
            if (ThereIsReferenceImage()) {
                throw std::runtime_error("Feature not implemented!");
            } else {
                error = std::get<filter::SavitzkyGolay>(sg_filter_).ApplyWrappedPhase(diff_op, &beta[2], *GetTRxPhase(0,0));
            }
        } else {
            if (ThereIsReferenceImage()) {
                error = std::get<filter::AnatomicalSavitzkyGolay>(sg_filter_).Apply(diff_op, &beta[2], *GetTRxPhase(0,0), *reference_image_);
            } else {
                error = std::get<filter::SavitzkyGolay>(sg_filter_).Apply(diff_op, &beta[2], *GetTRxPhase(0,0));
            }
        }
        if (error!=EPTlibError::Success) {
            return error;
        }
    }
    // select the degrees of freedom
    std::vector<int> dof(n_out);
    int idx_dof = 0;
    int n_dof = 0;
    int n_dop = 0;
    if (VolumeTomography()) {
        for (size_t i2 = 0; i2<nn_[2]; ++i2) {
            ::FillDoF(&dof,&idx_dof,&n_dof,&n_dop, n_dim,i2,step,nn_,beta[0]);
        }
    } else {
        ::FillDoF(&dof,&idx_dof,&n_dof,&n_dop, n_dim,slice_index_.value(),step,nn_,beta[0]);
    }
    // build coefficient matrix and forcing term
    Eigen::SparseMatrix<double> A(n_dof,n_dof);
    std::vector<Eigen::Triplet<double> > A_trip(0);
    Eigen::VectorXd b(n_dof);
    for (size_t idx_out = 0; idx_out<n_out; ++idx_out) {
        int idx = VolumeTomography() ? idx_out : idx_out+step[2]*slice_index_.value();
        int idof = dof[idx_out];
        if (idof > 0) {
            b[idof-1] = 0.0;
            // coefficient matrix
            for (size_t d = 0; d<n_dim; ++d) {
                // diffusion (central term)
                if (ThereIsArtificialDiffusion()) {
                    A_trip.push_back(Eigen::Triplet<double>(idof-1, idof-1, 2.0*artificial_diffusion_.value()/dd_[d]/dd_[d]));
                }
                std::array<int,2> sides{+1,-1};
                for (size_t s = 0; s<2; ++s) {
                    int jdx_out = idx_out+sides[s]*step[d];
                    int jdx = idx+sides[s]*step[d];
                    int jdof = dof[jdx_out];
                    // convection
                    if (beta[d](idx)*beta[d](jdx)>0) {
                        if (sides[s]*beta[d](idx)>0) {
                            double A_tmp = sides[s]*beta[d](idx)/dd_[d];
                            A_trip.push_back(Eigen::Triplet<double>(idof-1,idof-1,A_tmp));
                        } else {
                            double A_tmp = sides[s]*beta[d](jdx)/dd_[d];
                            if (jdof > 0) {
                                A_trip.push_back(Eigen::Triplet<double>(idof-1,jdof-1,A_tmp));
                            } else {
                                b[idof-1] += -dir_x*A_tmp;
                            }
                        }
                    } else {
                        double A_tmp = sides[s]*beta[d](idx)/dd_[d]/2.0;
                        A_trip.push_back(Eigen::Triplet<double>(idof-1,idof-1,A_tmp));
                        A_tmp = sides[s]*beta[d](jdx)/dd_[d]/2.0;
                        if (jdof > 0) {
                            A_trip.push_back(Eigen::Triplet<double>(idof-1,jdof-1,A_tmp));
                        } else {
                            b[idof-1] += -dir_x*A_tmp;
                        }
                    }
                    // diffusion
                    if (ThereIsArtificialDiffusion()) {
                        if (jdof > 0) {
                            A_trip.push_back(Eigen::Triplet<double>(idof-1, jdof-1, -artificial_diffusion_.value()/dd_[d]/dd_[d]));
                        } else {
                            b[idof-1] += dir_x*artificial_diffusion_.value()/dd_[d]/dd_[d];
                        }
                    }
                }
            }
            if (!VolumeTomography()) {
                A_trip.push_back(Eigen::Triplet<double>(idof-1,idof-1,beta[2](idx)));
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
    for (size_t idx = 0; idx<n_out; ++idx) {
        int idof = dof[idx];
        if (idof > 0) {
            sigma_->At(idx) = 1.0/x[idof-1];
        } else {
            sigma_->At(idx) = nand;
        }
    }
    return EPTlibError::Success;
}
