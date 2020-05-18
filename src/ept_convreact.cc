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

#include <iostream>
#include <fstream>

#include <Eigen/Sparse>

using namespace eptlib;

// EPTConvReact constructor
EPTConvReact::
EPTConvReact(const double freq, const std::array<int,NDIM> &nn,
    const std::array<double,NDIM> &dd, const Shape &shape) :
    EPTInterface(freq,nn,dd, 1,1), diff_coeff_(0.0), thereis_diff_(false), fd_filter_(shape) {
    return;
}

// EPTConvReact destructor
EPTConvReact::
~EPTConvReact() {
    return;
}

// EPTConvReact run
EPTlibError_t EPTConvReact::
Run() {
    if (thereis_tx_sens_.all()) {
        thereis_epsr_ = true;
        epsr_.resize(n_vox_);
    }
    if (thereis_trx_phase_.all()) {
        thereis_sigma_ = true;
        sigma_.resize(n_vox_);
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
EPTlibError_t EPTConvReact::
SetVolumeTomography() {
    is_volume_ = true;
    return EPTlibError::Success;
}
// EPTConvReact unset volume tomography
EPTlibError_t EPTConvReact::
UnsetVolumeTomography() {
    is_volume_ = false;
    return EPTlibError::Success;
}

// EPTConvReact set artificial diffusion
EPTlibError_t EPTConvReact::
SetArtificialDiffusion(const double diff_coeff) {
    diff_coeff_ = diff_coeff;
    thereis_diff_ = true;
    return EPTlibError::Success;
}
// EPTConvReact unset artificial diffusion
EPTlibError_t EPTConvReact::
UnsetArtificialDiffusion() {
    if (!thereis_diff_) {
        return EPTlibError::Unknown;
    }
    thereis_diff_ = false;
    return EPTlibError::Success;
}

// EPTConvReact complete EPT
EPTlibError_t EPTConvReact::
CompleteEPTConvReact() {
    double dir_epsr = 1.0;
    double dir_sigma = 0.0;
    std::complex<double> dir_x = 1.0/std::complex<double>(dir_epsr*EPS0,-dir_sigma/omega_);
    // compute the gradient
    std::vector<std::complex<double> > tx_sens_c(n_vox_);
    std::array<std::vector<std::complex<double> >,NDIM> beta;
    for (int idx = 0; idx<n_vox_; ++idx) {
        tx_sens_c[idx] = tx_sens_[0][idx]*std::exp(std::complex<double>(0.0,0.5*trx_phase_[0][idx]));
    }
    for (int d = 0; d<NDIM; ++d) {
        beta[d].resize(n_vox_,NAN);
        DifferentialOperator_t diff_op = static_cast<DifferentialOperator_t>(d);
        EPTlibError_t error = fd_filter_.Apply(diff_op, beta[0].data(),tx_sens_c.data(),nn_,dd_);
        if (error!=EPTlibError::Success) {
            return error;
        }
    }
    for (int idx = 0; idx<n_vox_; ++idx) {
        beta[0][idx] = beta[0][idx] - std::complex<double>(0.0,1.0)*beta[1][idx];
        beta[1][idx] = std::complex<double>(0.0,1.0)*beta[0][idx];
        beta[2][idx] = beta[2][idx];
    }
    // select the degrees of freedom
    std::vector<int> dof(n_vox_);
    int n_dof = 0;
    int n_dop = 0;

    // To be written in a universal way....
    {
        int idx = 0;
        for (int i2 = 0; i2<nn_[2]; ++i2) {
            for (int i1 = 0; i1<nn_[1]; ++i1) {
                for (int i0 = 0; i0<nn_[0]; ++i0) {
                    if (i0<2 || i0>nn_[0]-3 || i1<2 || i1>nn_[1]-3) {
                        dof[idx++] = --n_dop;
                    } else {
                        if (beta[0][idx]==beta[0][idx]) {
                            dof[idx++] = ++n_dof;
                        } else {
                            dof[idx++] = --n_dop;
                        }    
                    }
                }
            }
        }
    }

//    for (int idx = 0; idx<n_vox_; ++idx) {
//        if (beta[0][idx]==beta[0][idx]) {
//            dof[idx] = ++n_dof;
//        } else {
//            dof[idx] = --n_dop;
//        }        
//    }
    
    // build coefficient matrix and forcing term
    const std::array<int,NDIM> step{1,nn_[0],nn_[0]*nn_[1]};
    Eigen::SparseMatrix<std::complex<double> > A(n_dof,n_dof);
    std::vector<Eigen::Triplet<std::complex<double> > > A_trip(0);
    Eigen::VectorXcd b(n_dof);
    for (int idx = 0; idx<n_vox_; ++idx) {
        int idof = dof[idx];
        if (idof > 0) {
            b[idof-1] = 0.0;
            // coefficient matrix
            for (int d = 0; d<NDIM; ++d) {
                if (thereis_diff_) {
                    A_trip.push_back(Eigen::Triplet<std::complex<double> >(idof-1,idof-1,2.0*diff_coeff_/dd_[d]/dd_[d]));
                }
                int jdx = idx+step[d];
                int jdof = dof[jdx];
                if (jdof > 0) {
                    A_trip.push_back(Eigen::Triplet<std::complex<double> >(idof-1,jdof-1,beta[d][jdx]/2.0/dd_[d]));
                    if (thereis_diff_) {
                        A_trip.push_back(Eigen::Triplet<std::complex<double> >(idof-1,jdof-1,-diff_coeff_/dd_[d]/dd_[d]));
                    }
                } else {

                    if (d<2) {

                        b[idof-1] += -dir_x*beta[d][jdx]/2.0/dd_[d];
                        if (thereis_diff_) {
                            b[idof-1] += dir_x*diff_coeff_/dd_[d]/dd_[d];
                        }

                    }

                }
                jdx = idx-step[d];
                jdof = dof[jdx];
                if (jdof > 0) {
                    A_trip.push_back(Eigen::Triplet<std::complex<double> >(idof-1,jdof-1,-beta[d][jdx]/2.0/dd_[d]));
                    if (thereis_diff_) {
                        A_trip.push_back(Eigen::Triplet<std::complex<double> >(idof-1,jdof-1,-diff_coeff_/dd_[d]/dd_[d]));
                    }
                } else {

                    if (d<2) {

                        b[idof-1] += dir_x*beta[d][jdx]/2.0/dd_[d];
                        if (thereis_diff_) {
                            b[idof-1] += dir_x*diff_coeff_/dd_[d]/dd_[d];
                        }

                    }

                }
            }
            // forcing term
            b[idof-1] += -omega_*omega_*MU0*tx_sens_c[idx];
        }
    }
    A.setFromTriplets(A_trip.begin(),A_trip.end());
    // Solve the linear system
//    Eigen::BiCGSTAB<Eigen::SparseMatrix<std::complex<double> > > solver;
//    solver.compute(A);
    Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double> > > solver;
    A.makeCompressed();
    solver.analyzePattern(A);
    solver.factorize(A);
    Eigen::VectorXcd x = solver.solve(b);
    // Extract the electric properties from the result
    for (int idx = 0; idx<n_vox_; ++idx) {
        int idof = dof[idx];
        if (idof > 0) {
            std::complex<double> epsc = 1.0/x[idof-1];
            epsr_[idx] = epsc.real()/EPS0;
            sigma_[idx] = -epsc.imag()*omega_;
        } else {
            epsr_[idx] = dir_epsr;
            sigma_[idx] = dir_sigma;
        }
    }
    return EPTlibError::Success;
}

// EPTConvReact phase-based EPT
EPTlibError_t EPTConvReact::
PhaseEPTConvReact() {
    return EPTlibError::Success;
}
