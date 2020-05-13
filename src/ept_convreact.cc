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

#include <Eigen/Sparse>

using namespace eptlib;

// EPTConvReact constructor
EPTConvReact::
EPTConvReact(const double freq, const std::array<int,NDIM> &nn,
    const std::array<double,NDIM> &dd, const Shape &shape) :
    EPTInterface(freq,nn,dd, 1,1), fd_filter_(shape) {
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

// EPTConvReact complete EPT
EPTlibError_t EPTConvReact::
CompleteEPTConvReact() {
    // compute the gradient
    std::vector<std::complex<double> > tx_sens_c(n_vox_);
    std::array<std::vector<std::complex<double> >,NDIM> beta;
    for (int idx = 0; idx<n_vox_; ++idx) {
        tx_sens_c[idx] = tx_sens_[0][idx]*std::exp(std::complex<double>(0.0,0.5*trx_phase_[0][idx]));
    }
    for (int d = 0; d<NDIM; ++d) {
        beta[d].resize(n_vox_,NAN);
        EPTlibError_t error = fd_filter_.ComputeGradient(d,beta[d].data(),tx_sens_c.data(),nn_,dd_);
        if (error!=EPTlibError::Success) {
            return error;
        }
    }
    for (int idx = 0; idx<n_vox_; ++idx) {
        beta[0][idx] = beta[0][idx] - std::complex<double>(0.0,1.0)*beta[1][idx];
        beta[1][idx] = -std::complex<double>(0.0,1.0)*beta[0][idx];
        beta[2][idx] = beta[2][idx];
    }
    // select the degrees of freedom
    std::vector<int> dof(n_vox_);
    int n_dof = 0;
    int n_dop = 0;
    for (int idx = 0; idx<n_vox_; ++idx) {
        if (beta[0][idx]==beta[0][idx]) {
            dof[idx] = ++n_dof;
        } else {
            dof[idx] = --n_dop;
        }        
    }
    // build coefficient matrix and forcing term
    const std::array<int,NDIM> step{1,nn_[0],nn_[0]*nn_[1]};
    Eigen::SparseMatrix<std::complex<double> > A(n_dof,n_dof);
    std::vector<Eigen::Triplet<std::complex<double> > > A_trip(0);
    Eigen::VectorXcd b(n_dof);
    for (int idx = 0; idx<n_vox_; ++idx) {
        int idof = dof[idx];
        if (idof > 0) {
            // coefficient matrix
            for (int d = 0; d<NDIM; ++d) {
                int jdx = idx+step[d];
                int jdof = dof[jdx];
                if (jdof > 0) {
                    A_trip.push_back(Eigen::Triplet<std::complex<double> >(idof-1,jdof-1,beta[d][jdx]/2.0/dd_[d]));
                }
                jdx = idx-step[d];
                jdof = dof[jdx];
                if (jdof > 0) {
                    A_trip.push_back(Eigen::Triplet<std::complex<double> >(idof-1,jdof-1,-beta[d][jdx]/2.0/dd_[d]));
                }
            }
            // forcing term
            b[idof-1] = -omega_*omega_*MU0*tx_sens_c[idx];
        }
    }
    A.setFromTriplets(A_trip.begin(),A_trip.end());
    // Solve the linear system
    Eigen::BiCGSTAB<Eigen::SparseMatrix<std::complex<double> > > solver;
    solver.compute(A);
    Eigen::VectorXcd x = solver.solve(b);
    // Extract the electric properties from the result
    for (int idx = 0; idx<n_vox_; ++idx) {
        int idof = dof[idx];
        if (idof > 0) {
            std::complex<double> epsc = 1.0/x[idof-1];
            epsr_[idx] = epsc.real()/EPS0;
            sigma_[idx] = -epsc.imag()*omega_;
        } else {
            epsr_[idx] = 0.0;
            sigma_[idx] = 0.0;
        }
    }
    return EPTlibError::Success;
}

// EPTConvReact phase-based EPT
EPTlibError_t EPTConvReact::
PhaseEPTConvReact() {
    return EPTlibError::Success;
}
