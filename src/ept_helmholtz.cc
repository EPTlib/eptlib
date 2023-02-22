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

#include "eptlib/ept_helmholtz.h"

#include <complex>
#include <limits>
#include <vector>

using namespace eptlib;

namespace {
	static double nand = std::numeric_limits<double>::quiet_NaN();
	static std::complex<double> nancd = std::complex<double>(nand,nand);
}

// EPTHelmholtz constructor
EPTHelmholtz::
EPTHelmholtz(const size_t n0, const size_t n1, const size_t n2,
    const double d0, const double d1, const double d2,
    const double freq, const Shape &shape, const int degree,
    const bool trx_phase_is_wrapped, const bool compute_variance) :
    EPTInterface(n0,n1,n2, d0,d1,d2, freq, 1,1, trx_phase_is_wrapped),
    fd_lapl_(shape,degree),
    variance_(nullptr),
    compute_variance_(compute_variance) {
    return;
}

// EPTHelmholtz destructor
EPTHelmholtz::
~EPTHelmholtz() {
    return;
}

// EPTHelmholtz run
EPTlibError EPTHelmholtz::
Run() {
    // setup the output variables
    if (ThereIsTxSens(0)) {
        epsr_ = std::make_unique<Image<double> >(nn_[0],nn_[1],nn_[2]);
    }
    if (ThereIsTRxPhase(0,0)) {
        sigma_ = std::make_unique<Image<double> >(nn_[0],nn_[1],nn_[2]);
    }
    if (ComputeVariance() && !ThereIsEpsr()) {
        variance_ = std::make_unique<Image<double> >(nn_[0],nn_[1],nn_[2]);
    }
    // perform ept
    if (ThereIsEpsr() && ThereIsSigma()) {
        CompleteEPTHelm();
    } else if (ThereIsEpsr()) {
        MagnitudeEPTHelm();
    } else if (ThereIsSigma()) {
        PhaseEPTHelm();
    } else {
        return EPTlibError::MissingData;
    }
    return EPTlibError::Success;
}

// EPTHelmholtz complete method
void EPTHelmholtz::
CompleteEPTHelm() {
    auto n_vox = sigma_->GetNVox();
    std::array<int,N_DIM> nn{nn_[0],nn_[1],nn_[2]};
    // setup the input
    std::vector<std::complex<double> > tx_sens_c(n_vox);
    std::vector<std::complex<double> > eps_c(n_vox,::nancd);
    for (int idx = 0; idx<n_vox; ++idx) {
        std::complex<double> exponent = std::complex<double>(0.0, 0.5 * GetTRxPhase(0,0)->At(idx));
        tx_sens_c[idx] = GetTxSens(0)->At(idx) * std::exp(exponent);
    }
    // compute the laplacian
    DifferentialOperator diff_op = DifferentialOperator::Laplacian;
    fd_lapl_.Apply(diff_op, eps_c.data(), tx_sens_c.data(), nn, dd_);
    // extract the output
    for (int idx = 0; idx<n_vox; ++idx) {
        eps_c[idx] /= -MU0*omega_*omega_*tx_sens_c[idx];
        epsr_->At(idx) = std::real(eps_c[idx])/EPS0;
        sigma_->At(idx) = -std::imag(eps_c[idx])*omega_;
    }
    return;
}

// EPTHelmholtz magnitude-based method
void EPTHelmholtz::
MagnitudeEPTHelm() {
    auto n_vox = sigma_->GetNVox();
    std::array<int,N_DIM> nn{nn_[0],nn_[1],nn_[2]};
    // setup the input
    epsr_->GetData().assign(n_vox,::nand);
    // compute the laplacian
    DifferentialOperator diff_op = DifferentialOperator::Laplacian;
    fd_lapl_.Apply(diff_op, epsr_->GetData().data(), GetTxSens(0)->GetData().data(), nn, dd_);
    // extract the output
    for (int idx = 0; idx<n_vox; ++idx) {
        epsr_->At(idx) /= -EPS0*MU0*omega_*omega_ * GetTxSens(0)->At(idx);
    }
    return;
}

// EPTHelmholtz phase-based method
void EPTHelmholtz::
PhaseEPTHelm() {
    auto n_vox = sigma_->GetNVox();
    std::array<int,N_DIM> nn{nn_[0],nn_[1],nn_[2]};
    // setup the input
    sigma_->GetData().assign(n_vox,::nand);
    if(ComputeVariance()) {
        variance_->GetData().assign(n_vox,::nand);
    }
    // compute the laplacian
    DifferentialOperator diff_op = DifferentialOperator::Laplacian;
    if (PhaseIsWrapped() && !ComputeVariance()) {
        fd_lapl_.ApplyWrappedPhase(diff_op, sigma_->GetData().data(), GetTRxPhase(0,0)->GetData().data(), nn, dd_);
    } else if (!ComputeVariance()) {
        fd_lapl_.Apply(diff_op, sigma_->GetData().data(), GetTRxPhase(0,0)->GetData().data(), nn, dd_);
    } else {
        fd_lapl_.Apply(diff_op, sigma_->GetData().data(), variance_->GetData().data(), GetTRxPhase(0,0)->GetData().data(), nn, dd_);
    }
    // extract the output
    for (int idx = 0; idx<n_vox; ++idx) {
        sigma_->At(idx) /= 2.0*MU0*omega_;
        if (ComputeVariance()) {
            variance_->At(idx) /= 2.0*MU0*omega_;
        }
    }
    return;
}
