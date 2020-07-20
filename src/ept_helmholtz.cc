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

#include <complex>
#include <vector>

#include "eptlib/ept_helmholtz.h"

using namespace eptlib;

// EPTHelmholtz constructor
EPTHelmholtz::
EPTHelmholtz(const double freq, const std::array<int,NDIM> &nn,
    const std::array<double,NDIM> &dd, const Shape &shape) :
    EPTInterface(freq,nn,dd), fd_lapl_(shape) {
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
    if (thereis_tx_sens_.all()) {
        thereis_epsr_ = true;
        epsr_ = Image<double>(nn_[0],nn_[1],nn_[2]);
    }
    if (thereis_trx_phase_.all()) {
        thereis_sigma_ = true;
        sigma_ = Image<double>(nn_[0],nn_[1],nn_[2]);
    }
    if (thereis_epsr_ && thereis_sigma_) {
        // ...complete Helmholtz-based
        CompleteEPTHelm();
    } else if (thereis_epsr_) {
        // ...magnitude-based approximation
        MagnitudeEPTHelm();
    } else if (thereis_sigma_) {
        // ...phase-based approximation
        PhaseEPTHelm();
    } else {
        // ...not enough input data provided
        return EPTlibError::MissingData;
    }
    return EPTlibError::Success;
}

// EPTHelmholtz private methods
void EPTHelmholtz::
CompleteEPTHelm() {
    std::vector<std::complex<double> > tx_sens_c(n_vox_);
    std::vector<std::complex<double> > epsc(n_vox_);
    for (int idx = 0; idx<n_vox_; ++idx) {
        tx_sens_c[idx] = (*tx_sens_[0])[idx]*std::exp(std::complex<double>(0.0,0.5*(*trx_phase_[0])[idx]));
    }
    DifferentialOperator diff_op = DifferentialOperator::Laplacian;
    fd_lapl_.Apply(diff_op,epsc.data(),tx_sens_c.data(),nn_,dd_);
    for (int idx = 0; idx<n_vox_; ++idx) {
        epsc[idx] /= -MU0*omega_*omega_*tx_sens_c[idx];
        epsr_[idx] = std::real(epsc[idx])/EPS0;
        sigma_[idx] = -std::imag(epsc[idx])*omega_;
    }
    return;
}
void EPTHelmholtz::
MagnitudeEPTHelm() {
    DifferentialOperator diff_op = DifferentialOperator::Laplacian;
    fd_lapl_.Apply(diff_op,epsr_.GetData().data(),tx_sens_[0]->GetData().data(),nn_,dd_);
    for (int idx = 0; idx<n_vox_; ++idx) {
        epsr_[idx] /= -EPS0*MU0*omega_*omega_*(*tx_sens_[0])[idx];
    }
    return;
}
void EPTHelmholtz::
PhaseEPTHelm() {
    DifferentialOperator diff_op = DifferentialOperator::Laplacian;
    fd_lapl_.Apply(diff_op,sigma_.GetData().data(),trx_phase_[0]->GetData().data(),nn_,dd_);
    for (int idx = 0; idx<n_vox_; ++idx) {
        sigma_[idx] /= 2.0*MU0*omega_;
    }
    return;
}
