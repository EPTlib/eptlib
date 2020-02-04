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

#include "eptlib/ept_helmholtz.h"

#include <algorithm>
#include <complex>
#include <iostream>

using namespace eptlib;

// EPTHelmholtz destructor.
EPTHelmholtz::
~EPTHelmholtz() {
    return;
}

// EPTHelmholtz run.
EPTlibError_t EPTHelmholtz::
Run() {
    if (thereis_tx_sens_.all()) {
        thereis_epsr_ = true;
        epsr_.resize(n_vox_);
    }
    if (thereis_trx_phase_.all()) {
        thereis_sigma_ = true;
        sigma_.resize(n_vox_);
    }
    if (thereis_tx_sens_.all() && thereis_trx_phase_.all()) {
        // complete Helmholtz-based EPT
        std::vector<std::complex<real_t> > tx_sens_c(n_vox_);
        std::vector<std::complex<real_t> > epsc(n_vox_);
        for (int idx = 0; idx<n_vox_; ++idx) {
            tx_sens_c[idx] = tx_sens_[0][idx]*std::exp(std::complex<real_t>(0.0,0.5*trx_phase_[0][idx]));
        }
        fd_lapl_.ApplyFilter(epsc.data(),tx_sens_c.data(),nn_,dd_);
        std::transform(epsc.begin(),epsc.end(),tx_sens_c.begin(),epsc.begin(),
            [this](const std::complex<real_t> &a, const std::complex<real_t> &b) -> std::complex<real_t> {
                return -a/omega_/omega_/MU0/b;
            });
        std::transform(epsc.begin(),epsc.end(),epsr_.begin(),
            [](const std::complex<real_t> &a) -> real_t {
                return std::real(a)/EPS0;
            });
        std::transform(epsc.begin(),epsc.end(),sigma_.begin(),
            [this](const std::complex<real_t> &a) -> real_t {
                return -std::imag(a)*omega_;
            });

    } else if (thereis_tx_sens_.all()) {
        // magnitude-based approximation
        fd_lapl_.ApplyFilter(epsr_.data(),tx_sens_[0],nn_,dd_);
        std::transform(epsr_.begin(),epsr_.end(),tx_sens_[0],epsr_.begin(),
            [this](const real_t &a, const real_t &b) -> real_t {
                return -a/EPS0/MU0/omega_/omega_/b;
            });

    } else if (thereis_trx_phase_.all()) {
        // phase-based approximation
        fd_lapl_.ApplyFilter(sigma_.data(),trx_phase_[0],nn_,dd_);
        std::transform(sigma_.begin(),sigma_.end(),sigma_.begin(),
            [this](const real_t &a) -> real_t {
                return a/2.0/omega_/MU0;
            });

    } else {
        // not enough input data provided
        return EPTlibError::MissingData;

    }
    return EPTlibError::Success;
}
