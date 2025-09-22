/*****************************************************************************
*
*     Program: EPTlib
*     Author: Alessandro Arduino <a.arduino@inrim.it>
*
*  MIT License
*
*  Copyright (c) 2020-2025  Alessandro Arduino
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

#include "eptlib/ept_helmholtz_chi2.h"

#include <cmath>
#include <limits>

#include "eptlib/ept_helmholtz.h"

using namespace eptlib;

// EPTHelmholtzChi2 constructor
EPTHelmholtzChi2::
EPTHelmholtzChi2(const size_t n0, const size_t n1, const size_t n2,
    const double d0, const double d1, const double d2,
    const double freq, const std::vector<Shape> &windows,
    const int degree, const bool admit_unphysical_values,
    const double weight_param) :
    EPTInterface(n0,n1,n2, d0,d1,d2, freq, 1,1, false),
    windows_(windows),
    degree_(degree),
    variance_sigma_(nullptr),
    variance_epsr_(nullptr),
    index_sigma_(nullptr),
    index_epsr_(nullptr),
    admit_unphysical_values_(admit_unphysical_values),
    weight_param_(weight_param) {
    return;
}

// EPTHelmholtzChi2 destructor
EPTHelmholtzChi2::
~EPTHelmholtzChi2() {
    return;
}

// EPTHelmholtzChi2 run
EPTlibError EPTHelmholtzChi2::
Run() {
    if (!ThereIsTRxPhase(0,0)) {
        return EPTlibError::MissingData;
    }
    // setup the output variables
    bool trx_phase_is_wrapped = false;
    bool compute_variance = true;
    size_t n_vox = Prod(nn_);
    if (ThereIsTxSens(0)) {
        epsr_          = std::make_unique<Image<double> >(nn_[0], nn_[1], nn_[2]);
        variance_epsr_ = std::make_unique<Image<double> >(nn_[0], nn_[1], nn_[2]);
        index_epsr_    = std::make_unique<Image<int>    >(nn_[0], nn_[1], nn_[2]);
        epsr_         ->GetData().assign(n_vox, nand);
        variance_epsr_->GetData().assign(n_vox, nand);
        index_epsr_   ->GetData().assign(n_vox, -1);
    }
    if (ThereIsTRxPhase(0, 0)) {
        sigma_          = std::make_unique<Image<double> >(nn_[0],nn_[1],nn_[2]);
        variance_sigma_ = std::make_unique<Image<double> >(nn_[0],nn_[1],nn_[2]);
        index_sigma_    = std::make_unique<Image<int>    >(nn_[0],nn_[1],nn_[2]);
        sigma_         ->GetData().assign(n_vox,nand);
        variance_sigma_->GetData().assign(n_vox,nand);
        index_sigma_   ->GetData().assign(n_vox,-1);
    }
    // loop over the shapes
    auto freq = omega_/2.0/PI;
    for (int idx_s = 0; idx_s<windows_.size(); ++idx_s) {
        // setup Helmholtz-EPT
        EPTHelmholtz hept(nn_[0],nn_[1],nn_[2], dd_[0],dd_[1],dd_[2], freq, windows_[idx_s], degree_, trx_phase_is_wrapped, compute_variance, weight_param_);
        if (ThereIsTxSens(0)) {
            for (int id_tx = 0; id_tx < n_tx_ch_; ++id_tx) {
                hept.SetTxSensitivity(*GetTxSens(id_tx));
            }
        }
        if (ThereIsTRxPhase(0, 0)) {
            for (int id_tx = 0; id_tx < n_tx_ch_; ++id_tx) {
                for (int id_rx = 0; id_rx < n_rx_ch_; ++id_rx) {
                    hept.SetTRxPhase(*GetTRxPhase(id_tx, id_rx));
                }
            }
        }
        if (ThereIsReferenceImage()) {
            hept.SetReferenceImage(*GetReferenceImage());
        }
        // run Helmholtz-EPT
        EPTlibError error = hept.Run();
        if (error!=EPTlibError::Success) {
            return error;
        }
        // update output with the partial results
        for (int idx = 0; idx<n_vox; ++idx) {
            if (ThereIsTxSens(0)) {
                if (!std::isnan(hept.GetRelativePermittivityVariance()->At(idx)) &&
                    (hept.GetRelativePermittivityVariance()->At(idx) < variance_epsr_->At(idx) ||
                        std::isnan(variance_epsr_->At(idx)))) {
                    bool is_physical = hept.GetRelativePermittivity()->At(idx) > 1.0;
                    if (AdmitUnphysicalValues() || is_physical) {
                        epsr_         ->At(idx) = hept.GetRelativePermittivity()->At(idx);
                        variance_epsr_->At(idx) = hept.GetRelativePermittivityVariance()->At(idx);
                        index_epsr_   ->At(idx) = idx_s;
                    }
                }
            }
            if (ThereIsTRxPhase(0, 0)) {
                if (!std::isnan(hept.GetElectricConductivityVariance()->At(idx)) &&
                    (hept.GetElectricConductivityVariance()->At(idx) < variance_sigma_->At(idx) ||
                        std::isnan(variance_sigma_->At(idx)))) {
                    bool is_physical = hept.GetElectricConductivity()->At(idx) > 0.0;
                    if (AdmitUnphysicalValues() || is_physical) {
                        sigma_         ->At(idx) = hept.GetElectricConductivity()->At(idx);
                        variance_sigma_->At(idx) = hept.GetElectricConductivityVariance()->At(idx);
                        index_sigma_   ->At(idx) = idx_s;
                    }
                }
            }
        }
    }
    return EPTlibError::Success;
}
