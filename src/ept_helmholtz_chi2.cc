/*****************************************************************************
*
*     Program: EPTlib
*     Author: Alessandro Arduino <a.arduino@inrim.it>
*
*  MIT License
*
*  Copyright (c) 2020-2022  Alessandro Arduino
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

using namespace eptlib;

namespace {
    static double nand = std::numeric_limits<double>::quiet_NaN();
}

// EPTHelmholtzChi2 constructor
EPTHelmholtzChi2::
EPTHelmholtzChi2(const double freq,const std::array<int,NDIM> &nn,
    const std::array<double,NDIM> &dd,const std::vector<Shape> &shapes,
    const int degree) :
    EPTInterface(freq,nn,dd), freq_(freq), shapes_(shapes), degree_(degree),
    unphysical_values_(false) {
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
    // allocate
    if (thereis_trx_phase_.all()) {
        thereis_sigma_ = true;
        sigma_ = Image<double>(nn_[0],nn_[1],nn_[2]);
    } else {
        return EPTlibError::MissingData;
    }
    var_ = Image<double>(nn_[0],nn_[1],nn_[2]);
    shape_index_ = Image<int>(nn_[0],nn_[1],nn_[2]);
    // initialise
    for (int idx = 0; idx<n_vox_; ++idx) {
        var_[idx] = nand;
        shape_index_[idx] = -1;
        sigma_[idx] = nand;
    }
    // loop over the shapes
    for (int idx_s = 0; idx_s<shapes_.size(); ++idx_s) {
        // initialise the Helmholtz-EPT method
        EPTHelmholtz hept(freq_,nn_,dd_,shapes_[idx_s],degree_);
        // set-up Helmholtz-EPT
        hept.SetTRxPhase(trx_phase_[0]);
        hept.ToggleGetVar();
        // run Helmholtz-EPT
        EPTlibError error = hept.Run();
        if (error!=EPTlibError::Success) {
            return error;
        }
        // get the partial results
        for (int idx = 0; idx<n_vox_; ++idx) {
            if (hept.var_[idx]==hept.var_[idx] && (hept.var_[idx]<var_[idx] || !(var_[idx]==var_[idx]))) {
                bool is_physical = hept.sigma_[idx]>0.0;
                if (unphysical_values_||is_physical) {
                    var_[idx] = hept.var_[idx];
                    shape_index_[idx] = idx_s;
                    sigma_[idx] = hept.sigma_[idx];
                }
            }
        }
    }
    return Success;
}

// Get the result quality index
EPTlibError EPTHelmholtzChi2::
GetVar(Image<double> *var) {
    if (!thereis_sigma_) {
        return EPTlibError::MissingData;
    }
    *var = var_;
    return EPTlibError::Success;
}

// Get the result quality index
EPTlibError EPTHelmholtzChi2::
GetShapeIndex(Image<int> *shape_index) {
    if (!thereis_sigma_) {
        return EPTlibError::MissingData;
    }
    *shape_index = shape_index_;
    return EPTlibError::Success;
}

bool EPTHelmholtzChi2::
ToggleUnphysicalValues() {
	unphysical_values_ = !unphysical_values_;
	return unphysical_values_;
}
