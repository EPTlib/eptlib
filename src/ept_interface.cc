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

#include "eptlib/ept_interface.h"

#include <cstring>

using namespace eptlib;

// EPTInterface destructor
EPTInterface::
~EPTInterface() {
    return;
}

// EPTInterface setters
EPTlibError_t EPTInterface::
SetTxSensitivity(const real_t *tx_sens, const int j) {
    if (j<0 || j>=tx_ch_) {
        return EPTlibError::OutOfRange;
    }
    tx_sens_[j] = tx_sens;
    thereis_tx_sens_.set(j);
    return EPTlibError::Success;
}
EPTlibError_t EPTInterface::
SetTRxPhase(const real_t *trx_phase, const int j, const int k) {
    if (j<0 || j>=tx_ch_ || k<0 || k>=rx_ch_) {
        return EPTlibError::OutOfRange;
    }
    trx_phase_[j+k*tx_ch_] = trx_phase;
    thereis_trx_phase_.set(j+k*tx_ch_);
    return EPTlibError::Success;
}

// EPTInterface getters
EPTlibError_t EPTInterface::
GetElectricConductivity(real_t *sigma) {
    if (!thereis_sigma_) {
        return EPTlibError::MissingData;
    }
    std::memcpy(sigma,sigma_.data(),n_vox_*sizeof(real_t));
    return EPTlibError::Success;    
}
EPTlibError_t EPTInterface::
GetRelativePermittivity(real_t *epsr) {
    if (!thereis_epsr_) {
        return EPTlibError::MissingData;
    }
    std::memcpy(epsr,epsr_.data(),n_vox_*sizeof(real_t));
    return EPTlibError::Success;
}
