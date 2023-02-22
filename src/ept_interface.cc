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

#include "eptlib/ept_interface.h"

#include <cstring>

using namespace eptlib;

// EPTInterface constructor
EPTInterface::
EPTInterface(const size_t n0, const size_t n1, const size_t n2,
    const double d0, const double d1, const double d2,
    const double freq, const size_t n_tx_ch, const size_t n_rx_ch,
    const bool trx_phase_is_wrapped) :
    nn_({n0,n1,n2}),
    dd_({d0,d1,d2}),
    omega_(2.0*PI*freq),
    n_tx_ch_(n_tx_ch),
    n_rx_ch_(n_rx_ch),
    tx_sens_(n_tx_ch,nullptr),
    trx_phase_(n_tx_ch*n_rx_ch,nullptr),
    sigma_(nullptr),
    epsr_(nullptr),
    trx_phase_is_wrapped_(trx_phase_is_wrapped) {
    return;
}

// EPTInterface destructor
EPTInterface::
~EPTInterface() {
    return;
}

// EPTInterface set Tx sensitivity
EPTlibError EPTInterface::
SetTxSensitivity(const Image<double> &tx_sens, const size_t tx_ch) {
    if (tx_ch>=n_tx_ch_) {
        return EPTlibError::OutOfRange;
    }
    if (!CheckSizes(tx_sens.GetSize(),nn_)) {
        return EPTlibError::WrongDataFormat;
    }
    tx_sens_[tx_ch] = &tx_sens;
    return EPTlibError::Success;
}

// EPTInterface set TRx phase
EPTlibError EPTInterface::
SetTRxPhase(const Image<double> &trx_phase, const size_t tx_ch, const size_t rx_ch) {
    if (tx_ch>=n_tx_ch_ || rx_ch>=n_rx_ch_) {
        return EPTlibError::OutOfRange;
    }
    if (!CheckSizes(trx_phase.GetSize(),nn_)) {
        return EPTlibError::WrongDataFormat;
    }
    trx_phase_[tx_ch+rx_ch*n_tx_ch_] = &trx_phase;
    return EPTlibError::Success;
}
