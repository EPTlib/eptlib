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

// EPTInterface constructor
EPTInterface::
EPTInterface(const double freq, const std::array<int,NDIM> &nn,
    const std::array<double,NDIM> &dd, const int tx_ch, const int rx_ch) :
    omega_(2.0*PI*freq), tx_ch_(tx_ch), rx_ch_(rx_ch), nn_(nn), dd_(dd), 
    n_vox_(std::accumulate(nn.begin(),nn.end(),1,std::multiplies<int>())),
    tx_sens_(tx_ch,nullptr), trx_phase_(tx_ch*rx_ch,nullptr),
    thereis_tx_sens_(tx_ch,false), thereis_trx_phase_(tx_ch*rx_ch,false),
    sigma_(0), epsr_(0), thereis_sigma_(false), thereis_epsr_(false),
    postpro_(nullptr), thereis_postpro_(false) {
    return;
}

// EPTInterface destructor
EPTInterface::
~EPTInterface() {
    return;
}

// EPTInterface setters
EPTlibError_t EPTInterface::
SetTxSensitivity(const double *tx_sens, const int j) {
    if (j<0 || j>=tx_ch_) {
        return EPTlibError::OutOfRange;
    }
    tx_sens_[j] = tx_sens;
    thereis_tx_sens_.set(j);
    return EPTlibError::Success;
}
EPTlibError_t EPTInterface::
SetTRxPhase(const double *trx_phase, const int j, const int k) {
    if (j<0 || j>=tx_ch_ || k<0 || k>=rx_ch_) {
        return EPTlibError::OutOfRange;
    }
    trx_phase_[j+k*tx_ch_] = trx_phase;
    thereis_trx_phase_.set(j+k*tx_ch_);
    return EPTlibError::Success;
}

// EPTInterface getters
EPTlibError_t EPTInterface::
GetElectricConductivity(double *sigma) {
    if (!thereis_sigma_) {
        return EPTlibError::MissingData;
    }
    std::memcpy(sigma,sigma_.data(),n_vox_*sizeof(double));
    return EPTlibError::Success;    
}
EPTlibError_t EPTInterface::
GetRelativePermittivity(double *epsr) {
    if (!thereis_epsr_) {
        return EPTlibError::MissingData;
    }
    std::memcpy(epsr,epsr_.data(),n_vox_*sizeof(double));
    return EPTlibError::Success;
}

// EPTInterface post-processing
EPTlibError_t EPTInterface::
SetPostPro(const Shape &shape) {
    if (thereis_postpro_) {
        delete postpro_;
    }
    postpro_ = new MedianFilter(shape);
    thereis_postpro_ = true;
    return EPTlibError::Success;
}
EPTlibError_t EPTInterface::
UnsetPostPro() {
    if (thereis_postpro_) {
        delete postpro_;
        postpro_ = nullptr;
    }
    thereis_postpro_ = false;
    return EPTlibError::Success;
}
EPTlibError_t EPTInterface::
ApplyPostPro(const double *img) {
    if (!thereis_postpro_ || !(thereis_sigma_ || thereis_epsr_)) {
        return EPTlibError::MissingData;
    }
    std::vector<double> tmp(n_vox_);
    if (thereis_sigma_) {
        postpro_->ApplyFilter(tmp.data(),sigma_.data(),nn_,img);
        std::memcpy(sigma_.data(),tmp.data(),n_vox_*sizeof(double));
    }
    if (thereis_epsr_) {
        postpro_->ApplyFilter(tmp.data(),epsr_.data(),nn_,img);
        std::memcpy(epsr_.data(),tmp.data(),n_vox_*sizeof(double));
    }
    return EPTlibError::Success;
}