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

using namespace eptlib;

// EPTHelmholtz constructor
EPTHelmholtz::
EPTHelmholtz(const size_t n0, const size_t n1, const size_t n2,
    const double d0, const double d1, const double d2,
    const double freq, const Shape &window, const int degree,
    const bool trx_phase_is_wrapped, const bool compute_variance,
    const double weight_param) :
    EPTInterface(n0,n1,n2, d0,d1,d2, freq, 1,1, trx_phase_is_wrapped),
    sg_filter_(),
    sg_window_(window),
    sg_degree_(degree),
    variance_sigma_(nullptr),
    variance_epsr_(nullptr),
    compute_variance_(compute_variance),
    weight_param_(weight_param) {
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
    // check that the sg filter is correctly set up
    if (ThereIsReferenceImage() && !std::holds_alternative<filter::AnatomicalSavitzkyGolay>(sg_filter_)) {
        sg_filter_.emplace<filter::AnatomicalSavitzkyGolay>(dd_[0],dd_[1],dd_[2], sg_window_, sg_degree_, weight_param_);
    } else
    if (!ThereIsReferenceImage() && !std::holds_alternative<filter::SavitzkyGolay>(sg_filter_)) {
        sg_filter_.emplace<filter::SavitzkyGolay>(dd_[0],dd_[1],dd_[2], sg_window_, sg_degree_);
    }
    // setup the output variables
    if (ThereIsTxSens(0)) {
        epsr_ = std::make_unique<Image<double> >(nn_[0],nn_[1],nn_[2]);
    }
    if (ThereIsTRxPhase(0,0)) {
        sigma_ = std::make_unique<Image<double> >(nn_[0],nn_[1],nn_[2]);
    }
    if (ComputeVariance()) {
        if (ThereIsSigma()) {
            variance_sigma_ = std::make_unique<Image<double> >(nn_[0],nn_[1],nn_[2]);
        }
        if (ThereIsEpsr()) {
            variance_epsr_ = std::make_unique<Image<double> >(nn_[0],nn_[1],nn_[2]);
        }
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
    // setup the input
    Image<std::complex<double> > tx_sens_c(nn_[0], nn_[1], nn_[2]);
    Image<std::complex<double> > eps_c(nn_[0], nn_[1], nn_[2]);
    std::unique_ptr<Image<std::complex<double> > > variance_c = nullptr;
    for (int idx = 0; idx<eps_c.GetNVox(); ++idx) {
        std::complex<double> exponent = std::complex<double>(0.0, 0.5 * GetTRxPhase(0,0)->At(idx));
        tx_sens_c(idx) = GetTxSens(0)->At(idx) * std::exp(exponent);
        eps_c(idx) = nancd;
    }
    if (ComputeVariance()) {
        variance_c = std::make_unique<Image<std::complex<double> > >(nn_[0], nn_[1], nn_[2]);
        variance_c->GetData().assign(eps_c.GetNVox(), nand);
    }
    // compute the laplacian
    if (!ComputeVariance()) {
        if (ThereIsReferenceImage()) {
            std::get<filter::AnatomicalSavitzkyGolay>(sg_filter_).Apply(DifferentialOperator::Laplacian, &eps_c, tx_sens_c, *reference_image_);
        } else {
            std::get<filter::SavitzkyGolay>(sg_filter_).Apply(DifferentialOperator::Laplacian, &eps_c, tx_sens_c);
        }
    } else {
        if (ThereIsReferenceImage()) {
            std::get<filter::AnatomicalSavitzkyGolay>(sg_filter_).Apply(DifferentialOperator::Laplacian, &eps_c, variance_c.get(), tx_sens_c, *reference_image_);
        } else {
            std::get<filter::SavitzkyGolay>(sg_filter_).Apply(DifferentialOperator::Laplacian, &eps_c, variance_c.get(), tx_sens_c);
        }
    }
    // extract the output
    for (int idx = 0; idx<eps_c.GetNVox(); ++idx) {
        eps_c(idx) /= -MU0*omega_*omega_*tx_sens_c(idx);
        epsr_ ->At(idx) =  std::real(eps_c(idx))/EPS0;
        sigma_->At(idx) = -std::imag(eps_c(idx))*omega_;
        if (ComputeVariance()) {
            double coeff = MU0 * omega_*omega_ * GetTxSens(0)->At(idx)*GetTxSens(0)->At(idx);
            variance_epsr_->At(idx) = (std::real(variance_c->At(idx)) * std::real(tx_sens_c(idx))*std::real(tx_sens_c(idx)) + std::imag(variance_c->At(idx)) * std::imag(tx_sens_c(idx))*std::imag(tx_sens_c(idx))) / (EPS0*EPS0 * coeff*coeff);
            variance_sigma_->At(idx) = (std::real(variance_c->At(idx)) * std::imag(tx_sens_c(idx))*std::imag(tx_sens_c(idx)) + std::imag(variance_c->At(idx)) * std::real(tx_sens_c(idx))*std::real(tx_sens_c(idx))) * omega_*omega_ / (coeff*coeff);
        }
    }
    return;
}

// EPTHelmholtz magnitude-based method
void EPTHelmholtz::
MagnitudeEPTHelm() {
    // setup the input
    epsr_->GetData().assign(epsr_->GetNVox(), nand);
    if (ComputeVariance()) {
        variance_epsr_->GetData().assign(epsr_->GetNVox(), nand);
    }
    // compute the laplacian
    if (!ComputeVariance()) {
        if (ThereIsReferenceImage()) {
            std::get<filter::AnatomicalSavitzkyGolay>(sg_filter_).Apply(DifferentialOperator::Laplacian, epsr_.get(), *GetTxSens(0), *reference_image_);
        } else {
            std::get<filter::SavitzkyGolay>(sg_filter_).Apply(DifferentialOperator::Laplacian, epsr_.get(), *GetTxSens(0));
        }
    } else {
        if (ThereIsReferenceImage()) {
            std::get<filter::AnatomicalSavitzkyGolay>(sg_filter_).Apply(DifferentialOperator::Laplacian, epsr_.get(), variance_epsr_.get(), *GetTxSens(0), *reference_image_);
        } else {
            std::get<filter::SavitzkyGolay>(sg_filter_).Apply(DifferentialOperator::Laplacian, epsr_.get(), variance_epsr_.get(), *GetTxSens(0));
        }
    }
    // extract the output
    for (int idx = 0; idx<epsr_->GetNVox(); ++idx) {
        double coeff = EPS0*MU0*omega_*omega_ * GetTxSens(0)->At(idx);
        epsr_->At(idx) /= -coeff;
        if (ComputeVariance()) {
            variance_epsr_->At(idx) /= coeff*coeff;
        }
    }
    return;
}

// EPTHelmholtz phase-based method
void EPTHelmholtz::
PhaseEPTHelm() {
    // setup the input
    sigma_->GetData().assign(sigma_->GetNVox(), nand);
    if(ComputeVariance()) {
        variance_sigma_->GetData().assign(sigma_->GetNVox(), nand);
    }
    // compute the laplacian
    if (PhaseIsWrapped() && !ComputeVariance()) {
        if (ThereIsReferenceImage()) {
            throw std::runtime_error("Feature not implemented!");
        } else {
            std::get<filter::SavitzkyGolay>(sg_filter_).ApplyWrappedPhase(DifferentialOperator::Laplacian, sigma_.get(), *GetTRxPhase(0,0));
        }
    } else if (!ComputeVariance()) {
        if (ThereIsReferenceImage()) {
            std::get<filter::AnatomicalSavitzkyGolay>(sg_filter_).Apply(DifferentialOperator::Laplacian, sigma_.get(), *GetTRxPhase(0,0), *reference_image_);
        } else {
            std::get<filter::SavitzkyGolay>(sg_filter_).Apply(DifferentialOperator::Laplacian, sigma_.get(), *GetTRxPhase(0,0));
        }
    } else {
        if (ThereIsReferenceImage()) {
            std::get<filter::AnatomicalSavitzkyGolay>(sg_filter_).Apply(DifferentialOperator::Laplacian, sigma_.get(), variance_sigma_.get(), *GetTRxPhase(0,0), *reference_image_);
        } else {
            std::get<filter::SavitzkyGolay>(sg_filter_).Apply(DifferentialOperator::Laplacian, sigma_.get(), variance_sigma_.get(), *GetTRxPhase(0,0));
        }
    }
    // extract the output
    for (int idx = 0; idx<sigma_->GetNVox(); ++idx) {
        double coeff = 2.0*MU0*omega_;
        sigma_->At(idx) /= coeff;
        if (ComputeVariance()) {
            variance_sigma_->At(idx) /= coeff*coeff;
        }
    }
    return;
}
