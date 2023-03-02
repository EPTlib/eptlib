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

#include "gtest/gtest.h"

#include "eptlib/ept_interface.h"

#include <algorithm>
#include <array>
#include <memory>
#include <vector>

/*
 * Simple derived class to test the abstract interface.
 */
class FakeEPTInterface : public eptlib::EPTInterface {
    public:
        /*
         * Constructor
         */
        FakeEPTInterface(const size_t n0, const size_t n1, const size_t n2,
            const double d0, const double d1, const double d2) :
            EPTInterface(n0,n1,n2, d0,d1,d2, 64e6) {
            return;
        }
        /*
         * Fake implementation of the abstract method.
         * 
         * Copy the transmit sensitivity in sigma and the transceive phase in epsr.
         */
        eptlib::EPTlibError Run() override {
            using image_t = eptlib::Image<double>;
            sigma_ = std::make_unique<image_t>(nn_[0],nn_[1],nn_[2]);
            epsr_ = std::make_unique<image_t>(nn_[0],nn_[1],nn_[2]);
            std::copy(tx_sens_[0]->GetData().begin(),tx_sens_[0]->GetData().end(), sigma_->GetData().begin());
            std::copy(trx_phase_[0]->GetData().begin(),trx_phase_[0]->GetData().end(), epsr_->GetData().begin());
            return eptlib::EPTlibError::Success;
        }
};

TEST(EPTInterfaceGTest,SettersGetters) {
    const size_t n0 = 100;
    const size_t n1 = 100;
    const size_t n2 = 100;
    const double d0 = 0.002;
    const double d1 = 0.002;
    const double d2 = 0.002;
    eptlib::Image<double> dummy(1,1,1);
    eptlib::Image<double> tx_sens(n0,n1,n2);
    eptlib::Image<double> trx_phase(n0,n1,n2);
    std::iota(tx_sens.GetData().begin(),tx_sens.GetData().end(),1);
    std::transform(tx_sens.GetData().begin(),tx_sens.GetData().end(),trx_phase.GetData().begin(),
        [](const double a) -> double{
            return -a;
        });
    //
    FakeEPTInterface ept(n0,n1,n2, d0,d1,d2);
    ASSERT_FALSE(ept.PhaseIsWrapped());
    ASSERT_TRUE(ept.TogglePhaseIsWrapped());
    ASSERT_TRUE(ept.PhaseIsWrapped());
    ASSERT_FALSE(ept.TogglePhaseIsWrapped());
    //
    ASSERT_FALSE(ept.ThereIsTxSens(0));
    ASSERT_FALSE(ept.ThereAreAllTxSens());
    ASSERT_EQ(ept.SetTxSensitivity(tx_sens,1), eptlib::EPTlibError::OutOfRange);
    ASSERT_EQ(ept.SetTxSensitivity(dummy,0), eptlib::EPTlibError::WrongDataFormat);
    ASSERT_EQ(ept.SetTxSensitivity(tx_sens,0), eptlib::EPTlibError::Success);
    ASSERT_TRUE(ept.ThereIsTxSens(0));
    ASSERT_TRUE(ept.ThereAreAllTxSens());
    //
    ASSERT_FALSE(ept.ThereIsTRxPhase(0,0));
    ASSERT_FALSE(ept.ThereAreAllTRxPhase());
    ASSERT_EQ(ept.SetTRxPhase(trx_phase,1,0), eptlib::EPTlibError::OutOfRange);
    ASSERT_EQ(ept.SetTRxPhase(trx_phase,0,1), eptlib::EPTlibError::OutOfRange);
    ASSERT_EQ(ept.SetTRxPhase(dummy,0,0), eptlib::EPTlibError::WrongDataFormat);
    ASSERT_EQ(ept.SetTRxPhase(trx_phase,0,0), eptlib::EPTlibError::Success);
    ASSERT_TRUE(ept.ThereIsTRxPhase(0,0));
    ASSERT_TRUE(ept.ThereAreAllTRxPhase());
    //
    ASSERT_FALSE(ept.ThereIsSigma());
    ASSERT_FALSE(ept.ThereIsEpsr());
    ASSERT_EQ(ept.Run(), eptlib::EPTlibError::Success);
    ASSERT_TRUE(ept.ThereIsSigma());
    ASSERT_TRUE(ept.ThereIsEpsr());
    //
    auto tx_sens_ptr = ept.GetTxSens(0);
    auto trx_phase_ptr = ept.GetTRxPhase(0,0);
    ASSERT_EQ(tx_sens_ptr, &tx_sens);
    ASSERT_EQ(trx_phase_ptr, &trx_phase);
    //
    auto& sigma = ept.GetElectricConductivity();
    auto& epsr = ept.GetRelativePermittivity();
    for (int idx = 0; idx<tx_sens.GetNVox(); ++idx) {
        ASSERT_EQ((*sigma)(idx), tx_sens(idx));
        ASSERT_EQ((*epsr)(idx), trx_phase(idx));
    }
    return;
}
