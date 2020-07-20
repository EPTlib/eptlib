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

#include "gtest/gtest.h"

#include "eptlib/ept_interface.h"

#include <array>
#include <numeric>
#include <vector>

using namespace eptlib;

/*
 * Simple derived class to test the abstract interface.
 */
class SubEPTInterface : public EPTInterface {
    public:
        /*
         * Constructor
         */
        template <typename Tn,typename Td>
        SubEPTInterface(Tn nn,Td dd) :
            EPTInterface(64.0e6,nn,dd) {
            return;
        }
        /*
         * Fake implementation of the abstract method.
         * 
         * Copy the transmit sensitivity in sigma and the transceive phase in epsr.
         */
        EPTlibError Run() override {
            thereis_sigma_ = true;
            thereis_epsr_ = true;
            sigma_ = Image<double>(nn_[0],nn_[1],nn_[2]);
            epsr_ = Image<double>(nn_[0],nn_[1],nn_[2]);
            for (int idx = 0; idx<n_vox_; ++idx) {
                sigma_[idx] = (*tx_sens_[0])[idx];
                epsr_[idx] = (*trx_phase_[0])[idx];
            }
            return EPTlibError::Success;
        }
};

TEST(EPTInterfaceGTest,SettersGetters) {
    constexpr int n_dim = 3;
    const std::array<int,n_dim> nn = {100,100,100};
    const std::array<double,n_dim> dd = {0.002,0.002,0.002};
    const int n_vox = std::accumulate(nn.begin(),nn.end(),1,std::multiplies<int>());
    Image<double> tx_sens(nn[0],nn[1],nn[2]);
    Image<double> trx_phase(nn[0],nn[1],nn[2]);
    for (int idx = 0; idx<n_vox; ++idx) {
        tx_sens[idx] = idx+1;
        trx_phase[idx] = -idx-1;
    }
    //
    EPTlibError ierr;
    SubEPTInterface epti(nn,dd);
    //
    ierr = epti.SetTxSensitivity(&tx_sens);
    ASSERT_EQ(ierr,EPTlibError::Success);
    //
    ierr = epti.SetTRxPhase(&trx_phase);
    ASSERT_EQ(ierr,EPTlibError::Success);
    //
    ierr = epti.Run();
    ASSERT_EQ(ierr,EPTlibError::Success);
    //
    Image<double> sigma;
    Image<double> epsr;
    //
    ierr = epti.GetElectricConductivity(&sigma);
    ASSERT_EQ(ierr,EPTlibError::Success);
    real_t err_s = std::inner_product(tx_sens.GetData().begin(),tx_sens.GetData().end(),sigma.GetData().begin(),0.0,std::plus<double>(),std::minus<double>());
    ASSERT_DOUBLE_EQ(err_s,0.0);
    //
    ierr = epti.GetRelativePermittivity(&epsr);
    ASSERT_EQ(ierr,EPTlibError::Success);
    real_t err_e = std::inner_product(trx_phase.GetData().begin(),trx_phase.GetData().end(),epsr.GetData().begin(),0.0,std::plus<double>(),std::minus<double>());
    ASSERT_DOUBLE_EQ(err_e,0.0);
    //
    return;
}
