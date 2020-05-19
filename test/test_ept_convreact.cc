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

#include "eptlib/ept_convreact.h"

#include <complex>
#include <fstream>

#include "eptlib/util.h"

using namespace eptlib;

TEST(ConvReactEPTGTest,Run) {
    constexpr int n_dim = 3;
    constexpr std::complex<double> J = std::complex<double>(0.0,1.0);
    const std::array<int,n_dim> nn = {20,20,20};
    const std::array<double,n_dim> dd = {1.0e-3,1.0e-3,1.0e-3};
    const int n_vox = std::accumulate(nn.begin(),nn.end(),1,std::multiplies<int>());
    std::array<int,n_dim> ii;
    Image<double> tx_sens(nn[0],nn[1],nn[2]);
    Image<double> trx_phase(nn[0],nn[1],nn[2]);
    double freq = 64.0e6;
    double omega = 2.0*PI*freq;
    double epsr = 60.0;
    double sigma = 0.1;
    // define a plane wave through the medium
    std::complex<double> kappa = std::sqrt(omega*omega*MU0*std::complex<double>(epsr*EPS0,-sigma/omega));
    for (int idx = 0; idx<n_vox; ++idx) {
        IdxToMultiIdx(ii,idx,nn);
        std::complex<double> b1p = 0.5*J*std::exp(-J*kappa*(ii[2]*dd[2]));
        tx_sens[idx] = std::abs(b1p);
        trx_phase[idx] = std::arg(b1p)*2.0;
    }
    // define the kernel shape
    std::array<int,n_dim> rr = {1,1,1};
    Shape shape = shapes::Cross(rr);
    ASSERT_TRUE(shape.IsSymmetric());
    // define the convection-reaction EPT engine
    EPTConvReact ept_cr(freq,nn,dd,shape);
    ept_cr.SetTxSensitivity(&tx_sens);
    ept_cr.SetTRxPhase(&trx_phase);
    ept_cr.Run();
    // collect the results
    Image<double> epsr_ept;
    Image<double> sigma_ept;
    ept_cr.GetRelativePermittivity(&epsr_ept);
    ept_cr.GetElectricConductivity(&sigma_ept);

    std::ofstream ofile("test_sigma.raw",std::ios::binary);
    ofile.write(reinterpret_cast<char*>(sigma_ept.GetData().data()), sizeof(double)*n_vox);
    ofile.close();
    ofile.open("test_epsr.raw",std::ios::binary);
    ofile.write(reinterpret_cast<char*>(epsr_ept.GetData().data()), sizeof(double)*n_vox);
    ofile.close();

    //
    for (int idx = 0; idx<n_vox; ++idx) {
        IdxToMultiIdx(ii,idx,nn);
        bool kernel_in_domain = true;
        for (int d = 0; d<n_dim; ++d) {
            if (ii[d]==0 || ii[d]==nn[d]-1) {
                kernel_in_domain = false;
                break;
            }
        }
        if (kernel_in_domain) {
            ASSERT_NEAR(epsr_ept[idx],epsr,1e-3);
            ASSERT_NEAR(sigma_ept[idx],sigma,1e-3);
        } else {
            ASSERT_EQ(epsr_ept[idx],0.0);
            ASSERT_EQ(sigma_ept[idx],0.0);
        }
    }
    return;
}
