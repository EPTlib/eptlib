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

#include "eptlib/ept_helmholtz.h"

#include "eptlib/util.h"

using namespace eptlib;

TEST(HelmholtzEPTGTest,Run) {
    constexpr int n_dim = 3;
    constexpr std::complex<real_t> J = std::complex<real_t>(0.0,1.0);
    const std::array<int,n_dim> nn = {20,20,20};
    const std::array<real_t,n_dim> dd = {1.0e-3,1.0e-3,1.0e-3};
    const int n_vox = std::accumulate(nn.begin(),nn.end(),1,std::multiplies<int>());
    std::array<int,n_dim> ii;
    std::vector<real_t> tx_sens(n_vox);
    std::vector<real_t> trx_phase(n_vox);
    real_t freq = 64.0e6;
    real_t omega = 2.0*PI*freq;
    real_t epsr = 60.0;
    real_t sigma = 0.1;
    // define a plane wave through the medium
    std::complex<real_t> kappa = std::sqrt(omega*omega*MU0*std::complex<real_t>(epsr*EPS0,-sigma/omega));
    for (int idx = 0; idx<n_vox; ++idx) {
        IdxToMultiIdx(ii,idx,nn);
        std::complex<real_t> b1p = 0.5*J*std::exp(-J*kappa*(ii[2]*dd[2]));
        tx_sens[idx] = std::abs(b1p);
        trx_phase[idx] = std::arg(b1p)*2.0;
    }
    // define the kernel shape
    std::array<int,n_dim> rr = {1,1,1};
    Shape shape = shapes::Cross(rr);
    ASSERT_TRUE(shape.IsSymmetric());
    // define the Helmholtz EPT engine
    EPTHelmholtz ept_helm(freq,nn,dd,shape);
    ept_helm.SetTxSensitivity(tx_sens.data());
    ept_helm.SetTRxPhase(trx_phase.data());
    ept_helm.Run();
    // collect the results
    std::vector<real_t> epsr_ept(n_vox);
    std::vector<real_t> sigma_ept(n_vox);
    ept_helm.GetRelativePermittivity(epsr_ept.data());
    ept_helm.GetElectricConductivity(sigma_ept.data());
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
