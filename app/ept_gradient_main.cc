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

#include <array>
#include <algorithm>
#include <chrono>
#include <exception>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "eptlib/ept_gradient.h"
#include "eptlib/io/io_hdf5.h"
#include "eptlib/io/io_toml.h"
#include "eptlib/version.h"

#include "ept_gradient_main.h"

using namespace std;
using namespace eptlib;

constexpr char software::name[];

template <typename T>
EPTlibError_t ReadData(Image<T> &dst, const std::string &address);
template <typename T>
EPTlibError_t WriteData(const Image<T> &src, const std::array<int,NDIM> &nn, const std::string &address);

int main(int argc, char **argv) {
    // starting boilerplate
    cout<<software::str()<<endl;
    cout<<eptlib::project::str()<<" ("<<eptlib::build::str()<<") ["<<eptlib::compiler::str()<<"]\n"<<endl;
    cout<<LicenseBoilerplate()<<endl;

    // check if enough input are provided
    if (argc<2) {
        cout<<"Usage example: "<<software::str()<<" <configuration file>"<<endl;
        return 1;
    }

    // load the configuration file
    io::IOtoml io_toml(std::string(argv[1]), io::Mode::In);
    EPTlibError_t error;

    // declare the input variables
    string title;
    string description;
    array<int,NDIM> nn;
    int n_vox;
    array<double,NDIM> dd;
    double freq;
    int n_tx_ch;
    int n_rx_ch;
    string tx_sens_address_wc;
    string trx_phase_address_wc;
    vector<Image<double> > tx_sens(0);
    bool thereis_tx_sens = false;
    vector<Image<double> > trx_phase(0);
    bool thereis_trx_phase = false;
    string ref_img_address;
    Image<double> ref_img;
    bool thereis_ref_img;
    string sigma_address;
    bool thereis_sigma;
    string epsr_address;
    bool thereis_epsr;
    bool is_2d;
    double lambda;
    bool use_lcurve;
    bool thereis_lambda;
    int max_it;
    bool thereis_max_it;
    double tol;
    bool thereis_tol;

    // initialise the channel wildcard variables
    char chwc_tx = '>';
    char chwc_rx = '<';
    int chwc_start_from = 0;
    int chwc_step = 1;

    // read the input wildcards
    io_toml.GetValue<char>(chwc_tx, "input.wildcard.tx-character");
    io_toml.GetValue<char>(chwc_rx, "input.wildcard.rx-character");
    io_toml.GetValue<int>(chwc_start_from, "input.wildcard.start-from");
    io_toml.GetValue<int>(chwc_step, "input.wildcard.step");
    if (chwc_tx==chwc_rx) {
        cout<<"Fatal ERROR: the wildcard characters for Tx and Rx channels are equal"<<endl;
        return 2;
    }
    cout<<"# Wildcard characters\n"
        <<"# Tx channel char: '"<<chwc_tx<<"'\n"
        <<"# Rx channel char: '"<<chwc_rx<<"'\n"
        <<"# Start value: "<<chwc_start_from<<"\n"
        <<"# Step: "<<chwc_step<<"\n"
        <<endl;

    // read the input data
    // ...title
    error = io_toml.GetValue<string>(title,"title");
    assert(!error);
    cout<<"  "<<title<<"\n";
    error = io_toml.GetValue<string>(description,"description");
    assert(!error);
    cout<<"  "<<description<<"\n";
    // ...mesh
    error = io_toml.GetArrayOf<int>(nn,"mesh.size");
    assert(!error);
    cout<<"\n  mesh size: "
        <<"["<<nn[0]<<", "<<nn[1]<<", "<<nn[2]<<"]\n";
    n_vox = Prod(nn);
    cout<<"  number of voxels: "<<n_vox<<"\n";
    error = io_toml.GetArrayOf<double>(dd,"mesh.step");
    assert(!error);
    cout<<"  mesh step: "
        <<"["<<dd[0]<<", "<<dd[1]<<", "<<dd[2]<<"]\n";
    // ...input
    error = io_toml.GetValue<double>(freq,"input.frequency");
    assert(!error);
    cout<<"\n  frequency: "<<freq<<"\n";
    error = io_toml.GetValue<int>(n_tx_ch,"input.tx-channels");
    assert(!error);
    cout<<"  number of Tx channels: "<<n_tx_ch<<"\n";
    error = io_toml.GetValue<int>(n_rx_ch,"input.rx-channels");
    assert(!error);
    cout<<"  number of Rx channels: "<<n_rx_ch<<"\n";
    thereis_tx_sens = !io_toml.GetValue<string>(tx_sens_address_wc,"input.tx-sensitivity");
    if (thereis_tx_sens) {
        cout<<"  Tx sensitivity wildcard: '"<<tx_sens_address_wc<<"'\n";
        cout<<"    exp'd: [";
        for (int id_tx = 0; id_tx<n_tx_ch; ++id_tx) {
            string tx_sens_address(tx_sens_address_wc);
            replace(tx_sens_address.begin(),tx_sens_address.end(),chwc_tx,to_string(id_tx*chwc_step+chwc_start_from).c_str()[0]);
            Image<double> tmp(nn[0],nn[1],nn[2]);
            if (ReadData<double>(tmp, tx_sens_address)) {
                cout<<"Fatal ERROR: impossible read file '"<<tx_sens_address<<"'\n";
                return 2;
            }
            tx_sens.push_back(tmp);
            cout<<"'"<<tx_sens_address<<"', ";
        }
        cout<<"]\n";
    }
    thereis_trx_phase = !io_toml.GetValue<string>(trx_phase_address_wc,"input.trx-phase");
    if (thereis_trx_phase) {
        cout<<"  TRx phase wildcard: '"<<trx_phase_address_wc<<"'\n";
        cout<<"    exp'd: [ ";
        for (int id_rx = 0; id_rx<n_rx_ch; ++id_rx) {
            for (int id_tx = 0; id_tx<n_tx_ch; ++id_tx) {
                string trx_phase_address(trx_phase_address_wc);
                replace(trx_phase_address.begin(),trx_phase_address.end(),chwc_tx,to_string(id_tx*chwc_step+chwc_start_from).c_str()[0]);
                replace(trx_phase_address.begin(),trx_phase_address.end(),chwc_rx,to_string(id_rx*chwc_step+chwc_start_from).c_str()[0]);
                Image<double> tmp(nn[0],nn[1],nn[2]);
                if (ReadData<double>(tmp, trx_phase_address)) {
                    cout<<"Fatal ERROR: impossible read file '"<<trx_phase_address<<"'\n";
                    return 2;
                }
                trx_phase.push_back(tmp);
                cout<<"'"<<trx_phase_address<<"', ";
            }
        }
        cout<<"]\n";
    }
    thereis_ref_img = !io_toml.GetValue<string>(ref_img_address, "input.reference-img");
    if (thereis_ref_img) {
        ref_img = Image<double>(nn[0],nn[1],nn[2]);
        if (ReadData<double>(ref_img, ref_img_address)) {
            cout<<"Fatal ERROR: impossible read file '"<<ref_img_address<<"'\n";
            return 2;
        }
        cout<<"  Reference image: '"<<ref_img_address<<"'\n";
    }
    // ...output
    thereis_sigma = !io_toml.GetValue<string>(sigma_address, "output.electric-conductivity"); 
    if (thereis_sigma) {
        cout<<"\n  Output electric conductivity: '"<<sigma_address<<"'\n";
    }
    thereis_epsr = !io_toml.GetValue<string>(epsr_address, "output.relative-permittivity");
    if (thereis_epsr) {
        cout<<"  Output relative permittivity: '"<<epsr_address<<"'\n";
    }
    // ...parameters
    cout<<"\n";
    error = io_toml.GetValue<bool>(is_2d,"parameter.2d-assumption");
    if (!error) {
        cout<<"  2D assumption: '"<<(is_2d?"True":"False")<<"'\n";
    } else {
        is_2d = false;
    }
    error = io_toml.GetValue<bool>(use_lcurve,"parameter.use-lcurve");
    if (!error) {
        cout<<"  Use L-curve: '"<<(use_lcurve?"True":"False")<<"'\n";
    } else {
        use_lcurve = false;
    }
    error = io_toml.GetValue<double>(lambda,"parameter.estimate-weight");
    if (!error) {
        thereis_lambda = true;
        cout<<"  Estimate weight: "<<lambda<<"\n";
    } else {
        thereis_lambda = false;
    }
    error = io_toml.GetValue<int>(max_it,"parameter.maximum-iterations");
    if (!error) {
        thereis_max_it = true;
        cout<<"  Maximum iterations: "<<max_it<<"\n";
    } else {
        thereis_max_it = false;
    }
    error = io_toml.GetValue<double>(tol,"parameter.tolerance");
    if (!error) {
        thereis_tol = true;
        cout<<"  Tolerance: "<<tol<<"\n";
    } else {
        thereis_tol = false;
    }
    // load the files and perform the EPT
    std::array<int,NDIM> rr = {1,1,1};
    eptlib::Shape kernel_shape = eptlib::shapes::Cross(rr);
//    eptlib::Shape kernel_shape = eptlib::shapes::Ellipsoid(rr);
    eptlib::EPTGradient ept_grad(freq, nn, dd, n_tx_ch, kernel_shape, is_2d);
    if (thereis_tx_sens) {
        for (int id_tx = 0; id_tx<n_tx_ch; ++id_tx) {
            error = ept_grad.SetTxSensitivity(&(tx_sens[id_tx]), id_tx);
            if (error) {
                cout<<"\nFatal ERROR: Unable to set Tx sensitivity ch"<<id_tx<<" - "<<ToString(error)<<endl;
                return 2;
            }
        }
    }
    if (thereis_trx_phase) {
        for (int id_rx = 0; id_rx<n_rx_ch; ++id_rx) {
            for (int id_tx = 0; id_tx<n_tx_ch; ++id_tx) {
                error = ept_grad.SetTRxPhase(&(trx_phase[id_tx+n_tx_ch*id_rx]), id_tx,id_rx);
                if (error) {
                    cout<<"\nFatal ERROR: Unable to set TRx phase ch"<<id_tx<<","<<id_rx<<" - "<<ToString(error)<<endl;
                    return 2;
                }
            }
        }
    }
    if (thereis_lambda) {
        ept_grad.SetLambda(lambda);
    }

    SeedPoint sp;
    sp.ijk = {70,120,5};
    sp.epsr = 43.8;
    sp.sigma = 0.412;
    ept_grad.AddSeedPoint(sp);
    ept_grad.ToggleSeedPoints();

    ept_grad.SetRunMode(EPTGradientRun::FULL);

    cout<<"\nRun gradient EPT..."<<flush;
    auto start = std::chrono::system_clock::now();
    error = ept_grad.Run();
    auto end = std::chrono::system_clock::now();
    if (!error) {
        auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
        cout<<"done! ["<<elapsed.count()<<" s]\n"<<endl;
    } else {
        cout<<"execution failed\n"<<endl;
        return 1;
    }

    // write the result
    Image<double> sigma(nn[0],nn[1],nn[2]);
    Image<double> epsr(nn[0],nn[1],nn[2]);
    thereis_sigma = thereis_sigma&!ept_grad.GetElectricConductivity(&sigma);
    thereis_epsr = thereis_epsr&!ept_grad.GetRelativePermittivity(&epsr);
    if (thereis_sigma) {
        error = WriteData(sigma, nn, sigma_address);
        if (error!=EPTlibError::Success) {
            cout<<"Fatal ERROR: impossible write file '"<<sigma_address<<"'"<<endl;
            return 2;
        }
    }
    if (thereis_epsr) {
        error = WriteData(epsr, nn, epsr_address);
        if (error!=EPTlibError::Success) {
            cout<<"Fatal ERROR: impossible write file '"<<epsr_address<<"'"<<endl;
            return 2;
        }
    }
    return 0;
}

template <typename T>
EPTlibError_t ReadData(Image<T> &dst, const std::string &address) {
    std::string fname;
    std::string uri;
    io::GetAddress(address, fname,uri);
    size_t snip = fname.find_last_of(".");
    if (fname.substr(snip)==".raw") {
        std::ifstream ifile(fname,ios::binary);
        if (ifile.is_open()) {
            ifile.read(reinterpret_cast<char*>(dst.GetData().data()),dst.GetNVox()*sizeof(T));
            ifile.close();
        } else {
            return EPTlibError::MissingData;
        }
    } else {
        io::IOh5 ifile(fname, io::Mode::In);
        io::State_t io_state = ifile.ReadDataset(&dst, "/",uri);
        if (io_state!=io::State::Success) {
            return EPTlibError::MissingData;
        }
    }
    return EPTlibError::Success;
}
template <typename T>
EPTlibError_t WriteData(const Image<T> &src, const std::array<int,NDIM> &nn, const std::string &address) {
    std::string fname;
    std::string uri;
    io::GetAddress(address, fname,uri);
    size_t snip = fname.find_last_of(".");
    if (fname.substr(snip)==".raw") {
        std::ofstream ofile(fname,ios::binary);
        if (ofile.is_open()) {
            ofile.write(reinterpret_cast<const char*>(src.GetData().data()), src.GetNVox()*sizeof(T));
            ofile.close();
        } else {
            return EPTlibError::Unknown;
        }
    } else {
        io::IOh5 ofile(fname, io::Mode::Append);
        io::State_t io_state = ofile.WriteDataset(src, "/",uri);
        if (io_state!=io::State::Success) {
            return EPTlibError::Unknown;
        }
    }
    return EPTlibError::Success;
}
