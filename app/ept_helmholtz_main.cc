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
#include <exception>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "eptlib/ept_helmholtz.h"
#include "eptlib/io/io_hdf5.h"
#include "eptlib/io/io_toml.h"
#include "eptlib/version.h"

#include "ept_helmholtz_main.h"

using namespace std;
using namespace eptlib;

constexpr char software::name[];

void GetH5Address(const string &address, string &fname, string &url, string &urn);

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
    vector<vector<double> > tx_sens(0);
    bool thereis_tx_sens = false;
    vector<vector<double> > trx_phase(0);
    bool thereis_trx_phase = false;
    string ref_img_address;
    std::vector<double> ref_img;
    bool thereis_ref_img;
    string sigma_address;
    bool thereis_sigma;
    string epsr_address;
    bool thereis_epsr;

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
            std::string fname;
            std::string uri;
            eptlib::io::GetAddress(tx_sens_address,fname,uri);
            size_t snip = fname.find_last_of(".");
            if (fname.substr(snip)==".raw") {
                ifstream ifile(fname,ios::binary);
                if (ifile.is_open()) {
                    vector<double> tmp(n_vox);
                    ifile.read(reinterpret_cast<char*>(tmp.data()),n_vox*sizeof(double));
                    ifile.close();
                    tx_sens.push_back(tmp);
                } else {
                    cout<<"\n"
                        <<"ERROR: impossible read file '"<<tx_sens_address<<"'"<<endl;
                    return 2;
                }
            } else {
                io::IOh5 tx_sens_file(fname, io::Mode::In);
                vector<double> tmp(n_vox);
                io::State_t io_state = tx_sens_file.ReadDataset(tmp,nn, "/",uri);
                if (io_state!=io::State::Success) {
                    cout<<"\n"
                        <<"ERROR: impossible read file '"<<tx_sens_address<<"'\n"
                        <<"       "<<ToString(io_state)<<endl;
                    return 2;
                }
                tx_sens.push_back(tmp);
            }
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
                string fname;
                string uri;
                io::GetAddress(trx_phase_address, fname,uri);
                size_t snip = fname.find_last_of(".");
                if (fname.substr(snip)==".raw") {
                    ifstream ifile(fname,ios::binary);
                    if (ifile.is_open()) {
                        std::vector<double> tmp(n_vox);
                        ifile.read(reinterpret_cast<char*>(tmp.data()),n_vox*sizeof(double));
                        ifile.close();
                        trx_phase.push_back(tmp);
                    } else {
                        cout<<"\n"
                            <<"ERROR: impossible read file '"<<trx_phase_address<<"'"<<endl;
                        return 2;
                    }
                } else {
                    io::IOh5 trx_phase_file(fname, io::Mode::In);
                    vector<double> tmp(n_vox);
                    io::State_t io_state = trx_phase_file.ReadDataset(tmp,nn, "/",uri);
                    if (io_state!=io::State::Success) {
                        cout<<"\n"
                            <<"ERROR: impossible read file '"<<trx_phase_address<<"'\n"
                            <<"       "<<ToString(io_state)<<endl;
                        return 2;
                    }
                    trx_phase.push_back(tmp);
                }
                cout<<"'"<<trx_phase_address<<"', ";
            }
        }
        cout<<"]\n";
    }
    thereis_ref_img = !io_toml.GetValue<string>(ref_img_address, "input.reference-img");
    if (thereis_ref_img) {
        string fname;
        string uri;
        ref_img.resize(n_vox);
        io::GetAddress(ref_img_address, fname,uri);
        size_t snip = fname.find_last_of(".");
        if (fname.substr(snip)==".raw") {
            ifstream ifile(fname,ios::binary);
            if (ifile.is_open()) {
                ifile.read(reinterpret_cast<char*>(ref_img.data()),n_vox*sizeof(double));
                ifile.close();
            } else {
                cout<<"\n"
                    <<"ERROR: impossible read file '"<<ref_img_address<<"'"<<endl;
                return 2;
            }
        } else {
            io::IOh5 ref_img_file(fname, io::Mode::In);
            io::State_t io_state = ref_img_file.ReadDataset(ref_img,nn, "/",uri);
            if (io_state!=io::State::Success) {
                    cout<<"\n"
                        <<"ERROR: impossible read file '"<<ref_img_address<<"'\n"
                        <<"       "<<ToString(io_state)<<endl;
                    return 2;
                }
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

    // load the files and perform the EPT
    std::array<int,NDIM> rr = {1,1,1};
    eptlib::Shape kernel_shape = eptlib::shapes::Cross(rr);
    eptlib::EPTHelmholtz ept_helm(freq, nn, dd, kernel_shape);
    if (thereis_tx_sens) {
        for (int id_tx = 0; id_tx<n_tx_ch; ++id_tx) {
            ept_helm.SetTxSensitivity(tx_sens[id_tx].data(), id_tx);
        }
    }
    if (thereis_trx_phase) {
        for (int id_rx = 0; id_rx<n_rx_ch; ++id_rx) {
            for (int id_tx = 0; id_tx<n_tx_ch; ++id_tx) {
                ept_helm.SetTRxPhase(trx_phase[id_tx+n_tx_ch*id_rx].data(), id_tx,id_rx);
            }
        }
    }
    cout<<"\nRun Helmholtz-based EPT..."<<flush;
    error = ept_helm.Run();
    if (!error) {
        cout<<"done!\n"<<endl;
    } else {
        cout<<"execution failed\n"<<endl;
        return 1;
    }

    // apply post-processing filter
    cout<<"Post-processing..."<<flush;
    rr = {2,2,2};
    kernel_shape = eptlib::shapes::Ellipsoid(rr);
    ept_helm.SetPostPro(kernel_shape);
    EPTlibError_t eptlib_error = thereis_ref_img ? ept_helm.ApplyPostPro(ref_img.data()) : ept_helm.ApplyPostPro();
    if (eptlib_error==EPTlibError::Success) {
        cout<<"done!\n"<<endl;
    } else {
        cout<<"unknown filter\n"<<endl;
    }

    // write the result
    std::vector<double> sigma(n_vox);
    std::vector<double> epsr(n_vox);
    thereis_sigma = thereis_sigma&!ept_helm.GetElectricConductivity(sigma.data());
    thereis_epsr = thereis_epsr&!ept_helm.GetRelativePermittivity(epsr.data());
    if (thereis_sigma) {
        string fname;
        string uri;
        io::GetAddress(sigma_address, fname,uri);
        size_t snip = fname.find_last_of(".");
        if (fname.substr(snip)==".raw") {
            ofstream ofile(fname,ios::binary);
            if (ofile.is_open()) {
                ofile.write(reinterpret_cast<char*>(sigma.data()), n_vox*sizeof(double));
                ofile.close();
            } else {
                cout<<"ERROR: impossible write file '"<<fname<<"'"<<endl;
                return 2;
            }
        } else {
            io::IOh5 sigma_file(fname, io::Mode::Append);
            io::State_t io_state = sigma_file.WriteDataset(sigma,nn, "/",uri);
            if (io_state!=io::State::Success) {
                cout<<"ERROR: impossible write file '"<<sigma_address<<"'\n"
                    <<"       "<<ToString(io_state)<<endl;
                return 2;
            }
        }
    }
    if (thereis_epsr) {
        string fname;
        string uri;
        io::GetAddress(epsr_address, fname,uri);
        size_t snip = fname.find_last_of(".");
        if (fname.substr(snip)==".raw") {
            ofstream ofile(fname,ios::binary);
            if (ofile.is_open()) {
                ofile.write(reinterpret_cast<char*>(epsr.data()), n_vox*sizeof(double));
                ofile.close();
            } else {
                cout<<"ERROR: impossible write file '"<<fname<<"'"<<endl;
                return 2;
            }
        } else {
            io::IOh5 epsr_file(fname, io::Mode::Append);
            io::State_t io_state = epsr_file.WriteDataset(epsr,nn, "/",uri);
            if (io_state!=io::State::Success) {
                cout<<"ERROR: impossible write file '"<<epsr_address<<"'\n"
                    <<"       "<<ToString(io_state)<<endl;
                return 2;
            }
        }
    }

    return 0;
}
