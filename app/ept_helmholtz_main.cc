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

#include <toml/toml.h>

#include "eptlib/ept_helmholtz.h"
#include "eptlib/io_hdf5.h"
#include "eptlib/version.h"

#include "ept_helmholtz_main.h"

using namespace std;
using namespace eptlib;

constexpr char software::name[];

template <typename T, typename U>
void toml_get_array_of(const toml::Value &v, const string &key, U &vec);

void GetH5Address(const string &address, string &fname, string &url, string &urn);

int main(int argc, char **argv) {
    // starting boilerplate
    cout<<eptlib::project::str()<<" ("<<eptlib::build::str()<<") ["<<eptlib::compiler::str()<<"]\n"<<endl;
    cout<<"MIT License\n"
        <<"Copyright (c) 2020  Alessandro Arduino\n"
        <<"Istituto Nazionale di Ricerca Metrologica (INRiM)\n"<<endl;
    cout<<"THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR\n"
        <<"IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,\n"
        <<"FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE\n"
        <<"AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER\n"
        <<"LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,\n"
        <<"OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE\n"
        <<"SOFTWARE.\n"<<endl;

    // check if enough input are provided
    if (argc<2) {
        cout<<"Usage example: "<<software::str()<<" <configuration file>"<<endl;
        return 1;
    }

    // load and parse the configuration file
    string fname_config(argv[1]);
    ifstream ifile_config(fname_config,ios::in);
    toml::ParseResult config_pr = toml::parse(ifile_config);
    if (!config_pr.valid()) {
        cout<<"ERROR: bad toml format. "<<config_pr.errorReason<<endl;
        return 2;
    }
    const toml::Value& config_v = config_pr.value;

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
    string sigma_address;
    string epsr_address;

    // initialise the channel wildcard variables
    char chwc_tx = '>';
    char chwc_rx = '<';
    int chwc_start_from = 0;
    int chwc_step = 1;

    // read the input wildcards
    {
        const toml::Value* config_t = config_v.find("input.wildcard");
        if (config_t && config_t->is<toml::Table>()) {
            const toml::Value* config_x;
            config_x = config_t->find("tx-character");
            if (config_x && config_x->is<string>()) {
                chwc_tx = config_x->as<string>().c_str()[0];
            }
            config_x = config_t->find("rx-character");
            if (config_x && config_x->is<string>()) {
                chwc_rx = config_x->as<string>().c_str()[0];
            }
            config_x = config_t->find("start-from");
            if (config_x && config_x->is<int>()) {
                chwc_start_from = config_x->as<int>();
            }
            config_x = config_t->find("step");
            if (config_x && config_x->is<int>()) {
                chwc_step = config_x->as<int>();
            }
        }
        if (chwc_tx==chwc_rx) {
            cout<<"ERROR: the wildcard characters for Tx and Rx channels are equal"<<endl;
            return 2;
        }
        cout<<"# Wildcard characters\n"
            <<"# Tx channel char: '"<<chwc_tx<<"'\n"
            <<"# Rx channel char: '"<<chwc_rx<<"'\n"
            <<"# Start value: "<<chwc_start_from<<"\n"
            <<"# Step: "<<chwc_step<<"\n"
            <<endl;
    }

    // read the input data
    try {
        // ...title
        title = config_v.get<string>("title");
        cout<<"  "<<title<<"\n";
        description = config_v.get<string>("description");
        cout<<"  "<<description<<"\n";
        // ...mesh
        cout<<"\n  mesh size: ";
        toml_get_array_of<int>(config_v,"mesh.size",nn);
        cout<<"["<<nn[0]<<", "<<nn[1]<<", "<<nn[2]<<"]\n";
        n_vox = Prod(nn);
        cout<<"  number of voxels: "<<n_vox<<"\n";
        cout<<"  mesh step: ";
        toml_get_array_of<double>(config_v,"mesh.step",dd);
        cout<<"["<<dd[0]<<", "<<dd[1]<<", "<<dd[2]<<"]\n";
        // ...input
        freq = config_v.get<double>("input.frequency");
        cout<<"\n  frequency: "<<freq<<"\n";
        n_tx_ch = config_v.get<int>("input.tx-channels");
        cout<<"  number of Tx channels: "<<n_tx_ch<<"\n";
        n_rx_ch = config_v.get<int>("input.rx-channels");
        cout<<"  number of Rx channels: "<<n_rx_ch<<"\n";
        {
            const toml::Value* config_x;
            config_x = config_v.find("input.tx-sensitivity");
            if (config_x && config_x->is<string>()) {
                tx_sens_address_wc = config_x->as<string>();
                thereis_tx_sens = true;
                cout<<"  Tx sensitivity wildcard: '"<<tx_sens_address_wc<<"'\n";
                cout<<"    exp'd: [";
                for (int id_tx = 0; id_tx<n_tx_ch; ++id_tx) {
                    string tx_sens_address(tx_sens_address_wc);
                    replace(tx_sens_address.begin(),tx_sens_address.end(),chwc_tx,to_string(id_tx*chwc_step+chwc_start_from).c_str()[0]);
                    size_t snip = tx_sens_address.find_last_of(".");
                    if (tx_sens_address.substr(snip) == ".raw") {
                        ifstream ifile(tx_sens_address,ios::binary);
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
                        string tx_sens_fname, tx_sens_url, tx_sens_urn;
                        GetH5Address(tx_sens_address, tx_sens_fname, tx_sens_url, tx_sens_urn);
                        io::IOh5 tx_sens_file(tx_sens_fname, io::Mode::In);
                        vector<double> tmp(n_vox);
                        io::State_t io_state = tx_sens_file.ReadDataset(tmp,nn, tx_sens_url,tx_sens_urn);
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
        }
        {
            const toml::Value* config_x;
            config_x = config_v.find("input.trx-phase");
            if (config_x && config_x->is<string>()) {
                trx_phase_address_wc = config_x->as<string>();
                thereis_trx_phase = true;
                cout<<"  TRx phase wildcard: '"<<trx_phase_address_wc<<"'\n";
                cout<<"    exp'd: [ ";
                for (int id_rx = 0; id_rx<n_rx_ch; ++id_rx) {
                    for (int id_tx = 0; id_tx<n_tx_ch; ++id_tx) {
                        string trx_phase_address(trx_phase_address_wc);
                        replace(trx_phase_address.begin(),trx_phase_address.end(),chwc_tx,to_string(id_tx*chwc_step+chwc_start_from).c_str()[0]);
                        replace(trx_phase_address.begin(),trx_phase_address.end(),chwc_rx,to_string(id_rx*chwc_step+chwc_start_from).c_str()[0]);
                        size_t snip = trx_phase_address.find_last_of(".");
                        if (trx_phase_address.substr(snip) == ".raw") {
                            ifstream ifile(trx_phase_address,ios::binary);
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
                            string trx_phase_fname, trx_phase_url, trx_phase_urn;
                            GetH5Address(trx_phase_address, trx_phase_fname, trx_phase_url, trx_phase_urn);
                            io::IOh5 trx_phase_file(trx_phase_fname, io::Mode::In);
                            vector<double> tmp(n_vox);
                            io::State_t io_state = trx_phase_file.ReadDataset(tmp,nn, trx_phase_url,trx_phase_urn);
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
        }
        {
            const toml::Value* config_x;
            config_x = config_v.find("input.reference-img");
            if (config_x && config_x->is<string>()) {
                ref_img.resize(n_vox);
                ref_img_address = config_x->as<string>();
                size_t snip = ref_img_address.find_last_of(".");
                if (ref_img_address.substr(snip) == ".raw") {
                    ifstream ifile(ref_img_address,ios::binary);
                    if (ifile.is_open()) {
                        ifile.read(reinterpret_cast<char*>(ref_img.data()),n_vox*sizeof(double));
                        ifile.close();
                    } else {
                        cout<<"\n"
                            <<"ERROR: impossible read file '"<<ref_img_address<<"'"<<endl;
                        return 2;
                    }
                } else {
                    string ref_img_fname, ref_img_url, ref_img_urn;
                    GetH5Address(ref_img_address, ref_img_fname, ref_img_url, ref_img_urn);
                    io::IOh5 ref_img_file(ref_img_fname, io::Mode::In);
                    io::State_t io_state = ref_img_file.ReadDataset(ref_img,nn, ref_img_url,ref_img_urn);
                    if (io_state!=io::State::Success) {
                        cout<<"\n"
                            <<"ERROR: impossible read file '"<<ref_img_address<<"'\n"
                            <<"       "<<ToString(io_state)<<endl;
                        return 2;
                    }
                }
                cout<<"  Reference image: '"<<ref_img_address<<"'\n";
            } else {
                ref_img_address = "";
                ref_img.resize(0);
            }
        }
        // ...output
        sigma_address = config_v.get<string>("output.electric-conductivity");
        cout<<"\n  Output electric conductivity: '"<<sigma_address<<"'\n";
        epsr_address = config_v.get<string>("output.relative-permittivity");
        cout<<"  Output relative permittivity: '"<<epsr_address<<"'\n";
    } catch (runtime_error e) {
        cout<<"\n"
            <<"ERROR: "<<e.what()<<endl;
        return 2;
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
    ept_helm.Run();
    cout<<"done!\n"<<endl;

    // apply post-processing filter
    cout<<"Post-processing..."<<flush;
    rr = {2,2,2};
    kernel_shape = eptlib::shapes::Ellipsoid(rr);
    ept_helm.SetPostPro(kernel_shape);
    EPTlibError_t eptlib_error = ref_img.size()>0 ? ept_helm.ApplyPostPro(ref_img.data()) : ept_helm.ApplyPostPro();
    if (eptlib_error==EPTlibError::Success) {
        cout<<"done!\n"<<endl;
    } else {
        cout<<"unknown filter\n"<<endl;
    }

    // write the result
    std::vector<double> sigma(n_vox);
    std::vector<double> epsr(n_vox);
    ept_helm.GetElectricConductivity(sigma.data());
    ept_helm.GetRelativePermittivity(epsr.data());
    {
        size_t snip = sigma_address.find_last_of(".");
        if (sigma_address.substr(snip) == ".raw") {
            ofstream ofile(sigma_address,ios::binary);
            if (ofile.is_open()) {
                ofile.write(reinterpret_cast<char*>(sigma.data()), n_vox*sizeof(double));
                ofile.close();
            } else {
                cout<<"ERROR: impossible write file '"<<sigma_address<<"'"<<endl;
                return 2;
            }
        } else {
            string sigma_fname, sigma_url, sigma_urn;
            GetH5Address(sigma_address, sigma_fname,sigma_url,sigma_urn);
            io::IOh5 sigma_file(sigma_fname, io::Mode::Append);
            io::State_t io_state = sigma_file.WriteDataset(sigma,nn, sigma_url,sigma_urn);
            if (io_state!=io::State::Success) {
                cout<<"ERROR: impossible write file '"<<sigma_address<<"'\n"
                    <<"       "<<ToString(io_state)<<endl;
                return 2;
            }
        }
    }
    {
        size_t snip = epsr_address.find_last_of(".");
        if (epsr_address.substr(snip) == ".raw") {
            ofstream ofile(epsr_address,ios::binary);
            if (ofile.is_open()) {
                ofile.write(reinterpret_cast<char*>(epsr.data()), n_vox*sizeof(double));
                ofile.close();
            } else {
                cout<<"ERROR: impossible write file '"<<epsr_address<<"'"<<endl;
                return 2;
            }
        } else {
            string epsr_fname, epsr_url, epsr_urn;
            GetH5Address(epsr_address, epsr_fname,epsr_url,epsr_urn);
            io::IOh5 epsr_file(epsr_fname, io::Mode::Append);
            io::State_t io_state = epsr_file.WriteDataset(epsr,nn, epsr_url,epsr_urn);
            if (io_state!=io::State::Success) {
                cout<<"ERROR: impossible write file '"<<epsr_address<<"'\n"
                    <<"       "<<ToString(io_state)<<endl;
                return 2;
            }
        }
    }

    return 0;
}

template <typename T, typename U>
void toml_get_array_of(const toml::Value &v, const string &key, U &vec) {
    const toml::Array &a = v.get<toml::Array>(key);
    for (int d = 0; d<NDIM; ++d) {
        vec[d] = a[d].as<T>();
    }
    return;
};

void GetH5Address(const std::string &address, std::string &fname, std::string &url, std::string &urn) {
    size_t snip = 0;
    size_t snap = address.find_first_of(":", 0);
    fname = address.substr(snip,snap);
    snip = ++snap;
    snap = address.find_last_of("/", std::string::npos);
    url = address.substr(snip,snap-snip);
    snip = ++snap;
    urn = address.substr(snip);
    return;
};
