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

#include <algorithm>
#include <exception>
#include <fstream>
#include <iostream>
#include <string>

#include <toml/toml.h>

#include "eptlib/ept_helmholtz.h"
#include "eptlib/version.h"

#include "ept_helmholtz_main.h"

using namespace std;
using namespace eptlib_helm;

constexpr char software::name[];

template <typename T, typename U>
void toml_get_array_of(const toml::Value &v, const string &key, U &vec);

int main(int argc, char **argv) {
    // starting boilerplate
    cout<<eptlib::project::str()<<" ("<<eptlib::build::str()<<") ["<<eptlib::compiler::str()<<"]\n"<<std::endl;
    cout<<"MIT License\n"
        <<"Copyright (c) 2020  Alessandro Arduino\n"
        <<"Istituto Nazionale di Ricerca Metrologica (INRiM)\n"<<std::endl;
    cout<<"THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR\n"
        <<"IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,\n"
        <<"FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE\n"
        <<"AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER\n"
        <<"LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,\n"
        <<"OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE\n"
        <<"SOFTWARE.\n"<<std::endl;

    // check if enough input are provided
    if (argc<2) {
        cout<<"Usage: "<<software::str()<<" <configuration file>"<<std::endl;
        return 1;
    }

    // load and parse the configuration file
    ifstream config_ifs(argv[1]);
    toml::ParseResult config_pr = toml::parse(config_ifs);
    if (!config_pr.valid()) {
        cout<<"ERROR: bad toml format. "<<config_pr.errorReason<<std::endl;
        return 2;
    }
    const toml::Value& config_v = config_pr.value;

    // declare the input variables
    string title;
    string description;
    std::array<int,3> mesh_size;
    int n_vox;
    std::array<double,3> mesh_step;
    int n_dim;
    double freq;
    int n_tx_ch;
    int n_rx_ch;
    string tx_sens_fname_wc;
    string trx_phase_fname_wc;
    std::vector<std::vector<double> > tx_sens(0);
    std::vector<std::vector<double> > trx_phase(0);
    string reference_img_str;
    std::vector<double> reference_img;
    double *img;
    string sigma_fname;
    string epsr_fname;

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
        toml_get_array_of<int>(config_v,"mesh.size",mesh_size);
        cout<<"[ ";
        for (int n : mesh_size) {
            cout<<n<<", ";
        }
        cout<<"]\n";
        n_vox = std::accumulate(mesh_size.begin(),mesh_size.end(),1,std::multiplies<int>());
        cout<<"  number of voxels: "<<n_vox<<"\n";
        cout<<"  mesh step: ";
        toml_get_array_of<double>(config_v,"mesh.step",mesh_step);
        cout<<"[ ";
        for (double d : mesh_step) {
            cout<<d<<", ";
        }
        cout<<"]\n";
        n_dim = static_cast<int>(mesh_size.size());
        if (mesh_step.size() != n_dim) {
            throw runtime_error("format error: mesh.size and mesh.step must have the same size");
        }
        cout<<"  number of dimensions: "<<n_dim<<"\n";
        // ...input
        freq = config_v.get<double>("input.frequency");
        cout<<"\n  frequency: "<<freq<<"\n";
        n_tx_ch = config_v.get<int>("input.tx-channels");
        cout<<"  number of Tx channels: "<<n_tx_ch<<"\n";
        n_rx_ch = config_v.get<int>("input.rx-channels");
        cout<<"  number of Rx channels: "<<n_rx_ch<<"\n";
        tx_sens_fname_wc = config_v.get<string>("input.tx-sensitivity");
        cout<<"  Tx sensitivity wildcard: '"<<tx_sens_fname_wc<<"'\n";
        cout<<"    exp'd: [ ";
        for (int id_tx = 0; id_tx<n_tx_ch; ++id_tx) {
            string tx_sens_fname(tx_sens_fname_wc);
            replace(tx_sens_fname.begin(),tx_sens_fname.end(),chwc_tx,to_string(id_tx*chwc_step+chwc_start_from).c_str()[0]);
            ifstream ifile(tx_sens_fname,ios::binary);
            if (ifile.is_open()) {
                std::vector<double> tmp(n_vox);
                ifile.read(reinterpret_cast<char*>(tmp.data()),n_vox*sizeof(double));
                ifile.close();
                tx_sens.push_back(tmp);
            } else {
                cout<<"\n"
                    <<"ERROR: impossible read file '"<<tx_sens_fname<<"'"<<endl;
                return 2;
            }
            cout<<"'"<<tx_sens_fname<<"', ";
        }
        cout<<"]\n";
        trx_phase_fname_wc = config_v.get<string>("input.trx-phase");
        cout<<"  TRx phase wildcard: '"<<trx_phase_fname_wc<<"'\n";
        cout<<"    exp'd: [ ";
        for (int id_rx = 0; id_rx<n_rx_ch; ++id_rx) {
            for (int id_tx = 0; id_tx<n_tx_ch; ++id_tx) {
                string trx_phase_fname(trx_phase_fname_wc);
                replace(trx_phase_fname.begin(),trx_phase_fname.end(),chwc_tx,to_string(id_tx*chwc_step+chwc_start_from).c_str()[0]);
                replace(trx_phase_fname.begin(),trx_phase_fname.end(),chwc_rx,to_string(id_rx*chwc_step+chwc_start_from).c_str()[0]);
                ifstream ifile(trx_phase_fname,ios::binary);
                if (ifile.is_open()) {
                    std::vector<double> tmp(n_vox);
                    ifile.read(reinterpret_cast<char*>(tmp.data()),n_vox*sizeof(double));
                    ifile.close();
                    trx_phase.push_back(tmp);
                } else {
                    cout<<"\n"
                        <<"ERROR: impossible read file '"<<trx_phase_fname<<"'"<<endl;
                    return 2;
                }
                cout<<"'"<<trx_phase_fname<<"', ";
            }
        }
        cout<<"]\n";
        {
            const toml::Value* config_x;
            config_x = config_v.find("input.reference-img");
            if (config_x && config_x->is<string>()) {
                reference_img_str = config_x->as<string>().c_str();
                cout<<"  Reference image: '"<<reference_img_str<<"'\n";
                ifstream ifile(reference_img_str,ios::binary);
                if (ifile.is_open()) {
                    reference_img.resize(n_vox);
                    ifile.read(reinterpret_cast<char*>(reference_img.data()),n_vox*sizeof(double));
                    ifile.close();
                    img = reference_img.data();
                } else {
                    cout<<"\n"
                        <<"ERROR: impossible read file '"<<reference_img_str<<"'"<<endl;
                    return 2;
                }
            } else {
                reference_img_str = "";
                img = nullptr;
            }
        }
        // ...output
        sigma_fname = config_v.get<string>("output.electric-conductivity");
        cout<<"\n  Output electric conductivity: '"<<sigma_fname<<"'\n";
        epsr_fname = config_v.get<string>("output.relative-permittivity");
        cout<<"  Output relative permittivity: '"<<epsr_fname<<"'\n";
    } catch (runtime_error e) {
        cout<<"\n"
            <<"ERROR: "<<e.what()<<endl;
        return 2;
    }

    // load the files and perform the EPT
    std::array<int,3> rr = {1,1,1};
    eptlib::Shape kernel_shape = eptlib::shapes::Cross(rr);
    eptlib::EPTHelmholtz ept_helm(freq, mesh_size, mesh_step, kernel_shape);
    for (int id_tx = 0; id_tx<n_tx_ch; ++id_tx) {
        ept_helm.SetTxSensitivity(tx_sens[id_tx].data());
    }
    for (int id_rx = 0; id_rx<n_rx_ch; ++id_rx) {
        for (int id_tx = 0; id_tx<n_tx_ch; ++id_tx) {
            ept_helm.SetTRxPhase(trx_phase[id_tx+n_tx_ch*id_rx].data());
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
    EPTlibError_t eptlib_error = ept_helm.ApplyPostPro(img);
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
        ofstream ofile(sigma_fname,ios::binary);
        if (ofile.is_open()) {
            ofile.write(reinterpret_cast<char*>(sigma.data()), n_vox*sizeof(double));
            ofile.close();
        } else {
            cout<<"ERROR: impossible write file '"<<sigma_fname<<"'"<<endl;
            return 2;
        }
    }
    {
        ofstream ofile(epsr_fname,ios::binary);
        if (ofile.is_open()) {
            ofile.write(reinterpret_cast<char*>(epsr.data()), n_vox*sizeof(double));
            ofile.close();
        } else {
            cout<<"ERROR: impossible write file '"<<epsr_fname<<"'"<<endl;
            return 2;
        }
    }

    return 0;
}

template <typename T, typename U>
void toml_get_array_of(const toml::Value &v, const string &key, U &vec) {
    const toml::Array &a = v.get<toml::Array>(key);
    for (int d = 0; d<3; ++d) {
        vec[d] = a[d].as<T>();
    }
    return;
};
