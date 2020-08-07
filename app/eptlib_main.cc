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

#include <chrono>
#include <exception>
#include <iostream>
#include <memory>
#include <regex>
#include <utility>
#include <vector>

#include "eptlib/ept_methods.h"
#include "eptlib/io/io_hdf5.h"
#include "eptlib/io/io_toml.h"
#include "eptlib/version.h"

#include "eptlib_main.h"

#define LOADMANDATORY(what,io_toml,data,T) { \
    EPTlibError MACRO_error = io_toml->what<T>(data.first,data.second); \
    if (MACRO_error!=EPTlibError::Success) { \
        cout<<"FATAL ERROR in config file: "<<ToString(MACRO_error)+" '"+data.second+"'"<<endl; \
        return 1; \
    } \
}
#define LOADMANDATORYDATA(io_toml,data) LOADMANDATORY(GetValue,io_toml,data,decltype(data.first))
#define LOADMANDATORYLIST(io_toml,data) LOADMANDATORY(GetArrayOf,io_toml,data,decltype(data.first)::value_type)

#define LOADOPTIONAL(what,io_toml,data,T) { \
    EPTlibError MACRO_error = io_toml->what<T>(data.first,data.second); \
    if (WALL && MACRO_error!=EPTlibError::Success) { \
        cout<<"WARNING in config file: "<<ToString(MACRO_error)+" '"+data.second+"'"<<endl; \
    } \
}
#define LOADOPTIONALDATA(io_toml,data) LOADOPTIONAL(GetValue,io_toml,data,decltype(data.first))
#define LOADOPTIONALLIST(io_toml,data) LOADOPTIONAL(GetArrayOf,io_toml,data,decltype(data.first)::value_type)

#define LOADMAP(map,addr) { \
    string MACRO_fname; \
    string MACRO_uri; \
    io::GetAddress(addr,MACRO_fname,MACRO_uri); \
    io::IOh5 MACRO_ifile(MACRO_fname,io::Mode::In); \
    io::State MACRO_iostate = MACRO_ifile.ReadDataset(&map,"/",MACRO_uri); \
    if (MACRO_iostate!=io::State::Success) { \
        cout<<"FATAL ERROR: "<<ToString(MACRO_iostate)+" '"+addr+"'"<<endl; \
        return 1; \
    } \
}
#define SAVEMAP(map,addr) { \
    string MACRO_fname; \
    string MACRO_uri; \
    io::GetAddress(addr,MACRO_fname,MACRO_uri); \
    io::IOh5 MACRO_ofile(MACRO_fname,io::Mode::Append); \
    io::State MACRO_iostate = MACRO_ofile.WriteDataset(map,"/",MACRO_uri); \
    if (MACRO_iostate!=io::State::Success) { \
        cout<<"FATAL ERROR: "<<ToString(MACRO_iostate)+" '"+addr+"'"<<endl; \
        return 1; \
    } \
}

using namespace std;
using namespace eptlib;

constexpr char software::name[];

template <class T> using cfgdata = pair<T,string>;
template <class T> using cfglist = pair<array<T,NDIM>,string>;

int main(int argc, char **argv) {
    auto start = std::chrono::system_clock::now();
    // opening boilerplate
    cout<<project::str()<<" ("<<build::str()<<") ["<<compiler::str()<<"]\n"<<endl;
    cout<<LicenseBoilerplate()<<endl;
    // check the number of input
    if (argc<2) {
        cout<<"Usage example: "<<software::str()<<" <config file>"<<endl;
        return -1;
    }
    // load the config file
    unique_ptr<io::IOtoml> io_toml;
    try {
        io_toml.reset(new io::IOtoml(string(argv[1]),io::Mode::In));
    } catch(const ios_base::failure &e) {
        cout<<"FATAL ERROR in config file: "<<e.what()<<endl;
        return 1;
    }
    // declare the input variables
    //   mandatory input
    cfgdata<string> title; title.second = "title";
    cfgdata<string> descr; descr.second = "description";
    cfgdata<int> method; method.second = "method";
    cfglist<int> nn; nn.second = "mesh.size";
    cfglist<double> dd; dd.second = "mesh.step";
    cfgdata<double> freq; freq.second = "input.frequency";
    //   optional input
    cfgdata<int> n_txch(1,"input.tx-channels");
    cfgdata<int> n_rxch(1,"input.rx-channels");
    cfgdata<string> txsens_addr("","input.tx-sensitivity");
    cfgdata<string> trxphase_addr("","input.trx-phase");
    cfgdata<string> refimg_addr("","input.reference-img");
    cfgdata<string> sigma_addr("","output.electric-conductivity");
    cfgdata<string> epsr_addr("","output.relative-permittivity");
    // load the input data
    //   title
    LOADMANDATORYDATA(io_toml,title);
    LOADMANDATORYDATA(io_toml,descr);
    LOADMANDATORYDATA(io_toml,method);
    //   mesh
    LOADMANDATORYLIST(io_toml,nn);
    LOADMANDATORYLIST(io_toml,dd);
    //   input
    LOADMANDATORYDATA(io_toml,freq);
    LOADOPTIONALDATA(io_toml,n_txch);
    LOADOPTIONALDATA(io_toml,n_rxch);
    LOADOPTIONALDATA(io_toml,txsens_addr);
    LOADOPTIONALDATA(io_toml,trxphase_addr);
    LOADOPTIONALDATA(io_toml,refimg_addr);
    //   output
    LOADOPTIONALDATA(io_toml,sigma_addr);
    LOADOPTIONALDATA(io_toml,epsr_addr);
    cout<<endl;
    // check the provided data
    EPTMethod ept_method = static_cast<EPTMethod>(method.first);
    bool thereis_txsens = txsens_addr.first!="";
    bool thereis_trxphase = trxphase_addr.first!="";
    bool thereis_refimg = refimg_addr.first!="";
    bool thereis_sigma = sigma_addr.first!="";
    bool thereis_epsr = epsr_addr.first!="";
    //   EPT method
    if (method.first<0||ept_method>=EPTMethod::END) {
        cout<<"FATAL ERROR in config file: Wrong data format '"<<method.second<<"'"<<endl;
        return 1;
    }
    //   input addresses
    if (!thereis_txsens&&!thereis_trxphase) {
        cout<<"FATAL ERROR in config file: Neither '"<<txsens_addr.second<<"' nor '"<<trxphase_addr.second<<"' are provided"<<endl;
        return 1;
    }
    //   output addresses
    if (!thereis_sigma&&!thereis_epsr) {
        cout<<"FATAL ERROR in config file: Neither '"<<sigma_addr.second<<"' nor '"<<epsr_addr.second<<"' are provided"<<endl;
        return 1;
    }
    // report the readen values
    cout<<"  "<<title.first<<"\n";
    cout<<"  "<<descr.first<<"\n";
    cout<<"\n  Method: ("<<method.first<<") "<<ToString(ept_method)<<"\n";
    cout<<"\n  Mesh size: ["<<nn.first[0]<<", "<<nn.first[1]<<", "<<nn.first[2]<<"]\n";
    cout<<"  Mesh step: ["<<dd.first[0]<<", "<<dd.first[1]<<", "<<dd.first[2]<<"]\n";
    cout<<"\n  Frequency: "<<freq.first<<"\n";
    cout<<"  Tx channels: "<<n_txch.first<<"\n";
    cout<<"  Rx channels: "<<n_rxch.first<<"\n";
    cout<<"\n  Tx sensitivity addr.: '"<<txsens_addr.first<<"'\n";
    cout<<"  TRx phase addr.: '"<<trxphase_addr.first<<"'\n";
    cout<<"  Reference image addr.: '"<<refimg_addr.first<<"'\n";
    cout<<"\n  Output electric conductivity addr.: '"<<sigma_addr.first<<"'\n";
    cout<<"  Output relative permittivity addr.: '"<<epsr_addr.first<<"'\n";
    cout<<endl;
    // load the input maps
    vector<Image<double> > txsens(0);
    vector<Image<double> > trxphase(0);
    unique_ptr<Image<double> > refimg;
    //   look for alternative wildcards
    cfgdata<char> tx_wc('>',"input.wildcard.tx-character");
    cfgdata<char> rx_wc('<',"input.wildcard.rx-character");
    cfgdata<int> wc_start_from(0,"input.wildcard.start-from");
    cfgdata<int> wc_step(1,"input.wildcard.step");
    LOADOPTIONALDATA(io_toml,tx_wc);
    LOADOPTIONALDATA(io_toml,rx_wc);
    LOADOPTIONALDATA(io_toml,wc_start_from);
    LOADOPTIONALDATA(io_toml,wc_step);
    cout<<endl;
    //   check the provided wildcards
    if (tx_wc==rx_wc) {
        cout<<"FATAL ERROR in config file: wildcard characters for Tx and Rx are equal"<<endl;
        return 1;
    }
    //   report the wildcards
    cout<<"# Wildcard characters\n";
    cout<<"# Tx character: '"<<tx_wc.first<<"'\n";
    cout<<"# Rx character: '"<<rx_wc.first<<"'\n";
    cout<<"# Start value: "<<wc_start_from.first<<"\n";
    cout<<"# Step: "<<wc_step.first<<"\n";
    cout<<endl;
    //   load the maps
    if (thereis_txsens) {
        cout<<"Loading Tx sensitivity:\n"<<flush;
        for (int id_tx = 0; id_tx<n_txch.first; ++id_tx) {
            int tx = id_tx*wc_step.first+wc_start_from.first;
            string addr(txsens_addr.first);
            StringReplace(addr,string(1,tx_wc.first),to_string(tx));
            Image<double> map(nn.first[0],nn.first[1],nn.first[2]);
            LOADMAP(map,addr);
            txsens.push_back(map);
            cout<<"  '"<<addr<<"'\n"<<flush;
        }
        cout<<endl;
    }
    if (thereis_trxphase) {
        cout<<"Loading TRx phase:\n"<<flush;
        for (int id_rx = 0; id_rx<n_rxch.first; ++id_rx) {
            for (int id_tx = 0; id_tx<n_txch.first; ++id_tx) {
                int tx = id_tx*wc_step.first+wc_start_from.first;
                int rx = id_rx*wc_step.first+wc_start_from.first;
                string addr(trxphase_addr.first);
                StringReplace(addr,string(1,tx_wc.first),to_string(tx));
                StringReplace(addr,string(1,rx_wc.first),to_string(rx));
                Image<double> map(nn.first[0],nn.first[1],nn.first[2]);
                LOADMAP(map,addr);
                trxphase.push_back(map);
                cout<<"  '"<<addr<<"'\n"<<flush;
            }
        }
        cout<<endl;
    }
    if (thereis_refimg) {
        cout<<"Loading reference image:\n"<<flush;
        refimg.reset(new Image<double>(nn.first[0],nn.first[1],nn.first[2]));
        LOADMAP(*refimg,refimg_addr.first);
        cout<<"  '"<<refimg_addr.first<<"'\n"<<endl;
    }
    // set-up the specific parameters of EPT methods
    std::unique_ptr<EPTInterface> ept;
    switch (ept_method) {
        case EPTMethod::HELMHOLTZ: {
            // declare the parameters
            cfglist<int> rr({1,1,1},"parameter.savitzky-golay.size");
            cfgdata<int> shape(0,"parameter.savitzky-golay.shape");
            // load the parameters
            LOADOPTIONALLIST(io_toml,rr);
            LOADOPTIONALDATA(io_toml,shape);
            cout<<endl;
            // check the parameters
            KernelShape kernel_shape = static_cast<KernelShape>(shape.first);
            if (shape.first<0||kernel_shape>=KernelShape::END) {
                cout<<"FATAL ERROR in config file: Wrong data format '"<<shape.second<<"'"<<endl;
                return 1;
            }
            if (n_txch.first>1||n_rxch.first>1) {
                cout<<"FATAL ERROR in config file: 1 transmit/receive channel is needed by "<<ToString(ept_method)<<endl;
                return 1;
            }
            // report the parameters
            cout<<"Parameters:\n";
            cout<<"  Savitzky-Golay filter:\n";
            cout<<"    Kernel size: ["<<rr.first[0]<<", "<<rr.first[1]<<", "<<rr.first[2]<<"]\n";
            cout<<"    Kernel shape: ("<<shape.first<<") "<<ToString(kernel_shape)<<"\n";
            cout<<endl;
            // combine the parameters
            Shape kernel;
            switch (kernel_shape) {
                case KernelShape::CROSS:
                    kernel = shapes::Cross(rr.first);
                    break;
                case KernelShape::ELLIPSOID:
                    kernel = shapes::Ellipsoid(rr.first);
                    break;
                case KernelShape::CUBOID: {
                    kernel = shapes::CuboidR(rr.first);
                    break;
                }
            }
            // create the EPT method
            ept.reset(new EPTHelmholtz(freq.first,nn.first,dd.first,kernel));
            break;
        }
        case EPTMethod::CONVREACT: {
            // declare the parameters
            cfglist<int> rr({1,1,1},"parameter.savitzky-golay.size");
            cfgdata<int> shape(0,"parameter.savitzky-golay.shape");
            cfgdata<bool> is_3d(false,"parameter.volume-tomography");
            cfgdata<int> imaging_slice(nn.first[2]/2,"parameter.imaging-slice");
            cfgdata<double> dir_sigma(0.0,"parameter.dirichlet.electric-conductivity");
            cfgdata<double> dir_epsr(1.0,"parameter.dirichlet.relative-permittivity");
            cfgdata<bool> thereis_diff(false,"parameter.artificial-diffusion");
            cfgdata<double> diff_coeff(0.0,"parameter.artificial-diffusion-coefficient");
            // load the parameters
            LOADOPTIONALLIST(io_toml,rr);
            LOADOPTIONALDATA(io_toml,shape);
            LOADOPTIONALDATA(io_toml,is_3d);
            if (!is_3d.first) {
                LOADOPTIONALDATA(io_toml,imaging_slice);
            }
            LOADOPTIONALDATA(io_toml,dir_sigma);
            LOADOPTIONALDATA(io_toml,dir_epsr);
            LOADOPTIONALDATA(io_toml,thereis_diff);
            if (thereis_diff.first) {
                LOADOPTIONALDATA(io_toml,diff_coeff);
            }
            cout<<endl;
            // check the parameters
            KernelShape kernel_shape = static_cast<KernelShape>(shape.first);
            if (shape.first<0||kernel_shape>=KernelShape::END) {
                cout<<"FATAL ERROR in config file: Wrong data format '"<<shape.second<<"'"<<endl;
                return 1;
            }
            if (!is_3d.first) {
                if (imaging_slice.first<rr.first[2]||imaging_slice.first>nn.first[2]-1-rr.first[2]) {
                    cout<<"FATAL ERROR in config file: Wrong data format '"<<imaging_slice.second<<"'"<<endl;
                    return 1;
                }
            }
            if (dir_sigma.first<0.0) {
                cout<<"FATAL ERROR in config file: Wrong data format '"<<dir_sigma.second<<"'"<<endl;
                return 1;
            }
            if (dir_epsr.first<1.0) {
                cout<<"FATAL ERROR in config file: Wrong data format '"<<dir_epsr.second<<"'"<<endl;
                return 1;
            }
            if (n_txch.first>1||n_rxch.first>1) {
                cout<<"FATAL ERROR in config file: 1 transmit/receive channel is needed by "<<ToString(ept_method)<<endl;
                return 1;
            }
            // report the parameters
            cout<<"Parameters:\n";
            cout<<"  Savitzky-Golay filter:\n";
            cout<<"    Kernel size: ["<<rr.first[0]<<", "<<rr.first[1]<<", "<<rr.first[2]<<"]\n";
            cout<<"    Kernel shape: ("<<shape.first<<") "<<ToString(kernel_shape)<<"\n";
            cout<<"  Volume tomography: "<<(is_3d.first?"Yes":"No")<<"\n";
            if (!is_3d.first) {
                cout<<"  Imaging slice: "<<imaging_slice.first<<"\n";
            }
            cout<<"  Dirichlet:\n";
            cout<<"    Electric conductivity: "<<dir_sigma.first<<"\n";
            cout<<"    Relative permittivity: "<<dir_epsr.first<<"\n";
            cout<<"  Artificial diffusion: "<<(thereis_diff.first?"Yes":"No")<<"\n";
            if (thereis_diff.first) {
                cout<<"    Artificial diffusion coefficient: "<<diff_coeff.first<<"\n";
            }
            cout<<endl;
            // combine the parameters
            Shape kernel;
            switch (kernel_shape) {
                case KernelShape::CROSS:
                    kernel = shapes::Cross(rr.first);
                    break;
                case KernelShape::ELLIPSOID:
                    kernel = shapes::Ellipsoid(rr.first);
                    break;
                case KernelShape::CUBOID: {
                    kernel = shapes::CuboidR(rr.first);
                    break;
                }
            }
            // create the EPT method
            ept.reset(new EPTConvReact(freq.first,nn.first,dd.first,kernel));
            if (is_3d.first) {
                dynamic_cast<EPTConvReact*>(ept.get())->Toggle3D();
            } else {
                dynamic_cast<EPTConvReact*>(ept.get())->SelectSlice(imaging_slice.first);
            }
            dynamic_cast<EPTConvReact*>(ept.get())->SetDirichlet(dir_epsr.first,dir_sigma.first);
            if (thereis_diff.first) {
                dynamic_cast<EPTConvReact*>(ept.get())->SetArtificialDiffusion(diff_coeff.first);
            }
            break;
        }
        case EPTMethod::GRADIENT: {
                // declare the parameters
            cfglist<int> rr({1,1,1},"parameter.savitzky-golay.size");
            cfgdata<int> shape(0,"parameter.savitzky-golay.shape");
            cfgdata<bool> is_3d(false,"parameter.volume-tomography");
            cfgdata<int> imaging_slice(nn.first[2]/2,"parameter.imaging-slice");
            cfgdata<bool> full_run(true,"parameter.full-run");
            cfglist<int> seed_point_x0({0,0,0},"parameter.seed-point.coordinates");
            cfgdata<double> seed_point_sigma(0.0,"parameter.seed-point.electric-conductivity");
            cfgdata<double> seed_point_epsr(1.0,"parameter.seed-point.relative-permittivity");
            // load the parameters
            LOADOPTIONALLIST(io_toml,rr);
            LOADOPTIONALDATA(io_toml,shape);
            LOADOPTIONALDATA(io_toml,is_3d);
            if (!is_3d.first) {
                LOADOPTIONALDATA(io_toml,imaging_slice);
            }
            LOADOPTIONALDATA(io_toml,full_run);
            LOADOPTIONALLIST(io_toml,seed_point_x0);
            LOADOPTIONALDATA(io_toml,seed_point_sigma);
            LOADOPTIONALDATA(io_toml,seed_point_epsr);
            cout<<endl;
            // check the parameters
            KernelShape kernel_shape = static_cast<KernelShape>(shape.first);
            SeedPoint seed_point;
            seed_point.ijk = seed_point_x0.first;
            seed_point.sigma = seed_point_sigma.first;
            seed_point.epsr = seed_point_epsr.first;
            if (shape.first<0||kernel_shape>=KernelShape::END) {
                cout<<"FATAL ERROR in config file: Wrong data format '"<<shape.second<<"'"<<endl;
                return 1;
            }
            if (!is_3d.first) {
                if (imaging_slice.first<2*rr.first[2]||imaging_slice.first>nn.first[2]-1-2*rr.first[2]) {
                    cout<<"FATAL ERROR in config file: Wrong data format '"<<imaging_slice.second<<"'"<<endl;
                    return 1;
                }
            }
            for (int d = 0; d<NDIM; ++d) {
                if (seed_point.ijk[d]<0||seed_point.ijk[d]>=nn.first[d]) {
                    cout<<"FATAL ERROR in config file: Wrong data format '"<<seed_point_x0.second<<"'"<<endl;
                    return 1;
                }
            }
            if (seed_point.sigma<0.0) {
                cout<<"FATAL ERROR in config file: Wrong data format '"<<seed_point_sigma.second<<"'"<<endl;
                return 1;
            }
            if (seed_point.epsr<1.0) {
                cout<<"FATAL ERROR in config file: Wrong data format '"<<seed_point_epsr.second<<"'"<<endl;
                return 1;
            }
            if (n_txch.first<5) {
                cout<<"FATAL ERROR in config file: at least 5 transmit channels are needed by "<<ToString(ept_method)<<endl;
                return 1;
            }
            if (n_rxch.first>1) {
                cout<<"FATAL ERROR in config file: 1 receive channel is needed by "<<ToString(ept_method)<<endl;
                return 1;
            }
            // report the parameters
            cout<<"Parameters:\n";
            cout<<"  Savitzky-Golay filter:\n";
            cout<<"    Kernel size: ["<<rr.first[0]<<", "<<rr.first[1]<<", "<<rr.first[2]<<"]\n";
            cout<<"    Kernel shape: ("<<shape.first<<") "<<ToString(kernel_shape)<<"\n";
            cout<<"  Volume tomography: "<<(is_3d.first?"Yes":"No")<<"\n";
            if (!is_3d.first) {
                cout<<"  Imaging slice: "<<imaging_slice.first<<"\n";
            }
            cout<<"  Full run: "<<(full_run.first?"Yes":"No")<<"\n";
            if (full_run.first) {
                cout<<"  Seed point:\n";
                cout<<"    Coordinates: ["<<seed_point.ijk[0]<<", "<<seed_point.ijk[1]<<", "<<seed_point.ijk[2]<<"]\n";
                cout<<"    Electric conductivity: "<<seed_point.sigma<<"\n";
                cout<<"    Relative permittivity: "<<seed_point.epsr<<"\n";
            }
            cout<<endl;
            // combine the parameters
            Shape kernel;
            switch (kernel_shape) {
                case KernelShape::CROSS:
                    kernel = shapes::Cross(rr.first);
                    break;
                case KernelShape::ELLIPSOID:
                    kernel = shapes::Ellipsoid(rr.first);
                    break;
                case KernelShape::CUBOID: {
                    kernel = shapes::CuboidR(rr.first);
                    break;
                }
            }
            // create the EPT method
            ept.reset(new EPTGradient(freq.first,nn.first,dd.first,n_txch.first,kernel,!is_3d.first));
            if (!is_3d.first) {
                dynamic_cast<EPTGradient*>(ept.get())->SelectSlice(imaging_slice.first);
            }
            if (!full_run.first) {
                dynamic_cast<EPTGradient*>(ept.get())->SetRunMode(EPTGradientRun::LOCAL);
            }
            dynamic_cast<EPTGradient*>(ept.get())->ToggleSeedPoints();
            dynamic_cast<EPTGradient*>(ept.get())->AddSeedPoint(seed_point);
            break;
        }
    }
    // set-up the common parameters of the method
    if (thereis_txsens) {
        for (int id_tx = 0; id_tx<n_txch.first; ++id_tx) {
            ept->SetTxSensitivity(&(txsens[id_tx]),id_tx);
        }
    }
    if (thereis_trxphase) {
        for (int id_rx = 0; id_rx<n_rxch.first; ++id_rx) {
            for (int id_tx = 0; id_tx<n_txch.first; ++id_tx) {
                ept->SetTRxPhase(&(trxphase[id_tx+n_txch.first*id_rx]),id_tx,id_rx);
            }
        }
    }
    // run the method
    cout<<"Run "<<ToString(ept_method)<<"..."<<flush;
    auto run_start = std::chrono::system_clock::now();
    EPTlibError run_error = ept->Run();
    auto run_end = std::chrono::system_clock::now();
    if (run_error==EPTlibError::Success) {
        auto run_elapsed = std::chrono::duration_cast<std::chrono::seconds>(run_end-run_start);
        cout<<"done! ["<<run_elapsed.count()<<" s]\n";
    } else {
        cout<<"execution failed\n"<<endl;
        return 1;
    }
    cout<<endl;
    // get the results
    Image<double> map(nn.first[0],nn.first[1],nn.first[2]);
    if (thereis_sigma) {
        EPTlibError sigma_error = ept->GetElectricConductivity(&map);
        if (sigma_error==EPTlibError::Success) {
            SAVEMAP(map,sigma_addr.first);
        }
    }
    if (thereis_epsr) {
        EPTlibError epsr_error = ept->GetRelativePermittivity(&map);
        if (epsr_error==EPTlibError::Success) {
            SAVEMAP(map,epsr_addr.first);
        }
    }
    //
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    cout<<"Execution ended in "<<elapsed.count()<<" s\n";
    cout<<endl;
    return 0;
}
