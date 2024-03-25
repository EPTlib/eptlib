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

#include <omp.h>

#define LOADMANDATORY(what,io_toml,data,T) { \
    EPTlibError MACRO_error = io_toml->what<T>(&data.first,data.second); \
    if (MACRO_error!=EPTlibError::Success) { \
        cout<<"FATAL ERROR in config file: "<<ToString(MACRO_error)+" '"+data.second+"'"<<endl; \
        return 1; \
    } \
}
#define LOADMANDATORYDATA(io_toml,data) LOADMANDATORY(GetValue,io_toml,data,decltype(data.first))
#define LOADMANDATORYLIST(io_toml,data) LOADMANDATORY(GetArrayOf,io_toml,data,decltype(data.first)::value_type)

#define LOADOPTIONAL(what,io_toml,data,T) { \
    EPTlibError MACRO_error = io_toml->what<T>(&data.first,data.second); \
    if (WALL && MACRO_error!=EPTlibError::Success) { \
        cout<<"WARNING in config file: "<<ToString(MACRO_error)+" '"+data.second+"'"<<endl; \
    } \
}
#define LOADOPTIONALDATA(io_toml,data) LOADOPTIONAL(GetValue,io_toml,data,decltype(data.first))
#define LOADOPTIONALLIST(io_toml,data) LOADOPTIONAL(GetArrayOf,io_toml,data,decltype(data.first)::value_type)

#define LOADOPTIONALNOWARNING(what,io_toml,data,T) { \
    EPTlibError MACRO_error = io_toml->what<T>(&data.first,data.second); \
}
#define LOADOPTIONALNOWARNINGDATA(io_toml,data) LOADOPTIONALNOWARNING(GetValue,io_toml,data,decltype(data.first))
#define LOADOPTIONALNOWARNINGLIST(io_toml,data) LOADOPTIONALNOWARNING(GetArrayOf,io_toml,data,decltype(data.first)::value_type)

#define LOADMAP(map,addr) { \
    string MACRO_fname; \
    string MACRO_uri; \
    io::GetAddress(addr,MACRO_fname,MACRO_uri); \
    std::cout << "Opened file..." << std::endl; \
    io::IOh5 MACRO_ifile(MACRO_fname,io::Mode::In); \
    std::cout << "File opened" << std::endl; \
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
template <class T> using cfglist = pair<array<T,N_DIM>,string>;

EPTlibError TOMLGetSeedPoints(std::vector<SeedPoint> &seed_points,
    const io::IOtoml &io_toml, const std::string &uri) {
    EPTlibError error;
    toml::Array coordinates;
    toml::Array sigmas;
    toml::Array epsrs;
    error = io_toml.GetValue<toml::Array>(&coordinates, uri+".coordinates");
    if (error!=EPTlibError::Success) {
        return error;
    }
    error = io_toml.GetValue<toml::Array>(&sigmas, uri+".electric-conductivity");
    if (error!=EPTlibError::Success) {
        return error;
    }
    error = io_toml.GetValue<toml::Array>(&epsrs, uri+".relative-permittivity");
    if (error!=EPTlibError::Success) {
        return error;
    }
    int num_sp = coordinates.size();
    if (sigmas.size()!=num_sp || epsrs.size()!=num_sp) {
        return EPTlibError::WrongDataFormat;
    }
    seed_points.resize(num_sp);
    for (int idx = 0; idx<num_sp; ++idx) {
        const toml::Array &tmp = coordinates[idx].as<toml::Array>();
        if (!tmp[0].is<int>()) {
            return EPTlibError::WrongDataFormat;
        }
        for (int d = 0; d<N_DIM; ++d) {
            seed_points[idx].ijk[d] = tmp[d].as<int>();
        }
        seed_points[idx].sigma = sigmas[idx].as<double>();
        seed_points[idx].epsr = epsrs[idx].as<double>();
    }
    return EPTlibError::Success;
}

EPTlibError TOMLGetMultipleSGShapes(std::vector<int> &shapes,
    std::vector<std::array<int,N_DIM> > &sizes, const io::IOtoml &io_toml,
    const std::string &uri) {
    EPTlibError error;
    toml::Array toml_shapes;
    toml::Array toml_sizes;
    error = io_toml.GetValue<toml::Array>(&toml_shapes,uri+".shapes");
    if (error!=EPTlibError::Success) {
        return error;
    }
    error = io_toml.GetValue<toml::Array>(&toml_sizes,uri+".sizes");
    if (error!=EPTlibError::Success) {
        return error;
    }
    int num_sg = toml_shapes.size();
    if (toml_sizes.size()!=num_sg) {
        return EPTlibError::WrongDataFormat;
    }
    shapes.resize(num_sg);
    sizes.resize(num_sg);
    for (int idx = 0; idx<num_sg; ++idx) {
        shapes[idx] = toml_shapes[idx].as<int>();
        const toml::Array &tmp = toml_sizes[idx].as<toml::Array>();
        if (!tmp[0].is<int>()) {
            return EPTlibError::WrongDataFormat;
        }
        for (int d = 0; d<N_DIM; ++d) {
            sizes[idx][d] = tmp[d].as<int>();
        }
    }
    return EPTlibError::Success;
}

int PerformPostprocessing(const string &output_addr, const string &input_addr,
    const string &reference_addr, const string &variance_addr, const Shape &kernel,
    const double weight_param, const cfglist<int> &nn) {
    //
    Image<double> dst(nn.first[0], nn.first[1], nn.first[2]);
    Image<double> src(nn.first[0], nn.first[1], nn.first[2]);
    Image<double> ref;
    Image<double> unc;
    bool thereis_reference = false;
    bool thereis_uncertainty = false;
    //
    cout<<"  Loading source image:\n"<<flush;
    LOADMAP(src, input_addr);
    cout<<"    '"<<input_addr<<"'\n"<<endl;
    if (reference_addr!="") {
        ref = Image<double>(nn.first[0], nn.first[1], nn.first[2]);
        cout<<"  Loading reference image:\n"<<flush;
        LOADMAP(ref, reference_addr);
        cout<<"    '"<<reference_addr<<"'\n"<<endl;
        thereis_reference = true;
    }
    if (variance_addr!="") {
        unc = Image<double>(nn.first[0], nn.first[1], nn.first[2]);
        cout<<"  Loading variance map:\n"<<flush;
        LOADMAP(unc, variance_addr);
        cout<<"    '"<<variance_addr<<"'\n"<<endl;
        thereis_uncertainty = true;
    }
    //
    cout<<"Postprocessing ..."<<endl;
    auto postprocessing_start = std::chrono::system_clock::now();
    if (!thereis_reference && !thereis_uncertainty) {
        filter::Postprocessing(&dst, src, kernel);
    } else if (thereis_reference && !thereis_uncertainty) {
        filter::Postprocessing(&dst, src, kernel, &ref, nullptr, weight_param);
    } else if (!thereis_reference && thereis_uncertainty) {
        filter::Postprocessing(&dst, src, kernel, nullptr, &unc);
    } else {
        filter::Postprocessing(&dst, src, kernel, &ref, &unc, weight_param);
    }
    auto postprocessing_end = std::chrono::system_clock::now();
    auto postprocessing_elapsed = std::chrono::duration_cast<std::chrono::seconds>(postprocessing_end - postprocessing_start);
    cout<<"done! ["<<postprocessing_elapsed.count()<<" s]\n"<<endl;
    //
    SAVEMAP(dst, output_addr);
    return 0;
}

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
        io_toml = std::make_unique<io::IOtoml>(string(argv[1]), io::Mode::In);
    } catch(const ios_base::failure &e) {
        cout<<"FATAL ERROR in config file: "<<e.what()<<endl;
        return 1;
    }
    cout<<"Number of available threads: "<<omp_get_max_threads()<<endl;
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
    cfgdata<bool> wrapped_phase(false,"input.wrapped-phase");
    cfgdata<string> refimg_addr("","input.reference-image");
    cfgdata<string> sigma_addr("","output.electric-conductivity");
    cfgdata<string> epsr_addr("","output.relative-permittivity");
    cfgdata<string> postprocessing_output_sigma_addr("","postprocessing.output-electric-conductivity");
    cfgdata<string> postprocessing_output_epsr_addr("","postprocessing.output-relative-permittivity");
    cfgdata<bool> perform_only_postprocessing(false,"postprocessing.perform-only-postprocessing");
    cfgdata<string> postprocessing_input_sigma_addr("","postprocessing.input-electric-conductivity");
    cfgdata<string> postprocessing_input_epsr_addr("","postprocessing.input-relative-permittivity");
    cfgdata<string> postprocessing_input_reference_addr("","postprocessing.input-reference-image");
    cfgdata<string> postprocessing_input_sigma_variance_addr("","postprocessing.input-electric-conductivity-variance-map");
    cfgdata<string> postprocessing_input_epsr_variance_addr("","postprocessing.input-relative-permittivity-variance-map");
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
    LOADOPTIONALDATA(io_toml,wrapped_phase);
    LOADOPTIONALDATA(io_toml,refimg_addr);
    //   output
    LOADOPTIONALDATA(io_toml,sigma_addr);
    LOADOPTIONALDATA(io_toml,epsr_addr);
    //   postprocessing
    LOADOPTIONALDATA(io_toml,postprocessing_output_sigma_addr);
    LOADOPTIONALDATA(io_toml,postprocessing_output_epsr_addr);
    LOADOPTIONALDATA(io_toml,perform_only_postprocessing);
    if (perform_only_postprocessing.first) {
        LOADOPTIONALDATA(io_toml,postprocessing_input_sigma_addr);
        LOADOPTIONALDATA(io_toml,postprocessing_input_epsr_addr);
        LOADOPTIONALDATA(io_toml,postprocessing_input_reference_addr);
        LOADOPTIONALDATA(io_toml,postprocessing_input_sigma_variance_addr);
        LOADOPTIONALDATA(io_toml,postprocessing_input_epsr_variance_addr);
    }
    cout<<endl;
    // check the provided data
    EPTMethod ept_method = static_cast<EPTMethod>(method.first);
    bool thereis_txsens = txsens_addr.first!="";
    bool thereis_trxphase = trxphase_addr.first!="";
    bool thereis_refimg = refimg_addr.first!="";
    bool thereis_sigma = sigma_addr.first!="";
    bool thereis_epsr = epsr_addr.first!="";
    //   EPT method
    if (ept_method<=EPTMethod::BEGIN_STABLE||(ept_method>=EPTMethod::END_STABLE&&ept_method<=EPTMethod::BEGIN_EXPERIMENTAL)||ept_method>=EPTMethod::END_EXPERIMENTAL) {
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
    cout<<"  Phase is wrapped: "<<(wrapped_phase.first?"Yes":"No")<<"\n";
    cout<<"  Reference image addr.: '"<<refimg_addr.first<<"'\n";
    cout<<"\n  Output electric conductivity addr.: '"<<sigma_addr.first<<"'\n";
    cout<<"  Output relative permittivity addr.: '"<<epsr_addr.first<<"'\n";
    cout<<"\n  Postprocessing:\n";
    cout<<"    Output electric conductivity addr.: '"<<postprocessing_output_sigma_addr.first<<"'\n";
    cout<<"    Output relative permittivity addr.: '"<<postprocessing_output_epsr_addr.first<<"'\n";
    cout<<"    Perform only postprocessing: "<<(perform_only_postprocessing.first?"Yes":"No")<<"\n";
    if (perform_only_postprocessing.first) {
        cout<<"    Input electric conductivity addr.: '"<<postprocessing_input_sigma_addr.first<<"'\n";
        cout<<"    Input relative permittivity addr.: '"<<postprocessing_input_epsr_addr.first<<"'\n";
        cout<<"    Input reference image addr.: '"<<postprocessing_input_reference_addr.first<<"'\n";
        cout<<"    Input electric conductivity variance map addr.: '"<<postprocessing_input_sigma_variance_addr.first<<"'\n";
        cout<<"    Input relative permittivity variance map addr.: '"<<postprocessing_input_epsr_variance_addr.first<<"'\n";
    }
    cout<<endl;
    // perform EPT
    if (!perform_only_postprocessing.first) {
        // load the input maps
        vector<Image<double> > txsens(0);
        vector<Image<double> > trxphase(0);
        optional<Image<double> > refimg;
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
                StringReplace(&addr,string(1,tx_wc.first),to_string(tx));
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
                    StringReplace(&addr,string(1,tx_wc.first),to_string(tx));
                    StringReplace(&addr,string(1,rx_wc.first),to_string(rx));
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
            refimg.emplace(nn.first[0],nn.first[1],nn.first[2]);
            LOADMAP(*refimg,refimg_addr.first);
            cout<<"  '"<<refimg_addr.first<<"'\n"<<endl;
        }
        // set-up the specific parameters of EPT methods
        unique_ptr<EPTInterface> ept = nullptr;
        switch (ept_method) {
            case EPTMethod::HELMHOLTZ: {
                // declare the parameters
                cfglist<int> rr({1,1,1},"parameter.savitzky-golay.size");
                cfgdata<int> shape(0,"parameter.savitzky-golay.shape");
                cfgdata<int> degree(2,"parameter.savitzky-golay.degree");
                cfgdata<double> weight_param(0.05,"parameter.savitzky-golay.weight-param");
                cfgdata<string> output_var_addr("","parameter.output-variance");
                // load the parameters
                LOADOPTIONALLIST(io_toml,rr);
                LOADOPTIONALDATA(io_toml,shape);
                LOADOPTIONALDATA(io_toml,degree);
                LOADOPTIONALDATA(io_toml,output_var_addr);
                LOADOPTIONALDATA(io_toml,weight_param);
                cout<<endl;
                // check the parameters
                KernelShape kernel_shape = static_cast<KernelShape>(shape.first);
                if (shape.first<0||kernel_shape>=KernelShape::END) {
                    cout<<"FATAL ERROR in config file: Wrong data format '"<<shape.second<<"'"<<endl;
                    return 1;
                }
                if (degree.first<2) {
                    cout<<"FATAL ERROR in config file: Wrong data format '"<<degree.second<<"'"<<endl;
                    return 1;
                }
                if (n_txch.first>1||n_rxch.first>1) {
                    cout<<"FATAL ERROR in config file: 1 transmit/receive channel is needed by "<<ToString(ept_method)<<endl;
                    return 1;
                }
                if (output_var_addr.first!="" && wrapped_phase.first) {
                    cout<<"WARNING: Variance cannot be evaluated with wrapped phase. The phase will be assumed unwrapped."<<endl;
                }
                cout<<endl;
                // report the parameters
                cout<<"Parameters:\n";
                cout<<"  Savitzky-Golay filter:\n";
                cout<<"    Kernel size: ["<<rr.first[0]<<", "<<rr.first[1]<<", "<<rr.first[2]<<"]\n";
                cout<<"    Kernel shape: ("<<shape.first<<") "<<ToString(kernel_shape)<<"\n";
                cout<<"    Polynomial degree: "<<degree.first<<"\n";
                cout<<"    Weight parameter: "<<weight_param.first<<"\n";
                cout<<"  Output variance addr.: '"<<output_var_addr.first<<"'\n";
                cout<<endl;
                // combine the parameters
                Shape kernel;
                switch (kernel_shape) {
                    case KernelShape::CROSS:
                        kernel = shapes::Cross(rr.first[0],rr.first[1],rr.first[2]);
                        break;
                    case KernelShape::ELLIPSOID:
                        kernel = shapes::Ellipsoid(rr.first[0],rr.first[1],rr.first[2]);
                        break;
                    case KernelShape::CUBOID: {
                        kernel = shapes::CuboidR(rr.first[0],rr.first[1],rr.first[2]);
                        break;
                    }
                }
                // create the EPT method
                bool compute_variance = output_var_addr.first!="";
                ept = std::make_unique<EPTHelmholtz>(nn.first[0],nn.first[1],nn.first[2], dd.first[0],dd.first[1],dd.first[2], freq.first, kernel, degree.first, wrapped_phase.first, compute_variance, weight_param.first);
                break;
            }
            case EPTMethod::CONVREACT: {
                // declare the parameters
                cfglist<int> rr({1,1,1},"parameter.savitzky-golay.size");
                cfgdata<int> shape(0,"parameter.savitzky-golay.shape");
                cfgdata<int> degree(2,"parameter.savitzky-golay.degree");
                cfgdata<double> weight_param(0.05,"parameter.savitzky-golay.weight-param");
                cfgdata<bool> is_3d(false,"parameter.volume-tomography");
                cfgdata<int> imaging_slice(nn.first[2]/2,"parameter.imaging-slice");
                cfgdata<double> dir_sigma(0.0,"parameter.dirichlet.electric-conductivity");
                cfgdata<double> dir_epsr(1.0,"parameter.dirichlet.relative-permittivity");
                cfgdata<bool> thereis_diff(false,"parameter.artificial-diffusion");
                cfgdata<double> diff_coeff(0.0,"parameter.artificial-diffusion-coefficient");
                cfgdata<int> max_iterations(1000,"parameter.max-iterations");
                cfgdata<double> tolerance(1e-6,"parameter.tolerance");
                // load the parameters
                LOADOPTIONALLIST(io_toml,rr);
                LOADOPTIONALDATA(io_toml,shape);
                LOADOPTIONALDATA(io_toml,degree);
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
                LOADOPTIONALDATA(io_toml,max_iterations);
                LOADOPTIONALDATA(io_toml,tolerance);
                LOADOPTIONALDATA(io_toml,weight_param);
                cout<<endl;
                // check the parameters
                KernelShape kernel_shape = static_cast<KernelShape>(shape.first);
                if (shape.first<0||kernel_shape>=KernelShape::END) {
                    cout<<"FATAL ERROR in config file: Wrong data format '"<<shape.second<<"'"<<endl;
                    return 1;
                }
                if (degree.first<2) {
                    cout<<"FATAL ERROR in config file: Wrong data format '"<<degree.second<<"'"<<endl;
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
                cout<<"    Polynomial degree: "<<degree.first<<"\n";
                cout<<"    Weight parameter: "<<weight_param.first<<"\n";
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
                cout<<"  Max iterations: "<<max_iterations.first<<"\n";
                cout<<"  Tolerance: "<<tolerance.first<<"\n";
                cout<<endl;
                // combine the parameters
                Shape kernel;
                switch (kernel_shape) {
                    case KernelShape::CROSS:
                        kernel = shapes::Cross(rr.first[0], rr.first[1], rr.first[2]);
                        break;
                    case KernelShape::ELLIPSOID:
                        kernel = shapes::Ellipsoid(rr.first[0], rr.first[1], rr.first[2]);
                        break;
                    case KernelShape::CUBOID:
                        kernel = shapes::CuboidR(rr.first[0], rr.first[1], rr.first[2]);
                        break;
                }
                // create the EPT method
                ept = std::make_unique<EPTConvReact>(nn.first[0],nn.first[1],nn.first[2], dd.first[0],dd.first[1],dd.first[2], freq.first, kernel, degree.first, max_iterations.first, tolerance.first, weight_param.first);
                if (!is_3d.first) {
                    dynamic_cast<EPTConvReact*>(ept.get())->SelectSlice(imaging_slice.first);
                }
                if (thereis_diff.first) {
                    dynamic_cast<EPTConvReact*>(ept.get())->SetArtificialDiffusion(diff_coeff.first);
                }
                dynamic_cast<EPTConvReact*>(ept.get())->SetDirichlet(dir_epsr.first,dir_sigma.first);
                break;
            }
            case EPTMethod::GRADIENT: {
                // declare the parameters
                cfglist<int> rr({1,1,1},"parameter.savitzky-golay.size");
                cfgdata<int> shape(0,"parameter.savitzky-golay.shape");
                cfgdata<int> degree(2,"parameter.savitzky-golay.degree");
                cfgdata<double> weight_param(0.05,"parameter.savitzky-golay.weight-param");
                cfgdata<bool> is_3d(false,"parameter.volume-tomography");
                cfgdata<int> imaging_slice(nn.first[2]/2,"parameter.imaging-slice");
                cfgdata<bool> full_run(true,"parameter.full-run");
                cfgdata<string> output_gradient_addr("","parameter.output-gradient");
                cfgdata<bool> use_seed_point(false, "parameter.seed-point.use-seed-point");
                std::string seed_point_url = "parameter.seed-point";
                cfgdata<double> regularization_coefficient(1.0,"parameter.regularization.regularization-coefficient");
                cfgdata<double> regularization_gradient_tolerance(0.0,"parameter.regularization.gradient-tolerance");
                cfgdata<string> regularization_output_mask_addr("","parameter.regularization.output-mask");
                // load the parameters
                LOADOPTIONALLIST(io_toml,rr);
                LOADOPTIONALDATA(io_toml,shape);
                LOADOPTIONALDATA(io_toml,degree);
                LOADOPTIONALDATA(io_toml,is_3d);
                if (!is_3d.first) {
                    LOADOPTIONALDATA(io_toml,imaging_slice);
                }
                LOADOPTIONALDATA(io_toml,full_run);
                LOADOPTIONALDATA(io_toml,output_gradient_addr);
                LOADOPTIONALDATA(io_toml,use_seed_point);
                LOADOPTIONALDATA(io_toml,regularization_coefficient);
                LOADOPTIONALDATA(io_toml,regularization_gradient_tolerance);
                LOADOPTIONALDATA(io_toml,regularization_output_mask_addr);
                LOADOPTIONALDATA(io_toml,weight_param);
                cout<<endl;
                // check the parameters
                KernelShape kernel_shape = static_cast<KernelShape>(shape.first);
                std::vector<SeedPoint> seed_points(0);
                if (use_seed_point.first) {
                    EPTlibError error = TOMLGetSeedPoints(seed_points,*io_toml,seed_point_url);
                    if (error != EPTlibError::Success) {
                        cout<<"FATAL ERROR in config file: '"<<seed_point_url<<"'"<<endl;
                        return 1;
                    }
                }
                if (shape.first<0||kernel_shape>=KernelShape::END) {
                    cout<<"FATAL ERROR in config file: Wrong data format '"<<shape.second<<"'"<<endl;
                    return 1;
                }
                if (degree.first<2) {
                    cout<<"FATAL ERROR in config file: Wrong data format '"<<degree.second<<"'"<<endl;
                    return 1;
                }
                if (!is_3d.first) {
                    if (imaging_slice.first<2*rr.first[2]||imaging_slice.first>nn.first[2]-1-2*rr.first[2]) {
                        cout<<"FATAL ERROR in config file: Wrong data format '"<<imaging_slice.second<<"'"<<endl;
                        return 1;
                    }
                }
                for (int idx_sp = 0; idx_sp<seed_points.size(); ++idx_sp) {
                    for (int d = 0; d<N_DIM; ++d) {
                        if (seed_points[idx_sp].ijk[d]<0||seed_points[idx_sp].ijk[d]>=nn.first[d]) {
                            cout<<"FATAL ERROR in config file: Wrong data format '"<<seed_point_url<<"'"<<endl;
                            return 1;
                        }
                    }
                    if (seed_points[idx_sp].sigma<0.0) {
                        cout<<"FATAL ERROR in config file: Wrong data format '"<<seed_point_url<<"'"<<endl;
                        return 1;
                    }
                    if (seed_points[idx_sp].epsr<1.0) {
                        cout<<"FATAL ERROR in config file: Wrong data format '"<<seed_point_url<<"'"<<endl;
                        return 1;
                    }
                }
                if (regularization_coefficient.first<0.0) {
                    cout<<"FATAL ERROR in config file: Wrong data format '"<<regularization_coefficient.second<<"'"<<endl;
                    return 1;
                }
                if (regularization_gradient_tolerance.first<0.0 || regularization_gradient_tolerance.first>1.0) {
                    cout<<"FATAL ERROR in config file: Wrong data format '"<<regularization_gradient_tolerance.second<<"'"<<endl;
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
                cout<<"    Polynomial degree: "<<degree.first<<"\n";
                cout<<"  Volume tomography: "<<(is_3d.first?"Yes":"No")<<"\n";
                if (!is_3d.first) {
                    cout<<"  Imaging slice: "<<imaging_slice.first<<"\n";
                }
                cout<<"  Full run: "<<(full_run.first?"Yes":"No")<<"\n";
                cout<<"  Output gradient addr.: '"<<output_gradient_addr.first<<"'\n";
                if (full_run.first) {
                    cout<<"  Seed point:\n";
                    cout<<"    Use seed point: "<<(use_seed_point.first?"Yes":"No")<<"\n";
                    if (use_seed_point.first) {
                        cout<<"    Coordinates: [";
                        for (int idx_sp = 0; idx_sp<seed_points.size(); ++idx_sp) {
                            cout<<"["<<seed_points[idx_sp].ijk[0]<<", "<<seed_points[idx_sp].ijk[1]<<", "<<seed_points[idx_sp].ijk[2]<<"], ";
                        }
                        cout<<"]\n";
                        cout<<"    Electric conductivity: [";
                        for (int idx_sp = 0; idx_sp<seed_points.size(); ++idx_sp) {
                            cout<<seed_points[idx_sp].sigma<<", ";
                        }
                        cout<<"]\n";
                        cout<<"    Relative permittivity: [";
                        for (int idx_sp = 0; idx_sp<seed_points.size(); ++idx_sp) {
                            cout<<seed_points[idx_sp].epsr<<", ";
                        }
                        cout<<"]\n";
                    } else {
                        cout<<"  Regularization:\n";
                        cout<<"    Regularization coefficient: "<<regularization_coefficient.first<<"\n";
                        cout<<"    Gradient tolerance: "<<regularization_gradient_tolerance.first<<"\n";
                        cout<<"    Output mask addr.: '"<<regularization_output_mask_addr.first<<"'\n";
                    }
                }
                cout<<"  Weight parameter: "<<weight_param.first<<"\n";
                cout<<endl;
                // combine the parameters
                Shape kernel;
                switch (kernel_shape) {
                    case KernelShape::CROSS:
                        kernel = shapes::Cross(rr.first[0], rr.first[1], rr.first[2]);
                        break;
                    case KernelShape::ELLIPSOID:
                        kernel = shapes::Ellipsoid(rr.first[0], rr.first[1], rr.first[2]);
                        break;
                    case KernelShape::CUBOID: {
                        kernel = shapes::CuboidR(rr.first[0], rr.first[1], rr.first[2]);
                        break;
                    }
                }
                // create the EPT method
                EPTGradientRun run_mode = full_run.first ? EPTGradientRun::FULL : EPTGradientRun::LOCAL;
                ept = std::make_unique<EPTGradient>(nn.first[0],nn.first[1],nn.first[2], dd.first[0],dd.first[1],dd.first[2], freq.first, kernel, n_txch.first, degree.first, run_mode, weight_param.first);
                if (!is_3d.first) {
                    dynamic_cast<EPTGradient*>(ept.get())->SelectSlice(imaging_slice.first);
                }
                if (use_seed_point.first) {
                    for (int idx_sp = 0; idx_sp<seed_points.size(); ++idx_sp) {
                        dynamic_cast<EPTGradient*>(ept.get())->AddSeedPoint(seed_points[idx_sp]);
                    }
                } else {
                    dynamic_cast<EPTGradient*>(ept.get())->SetRegularizationCoefficient(regularization_coefficient.first);
                    dynamic_cast<EPTGradient*>(ept.get())->SetGradientTolerance(regularization_gradient_tolerance.first);
                }
                break;
            }
            case EPTMethod::HELMHOLTZ_CHI2: {
                // declare the parameters
                string savitzky_golay_url = "parameter.savitzky-golay";
                cfgdata<int> degree(2,"parameter.savitzky-golay.degree");
                cfgdata<string> output_sg_index_addr("","parameter.savitzky-golay.output-index");
                cfgdata<bool> admit_unphysical_values(false,"parameter.unphysical-values");
                cfgdata<string> output_var_addr("","parameter.output-variance");
                // load the parameters
                LOADOPTIONALDATA(io_toml,degree);
                LOADOPTIONALDATA(io_toml,output_sg_index_addr);
                LOADOPTIONALDATA(io_toml,admit_unphysical_values);
                LOADOPTIONALDATA(io_toml,output_var_addr);
                cout<<endl;
                std::vector<int> shapes(0);
                std::vector<std::array<int,N_DIM> > sizes(0);
                EPTlibError error = TOMLGetMultipleSGShapes(shapes,sizes,*io_toml,savitzky_golay_url);
                if (error!=EPTlibError::Success) {
                    cout<<"FATAL ERROR in config file: '"<<savitzky_golay_url<<"'"<<endl;
                    return 1;
                }
                // check the parameters
                std::vector<KernelShape> kernel_shapes(shapes.size());
                for (int idx = 0; idx<shapes.size(); ++idx) {
                    kernel_shapes[idx] = static_cast<KernelShape>(shapes[idx]);
                    if (shapes[idx]<0||kernel_shapes[idx]>=KernelShape::END) {
                        cout<<"FATAL ERROR in config file: Wrong data format '"<<savitzky_golay_url<<"'"<<endl;
                        return 1;
                    }
                }
                if (n_txch.first>1||n_rxch.first>1) {
                    cout<<"FATAL ERROR in config file: 1 transmit/receive channel is needed by "<<ToString(ept_method)<<endl;
                    return 1;
                }
                if (!thereis_trxphase) {
                    cout<<"FATAL ERROR in config file: The transceive phase address is needed by "<<ToString(ept_method)<<endl;
                    return 1;
                }
                if (degree.first<2) {
                    cout<<"FATAL ERROR in config file: Wrong data format '"<<degree.second<<"'"<<endl;
                    return 1;
                }
                if (thereis_txsens) {
                    cout<<"WARNING: This method works only with the phase-based approximation. Relative permittivity will not be computed and the Tx sensitivity will not be used."<<endl;
                }
                if (wrapped_phase.first) {
                    cout<<"WARNING: Variance cannot be evaluated with wrapped phase. The phase will be assumed unwrapped."<<endl;
                }
                cout<<endl;
                // report the parameters
                cout<<"Parameters:\n";
                cout<<"  Savitzky-Golay filter:\n";
                cout<<"    Kernel shapes: ["<<shapes[0];
                for (int idx = 1; idx<shapes.size(); ++idx) {
                    cout<<", "<<shapes[idx];
                }
                cout<<"]\n";
                cout<<"    Kernel sizes: [";
                cout<<"["<<sizes[0][0]<<", "<<sizes[0][1]<<", "<<sizes[0][2]<<"]";
                for (int idx = 1; idx<sizes.size(); ++idx) {
                    cout<<", "<<"["<<sizes[idx][0]<<", "<<sizes[idx][1]<<", "<<sizes[idx][2]<<"]";
                }
                cout<<"]\n";
                cout<<"    Polynomial degree: "<<degree.first<<"\n";
                cout<<"    Output index addr.: '"<<output_sg_index_addr.first<<"'\n";
                cout<<"  Admit unphysical values: "<<(admit_unphysical_values.first?"Yes":"No")<<"\n";
                cout<<"  Output variance addr.: '"<<output_var_addr.first<<"'\n";
                cout<<endl;
                // combine the parameters
                std::vector<Shape> kernels(shapes.size());
                for (int idx = 0; idx<kernels.size(); ++idx) {
                    switch (kernel_shapes[idx]) {
                        case KernelShape::CROSS:
                            kernels[idx] = shapes::Cross(sizes[idx][0],sizes[idx][1],sizes[idx][2]);
                            break;
                        case KernelShape::ELLIPSOID:
                            kernels[idx] = shapes::Ellipsoid(sizes[idx][0],sizes[idx][1],sizes[idx][2]);
                            break;
                        case KernelShape::CUBOID: {
                            kernels[idx] = shapes::CuboidR(sizes[idx][0],sizes[idx][1],sizes[idx][2]);
                            break;
                        }
                    }
                }
                // create the EPT method
                ept = std::make_unique<EPTHelmholtzChi2>(nn.first[0],nn.first[1],nn.first[2], dd.first[0],dd.first[1],dd.first[2], freq.first, kernels, degree.first, admit_unphysical_values.first);
                break;
            }
        }
        // set-up the common parameters of the method
        if (thereis_txsens) {
            for (int id_tx = 0; id_tx<n_txch.first; ++id_tx) {
                ept->SetTxSensitivity(txsens[id_tx], id_tx);
            }
        }
        if (thereis_trxphase) {
            for (int id_rx = 0; id_rx<n_rxch.first; ++id_rx) {
                for (int id_tx = 0; id_tx<n_txch.first; ++id_tx) {
                    ept->SetTRxPhase(trxphase[id_tx+n_txch.first*id_rx], id_tx,id_rx);
                }
            }
        }
        if (thereis_refimg) {
            ept->SetReferenceImage(refimg.value());
        }
        // run the method
        cout<<"Run "<<ToString(ept_method)<<"..."<<flush;
        auto run_start = std::chrono::system_clock::now();
        EPTlibError run_error = ept->Run();
        auto run_end = std::chrono::system_clock::now();
        if (run_error==EPTlibError::Success) {
            auto run_elapsed = std::chrono::duration_cast<std::chrono::seconds>(run_end-run_start);
            cout<<"done! ["<<run_elapsed.count()<<" s]\n"<<flush;
        } else {
            cout<<"execution failed\n"<<endl;
            return 1;
        }
        // report possible output parameters
        if (ept_method==EPTMethod::HELMHOLTZ) {
            // Save the variance of the recovered maps
            cfgdata<string> output_var_addr("","parameter.output-variance");
            LOADOPTIONALNOWARNINGDATA(io_toml,output_var_addr);
            bool thereis_var = output_var_addr.first!="";
            if (thereis_var) {
                if (dynamic_cast<EPTHelmholtz*>(ept.get())->ThereIsElectricConductivityVariance()) {
                    auto& var = dynamic_cast<EPTHelmholtz*>(ept.get())->GetElectricConductivityVariance();
                    SAVEMAP(*var, output_var_addr.first + "/electric-conductivity");
                }
                if (dynamic_cast<EPTHelmholtz*>(ept.get())->ThereIsRelativePermittivityVariance()) {
                    auto& var = dynamic_cast<EPTHelmholtz*>(ept.get())->GetRelativePermittivityVariance();
                    SAVEMAP(*var, output_var_addr.first + "/relative-permittivity");
                }
            }
        } else if (ept_method==EPTMethod::CONVREACT) {
            cout<<"  Iterative solver parameters:\n";
            cout<<"    Iterations: "<<dynamic_cast<EPTConvReact*>(ept.get())->GetSolverIterations()<<"\n";
            cout<<"    Estimated error: "<<dynamic_cast<EPTConvReact*>(ept.get())->GetSolverResidual()<<std::endl;
        } else if (ept_method==EPTMethod::GRADIENT) {
            if (dynamic_cast<EPTGradient*>(ept.get())->GetRunMode()==EPTGradientRun::FULL) {
                cout<<"  Gradient inversion parameters:\n";
                cout<<"    Cost functional: "<<dynamic_cast<EPTGradient*>(ept.get())->GetCostFunctional()<<"\n";
                if (!dynamic_cast<EPTGradient*>(ept.get())->SeedPointsAreUsed()) {
                    cout<<"    Cost regularization: "<<dynamic_cast<EPTGradient*>(ept.get())->GetCostRegularization()<<"\n";
                    // Save the mask
                    {
                        cfgdata<string> regularization_output_mask_addr("","parameter.regularization.output-mask");
                        LOADOPTIONALNOWARNINGDATA(io_toml,regularization_output_mask_addr);
                        bool thereis_mask = regularization_output_mask_addr.first!="";
                        if (thereis_mask) {
                            const boost::dynamic_bitset<>& mask = dynamic_cast<EPTGradient*>(ept.get())->GetMask();
                            std::vector<int> nn_mask{nn.first[0]-1, nn.first[1]-1, 1};
                            if (mask.size()!=nn_mask[0]*nn_mask[1]) {
                                nn_mask[2] = nn.first[2]-1;
                            }
                            Image<int> mask_img(nn_mask[0], nn_mask[1], nn_mask[2]);
                            for (int idx = 0; idx<mask_img.GetNVox(); ++idx) {
                                if (mask[idx]) {
                                    mask_img(idx) = 1;
                                } else {
                                    mask_img(idx) = 0;
                                }
                            }
                            SAVEMAP(mask_img,regularization_output_mask_addr.first);
                        }
                    }
                }
            }
            // Save the gradient distribution 
            {
                cfgdata<bool> is_3d(false,"parameter.volume-tomography");
                cfgdata<string> output_gradient_addr("","parameter.output-gradient");
                LOADOPTIONALNOWARNINGDATA(io_toml,is_3d);
                LOADOPTIONALNOWARNINGDATA(io_toml,output_gradient_addr);
                bool thereis_grad = output_gradient_addr.first!="";
                if (thereis_grad) {
                    std::vector<int> nn_grad{nn.first[0], nn.first[1], 1};
                    if (is_3d.first) {
                        nn_grad[2] = nn.first[2];
                    }
                    Image<double> grad_real(nn_grad[0], nn_grad[1], nn_grad[2]);
                    Image<double> grad_imag(nn_grad[0], nn_grad[1], nn_grad[2]);
                    // GPlus
                    auto& grad_p = dynamic_cast<EPTGradient*>(ept.get())->GetGPlus();
                    if (grad_p.size() > 0) {
                        for (int idx = 0; idx<grad_imag.GetNVox(); ++idx) {
                            grad_real(idx) = grad_p[idx].real();
                            grad_imag(idx) = grad_p[idx].imag();
                        }
                        string real_addr = output_gradient_addr.first+"/gplus/real";
                        string imag_addr = output_gradient_addr.first+"/gplus/imag";
                        SAVEMAP(grad_real,real_addr);
                        SAVEMAP(grad_imag,imag_addr);
                    }
                    // GZ
                    auto& grad_z = dynamic_cast<EPTGradient*>(ept.get())->GetGZ();
                    if (grad_z.size() > 0) {
                        for (int idx = 0; idx<grad_imag.GetNVox(); ++idx) {
                            grad_real(idx) = grad_z[idx].real();
                            grad_imag(idx) = grad_z[idx].imag();
                        }
                        string real_addr = output_gradient_addr.first+"/gz/real";
                        string imag_addr = output_gradient_addr.first+"/gz/imag";
                        SAVEMAP(grad_real,real_addr);
                        SAVEMAP(grad_imag,imag_addr);
                    }
                }
            }
        } else if (ept_method==EPTMethod::HELMHOLTZ_CHI2) {
            // Save the quality index chi2 distribution
            cfgdata<string> output_sg_index_addr("","parameter.savitzky-golay.output-index");
            LOADOPTIONALNOWARNINGDATA(io_toml,output_sg_index_addr);
            bool thereis_index = output_sg_index_addr.first!="";
            if (thereis_index) {
                auto& index = dynamic_cast<EPTHelmholtzChi2*>(ept.get())->GetIndex();
                SAVEMAP(*index,output_sg_index_addr.first);
            }
            cfgdata<string> output_var_addr("","parameter.output-variance");
            LOADOPTIONALNOWARNINGDATA(io_toml,output_var_addr);
            bool thereis_var = output_var_addr.first!="";
            if (thereis_var) {
                auto& var = dynamic_cast<EPTHelmholtzChi2*>(ept.get())->GetElectricConductivityVariance();
                SAVEMAP(*var,output_var_addr.first + "/electric-conductivity");
            }
        }
        cout<<endl;
        // get the results
        if (thereis_sigma) {
            auto& map = ept->GetElectricConductivity();
            SAVEMAP(*map,sigma_addr.first);
        }
        if (thereis_epsr) {
            auto& map = ept->GetRelativePermittivity();
            SAVEMAP(*map,epsr_addr.first);
        }
    }
    //
    if (postprocessing_output_sigma_addr.first!="") {
        string output = postprocessing_output_sigma_addr.first;
        string input;
        string reference;
        string variance;
        bool postprocess_flag = true;
        if (perform_only_postprocessing.first) {
            if (postprocessing_input_sigma_addr.first!="") {
                input = postprocessing_input_sigma_addr.first;
                reference = postprocessing_input_reference_addr.first;
                variance = postprocessing_input_sigma_variance_addr.first;
            } else {
                cout<<"WARNING: No input conductivity to be postprocessed!"<<endl;
                postprocess_flag = false;
            }
        } else {
            if (sigma_addr.first!="") {
                cfgdata<string> output_var_addr("","parameter.output-variance");
                LOADOPTIONALNOWARNINGDATA(io_toml, output_var_addr);
                input = sigma_addr.first;
                reference = refimg_addr.first;
                variance = output_var_addr.first + "/electric-conductivity";
            } else {
                cout<<"WARNING: No input conductivity to be postprocessed!"<<endl;
                postprocess_flag = false;
            }
        }
        if (postprocess_flag) {
            cfglist<int> rr({1,1,1},"postprocessing.kernel-size");
            cfgdata<int> shape(2,"postprocessing.kernel-shape");
            cfgdata<double> weight_param(0.10,"postprocessing.weight-param");
            LOADOPTIONALNOWARNINGLIST(io_toml,rr);
            LOADOPTIONALNOWARNINGDATA(io_toml,shape);
            LOADOPTIONALNOWARNINGDATA(io_toml,weight_param);
            KernelShape kernel_shape = static_cast<KernelShape>(shape.first);
            Shape kernel;
            switch (kernel_shape) {
                case KernelShape::CROSS:
                    kernel = shapes::Cross(rr.first[0],rr.first[1],rr.first[2]);
                    break;
                case KernelShape::ELLIPSOID:
                    kernel = shapes::Ellipsoid(rr.first[0],rr.first[1],rr.first[2]);
                    break;
                case KernelShape::CUBOID: {
                    kernel = shapes::CuboidR(rr.first[0],rr.first[1],rr.first[2]);
                    break;
                }
            }
            // Report postprocessing parameters
            cout<<"Postprocessing of the conductivity map\n";
            cout<<"Parameters:\n";
            cout<<"  Kernel size: ["<<rr.first[0]<<", "<<rr.first[1]<<", "<<rr.first[2]<<"]\n";
            cout<<"  Kernel shape: ("<<shape.first<<") "<<ToString(kernel_shape)<<"\n";
            cout<<"  Weight parameter: "<<weight_param.first<<"\n";
            cout<<"  Input address: "<<input<<"\n";
            cout<<"  Reference address: "<<reference<<"\n";
            cout<<"  Variance address: "<<variance<<"\n";
            cout<<"  Output address: "<<output<<"\n";
            cout<<endl;
            // Perform the postprocessing
            int postprocessing_error = PerformPostprocessing(output, input, reference, variance, kernel, weight_param.first, nn);
            if (postprocessing_error != 0) {
                cout<<"execution failed\n"<<endl;
                return 1;
            }
        }
    }
    if (postprocessing_output_epsr_addr.first!="") {
        string output = postprocessing_output_epsr_addr.first;
        string input;
        string reference;
        string variance;
        bool postprocess_flag = true;
        if (perform_only_postprocessing.first) {
            if (postprocessing_output_epsr_addr.first!="") {
                input = postprocessing_input_epsr_addr.first;
                reference = postprocessing_input_reference_addr.first;
                variance = postprocessing_input_epsr_variance_addr.first;
            } else {
                cout<<"WARNING: No input permittivity to be postprocessed!"<<endl;
                postprocess_flag = false;
            }
        } else {
            if (epsr_addr.first!="") {
                cfgdata<string> output_var_addr("","parameter.output-variance");
                LOADOPTIONALNOWARNINGDATA(io_toml, output_var_addr);
                input = epsr_addr.first;
                reference = refimg_addr.first;
                variance = output_var_addr.first + "/relative-permittivity";
            } else {
                cout<<"WARNING: No input permittivity to be postprocessed!"<<endl;
                postprocess_flag = false;
            }
        }
        if (postprocess_flag) {
            cfglist<int> rr({1,1,1},"postprocessing.kernel-size");
            cfgdata<int> shape(2,"postprocessing.kernel-shape");
            cfgdata<double> weight_param(0.10,"postprocessing.weight-param");
            LOADOPTIONALNOWARNINGLIST(io_toml,rr);
            LOADOPTIONALNOWARNINGDATA(io_toml,shape);
            LOADOPTIONALNOWARNINGDATA(io_toml,weight_param);
            KernelShape kernel_shape = static_cast<KernelShape>(shape.first);
            Shape kernel;
            switch (kernel_shape) {
                case KernelShape::CROSS:
                    kernel = shapes::Cross(rr.first[0],rr.first[1],rr.first[2]);
                    break;
                case KernelShape::ELLIPSOID:
                    kernel = shapes::Ellipsoid(rr.first[0],rr.first[1],rr.first[2]);
                    break;
                case KernelShape::CUBOID: {
                    kernel = shapes::CuboidR(rr.first[0],rr.first[1],rr.first[2]);
                    break;
                }
            }
            // Report postprocessing parameters
            cout<<"Postprocessing of the permittivity map\n";
            cout<<"Parameters:\n";
            cout<<"  Kernel size: ["<<rr.first[0]<<", "<<rr.first[1]<<", "<<rr.first[2]<<"]\n";
            cout<<"  Kernel shape: ("<<shape.first<<") "<<ToString(kernel_shape)<<"\n";
            cout<<"  Weight parameter: "<<weight_param.first<<"\n";
            cout<<"  Input address: "<<input<<"\n";
            cout<<"  Reference address: "<<reference<<"\n";
            cout<<"  Variance address: "<<variance<<"\n";
            cout<<"  Output address: "<<output<<"\n";
            cout<<endl;
            // Perform the postprocessing
            int postprocessing_error = PerformPostprocessing(output, input, reference, variance, kernel, weight_param.first, nn);
            if (postprocessing_error != 0) {
                cout<<"execution failed\n"<<endl;
                return 1;
            }
        }
    }
    //
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    cout<<"Execution ended in "<<elapsed.count()<<" s\n";
    cout<<endl;
    return 0;
}
