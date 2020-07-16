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
#include <utility>

#include "eptlib/ept_methods.h"
#include "eptlib/io/io_hdf5.h"
#include "eptlib/io/io_toml.h"
#include "eptlib/version.h"

#include "eptlib_main.h"

#define LOADMANDATORY(what,io_toml,data,T) { \
    EPTlibError_t error = io_toml->what<T>(data.first,data.second); \
    if (error!=EPTlibError::Success) { \
        string message = ToString(error)+" '"+data.second+"'"; \
        cout<<"FATAL ERROR in config file: "<<message<<endl; \
        delete io_toml; \
        return 1; \
    } \
}
#define LOADMANDATORYDATA(io_toml,data) LOADMANDATORY(GetValue,io_toml,data,decltype(data.first))
#define LOADMANDATORYLIST(io_toml,data) LOADMANDATORY(GetArrayOf,io_toml,data,decltype(data.first)::value_type)

#define LOADOPTIONAL(what,io_toml,data,T) { \
    EPTlibError_t error = io_toml->what<T>(data.first,data.second); \
    if (WALL && error!=EPTlibError::Success) { \
        string message = ToString(error)+" '"+data.second+"'"; \
        cout<<"WARNING in config file: "<<message<<endl; \
    } \
}
#define LOADOPTIONALDATA(io_toml,data) LOADOPTIONAL(GetValue,io_toml,data,decltype(data.first))
#define LOADOPTIONALLIST(io_toml,data) LOADOPTIONAL(GetArrayOf,io_toml,data,decltype(data.first[0]))

using namespace std;
using namespace eptlib;

constexpr char software::name[];

template <class T> using data = pair<T,string>;
template <class T> using list = pair<array<T,NDIM>,string>;

int main(int argc, char **argv) {
    // opening boilerplate
    cout<<project::str()<<" ("<<build::str()<<") ["<<compiler::str()<<"]\n"<<endl;
    cout<<LicenseBoilerplate()<<endl;
    // check the number of input
    if (argc<2) {
        cout<<"Usage example: "<<software::str()<<" <config file>"<<endl;
        return -1;
    }
    // load the config file
    io::IOtoml *io_toml;
    try {
        io_toml = new io::IOtoml(std::string(argv[1]),io::Mode::In);
    } catch(const ios_base::failure &e) {
        cout<<"FATAL ERROR in config file: "<<e.what()<<endl;
        delete io_toml;
        return 1;
    }
    // declare the input variables
    //   mandatory input
    data<string> title; title.second = "title";
    data<string> descr; descr.second = "description";
    data<int> method; method.second = "method";
    list<int> nn; nn.second = "mesh.size";
    list<double> dd; dd.second = "mesh.step";
    data<double> freq; freq.second = "input.frequency";
    //   optional input
    data<int> n_txch; n_txch.first = 1; n_txch.second = "input.tx-channels";
    data<int> n_rxch; n_rxch.first = 1; n_rxch.second = "input.rx-channels";
    data<string> txsens_addr; txsens_addr.first = ""; txsens_addr.second = "input.tx-sensitivity";
    data<string> trxphase_addr; trxphase_addr.first = ""; trxphase_addr.second = "input.trx-phase";
    data<string> refimg_addr; refimg_addr.first = ""; refimg_addr.second = "input.reference-img";
    data<string> sigma_addr; sigma_addr.first = ""; sigma_addr.second = "output.electric-conductivity";
    data<string> epsr_addr; epsr_addr.first = ""; epsr_addr.second = "output.relative-permittivity";
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
    // Check the EPT method
    EPTMethod_t ept_method = static_cast<EPTMethod_t>(method.first);
    if (ept_method<0||ept_method>=EPTMethod::END) {
        cout<<"FATAL ERROR in config file: Wrong data format '"<<method.second<<"'"<<endl;
        delete io_toml;
        return 1;
    }
    // report the readen values
    //   head
    cout<<"\n  "<<title.first<<"\n";
    cout<<"  "<<descr.first<<"\n";
    //   body
    cout<<"\n  Method: "<<ToString(ept_method)<<" ("<<method.first<<")\n";
    cout<<"\n  Mesh size: ["<<nn.first[0]<<", "<<nn.first[1]<<", "<<nn.first[2]<<"]\n";
    cout<<"  Mesh step: ["<<dd.first[0]<<", "<<dd.first[1]<<", "<<dd.first[2]<<"]\n";
    cout<<"\n  Frequency: "<<freq.first<<"\n";
    cout<<"  Tx channels: "<<n_txch.first<<"\n";
    cout<<"  Rx channels: "<<n_rxch.first<<"\n";
    //   tail

    // free the dynamically allocated memory
    delete io_toml;
    return 0;
}
