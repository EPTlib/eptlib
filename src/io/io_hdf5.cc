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

#include "eptlib/io/io_hdf5.h"

#include <algorithm>

using namespace eptlib;
using namespace eptlib::io;

namespace {

    /**
     * Provide the uri, given url and urn.
     * 
     * @param url url
     * @param urn urn
     * 
     * @return uri
     */
    inline std::string URI(const std::string &url, const std::string &urn) {
        std::string uri("/"+url+"/"+urn);
        size_t pos = uri.find("//");
        while (pos != std::string::npos) {
            uri.replace(pos, 2, "/");
            pos = uri.find("//");
        }
        return uri;
    }

    /**
     * Create an hdf5 group, given an url.
     * 
     * @param file hdf5 file
     * @param url url
     * 
     * @return hdf5 group
     */
    inline H5::Group CreateGroup(const H5::H5File &file, const std::string &url) {
        size_t depth = 0;
        size_t snip = url.find_first_of("/", 0) ? 0 : 1;
        size_t snap = url.find_first_of("/", snip);
        std::string subpath = url.substr(snip,snap);
        H5::Group* tmp, group;
        while (!subpath.empty()) {
            try {
                group = depth ? tmp->openGroup(subpath) : file.openGroup(subpath);
            } catch (const H5::Exception& e) {
                group = depth ? tmp->createGroup(subpath) : file.createGroup(subpath);
            }
            snip = ++snap;
            snap = url.find_first_of("/",snip);
            subpath = url.substr(snip,snap);
            tmp = &group;
            depth++;
        }
        return *tmp;
    }

    // HDF5 types traits
    template <typename T>
    struct HDF5Types;
    // traits specialisations
    template <>
    struct HDF5Types<size_t> {
        static const H5::DataType Type() {
            return H5::PredType::NATIVE_ULONG;
        }
    };
    template<>
    struct HDF5Types<double> {
        static const H5::DataType Type() {
            return H5::PredType::NATIVE_DOUBLE;
        }
    };
    template<>
    struct HDF5Types<float> {
        static const H5::DataType Type() {
            return H5::PredType::NATIVE_FLOAT;
        }
    };
    template<>
    struct HDF5Types<int> {
        static const H5::DataType Type()  {
            return H5::PredType::NATIVE_INT;
        }
    };
    template<>
    struct HDF5Types<long> {
        static const H5::DataType Type() {
            return H5::PredType::NATIVE_LONG;
        }
    };

}  //

// IOh5 constructor
IOh5::
IOh5(const std::string &fname, const Mode_t mode) :
    fname_(fname), mode_(mode) {
    H5::Exception::dontPrint();
    switch (mode_) {
        case Mode::In:
            file_ = H5::H5File(fname_, H5F_ACC_RDONLY);
            break;
        case Mode::Out:
            file_ = H5::H5File(fname_, H5F_ACC_TRUNC);
            break;
        case Mode::Append:
            try {
                file_ = H5::H5File(fname_, H5F_ACC_RDWR);
            } catch (const H5::FileIException &e) {
                file_ = H5::H5File(fname_, H5F_ACC_TRUNC);
            }
            break;
    }
    return;
}

// IOh5 destructor
IOh5::
~IOh5() {
    file_.close();
    return;
}

// IOh5 read dataset
template <typename T>
State_t IOh5::
ReadDataset(std::vector<T> &data, std::array<int,NDIM> &nn, const std::string &url, const std::string &urn) {
    H5::Exception::dontPrint();
    try {
        // locate the dataset
        H5::DataSet dset = file_.openDataSet(URI(url,urn));
        H5::DataSpace dspace = dset.getSpace();
        // read the data
        std::array<hsize_t,NDIM> dims;
        size_t ndim = dspace.getSimpleExtentDims(dims.data(),NULL);
        std::reverse_copy(dims.begin(),dims.end(),nn.begin());
        data.resize(Prod(nn));
        dset.read(data.data(),::HDF5Types<T>::Type());
    } catch (const H5::FileIException& e) {
        return State::HDF5FileException;
    } catch (const H5::DataSetIException& e) {
        return State::HDF5DatasetException;
    } catch (const H5::DataSpaceIException& e) {
        return State::HDF5DataspaceException;
    } catch (const H5::DataTypeIException& e) {
        return State::HDF5DatatypeException;
    }
    return State::Success;
}

// IOh5 write dataset
template <typename T>
State_t IOh5::
WriteDataset(std::vector<T> &data, std::array<int,NDIM> &nn, const std::string &url, const std::string &urn) {
    H5::Exception::dontPrint();
    try {
        // open or create the group
        H5::Group group;
        try {
            group = file_.openGroup(url);
        } catch (const H5::Exception& e) {
            group = CreateGroup(file_,url);
        }
        // create the dataset
        std::array<hsize_t,NDIM> dims;
        std::reverse_copy(nn.begin(),nn.end(),dims.begin());
        H5::DataSpace dspace(NDIM,dims.data());
        H5::DataType dtype(::HDF5Types<T>::Type());
        H5::DataSet dset = group.createDataSet(urn,dtype,dspace);
        // write the data in the dataset
        dset.write(data.data(),dtype);
    } catch (const H5::FileIException& e) {
        return State::HDF5FileException;
    } catch (const H5::GroupIException& e) {
        return State::HDF5FileException;
    } catch (const H5::DataSetIException& e) {
        return State::HDF5DatasetException;
    } catch (const H5::DataSpaceIException& e) {
        return State::HDF5DataspaceException;
    } catch (const H5::DataTypeIException& e) {
        return State::HDF5DatatypeException;
    }
    return State::Success;
}

// Template specialisations
// ReadDataset
template State_t IOh5::ReadDataset<size_t>(std::vector<size_t> &data, std::array<int,NDIM> &nn, const std::string &url, const std::string &urn);
template State_t IOh5::ReadDataset<float>(std::vector<float> &data, std::array<int,NDIM> &nn, const std::string &url, const std::string &urn);
template State_t IOh5::ReadDataset<double>(std::vector<double> &data, std::array<int,NDIM> &nn, const std::string &url, const std::string &urn);
template State_t IOh5::ReadDataset<int>(std::vector<int> &data, std::array<int,NDIM> &nn, const std::string &url, const std::string &urn);
template State_t IOh5::ReadDataset<long>(std::vector<long> &data, std::array<int,NDIM> &nn, const std::string &url, const std::string &urn);
// WriteDataset
template State_t IOh5::WriteDataset<size_t>(std::vector<size_t> &data, std::array<int,NDIM> &nn, const std::string &url, const std::string &urn);
template State_t IOh5::WriteDataset<float>(std::vector<float> &data, std::array<int,NDIM> &nn, const std::string &url, const std::string &urn);
template State_t IOh5::WriteDataset<double>(std::vector<double> &data, std::array<int,NDIM> &nn, const std::string &url, const std::string &urn);
template State_t IOh5::WriteDataset<int>(std::vector<int> &data, std::array<int,NDIM> &nn, const std::string &url, const std::string &urn);
template State_t IOh5::WriteDataset<long>(std::vector<long> &data, std::array<int,NDIM> &nn, const std::string &url, const std::string &urn);
