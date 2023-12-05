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

#include "eptlib/io/io_hdf5.h"

namespace {

    // HDF5 types traits
    template <typename T> struct HDF5Types;
    // traits specialisations
    template <> struct HDF5Types<float> {
        static const H5::DataType Type() {
            return H5::PredType::NATIVE_FLOAT;
        }
    };
    template <> struct HDF5Types<double> {
        static const H5::DataType Type() {
            return H5::PredType::NATIVE_DOUBLE;
        }
    };
    template <> struct HDF5Types<int> {
        static const H5::DataType Type()  {
            return H5::PredType::NATIVE_INT;
        }
    };
    template <> struct HDF5Types<long> {
        static const H5::DataType Type() {
            return H5::PredType::NATIVE_LONG;
        }
    };

}

// IOh5 constructor
eptlib::io::IOh5::
IOh5(const std::string &fname, const eptlib::io::Mode mode) {
    H5::Exception::dontPrint();
    switch (mode) {
        case eptlib::io::Mode::In:
            file_ = H5::H5File(fname, H5F_ACC_RDONLY);
            break;
        case eptlib::io::Mode::Out:
            file_ = H5::H5File(fname, H5F_ACC_TRUNC);
            break;
        case eptlib::io::Mode::Append:
            try {
                file_ = H5::H5File(fname, H5F_ACC_RDWR);
            } catch (const H5::FileIException&) {
                file_ = H5::H5File(fname, H5F_ACC_TRUNC);
            }
            break;
    }
    return;
}

// IOh5 destructor
eptlib::io::IOh5::
~IOh5() {
    file_.close();
    return;
}

// IOh5 get the address of the file
std::string eptlib::io::IOh5::
IOh5::GetFileName() const {
    return file_.getFileName();
}

// IOh5 get the size of the file in byte
size_t eptlib::io::IOh5::
IOh5::GetFileSize() const {
    return file_.getFileSize();
}

// IOh5 get a description of the file with address and size
std::string eptlib::io::IOh5::
IOh5::GetFileDescription() const {
    std::string fname = this->GetFileName();
    size_t fsize = this->GetFileSize();
    return fname+" ["+BytesWithSuffix(fsize)+"]";
}

// IOh5 provide the uri, given url and urn
std::string eptlib::io::IOh5::
URI(const std::string &url, const std::string &urn) const {
    std::string uri("/"+url+"/"+urn);
    size_t pos = uri.find("//");
    while (pos != std::string::npos) {
        uri.replace(pos, 2, "/");
        pos = uri.find("//");
    }
    return uri;
}

// IOh5 create an hdf5 group, given an url
H5::Group eptlib::io::IOh5::
CreateGroup(const std::string &url) const {
    size_t depth = 0;
    size_t snip = url.find_first_of("/", 0) ? 0 : 1;
    size_t snap = url.find_first_of("/", snip);
    std::string subpath = url.substr(snip,snap-snip);
    H5::Group group;
    while (!subpath.empty()) {
        try {
            group = depth ? group.openGroup(subpath) : file_.openGroup(subpath);
        } catch (const H5::Exception&) {
            group = depth ? group.createGroup(subpath) : file_.createGroup(subpath);
        }
        snip = ++snap;
        snap = url.find_first_of("/",snip);
        subpath = url.substr(snip,snap-snip);
        depth++;
    }
    return group;
}

//IOh5 read a dataset from the .h5 file
template <typename T>
eptlib::io::State eptlib::io::IOh5::
ReadDataset(eptlib::Image<T> *img, const std::string &url, const std::string &urn) const {
    H5::Exception::dontPrint();
    try {
        // locate the dataset
        H5::DataSet dset = file_.openDataSet(URI(url,urn));
        H5::DataSpace dspace = dset.getSpace();
        // read the image size
        std::array<hsize_t,N_DIM> nn;
        auto hdim = dspace.getSimpleExtentDims(nn.data());
        std::reverse(nn.begin(),nn.begin()+hdim);
        for (auto d=hdim; d<N_DIM; ++d) {
            nn[d] = 1;
        }
        // read the image data
        *img = Image<T>(nn[0],nn[1],nn[2]);
        dset.read(img->GetData().data(),HDF5Types<T>::Type());
    } catch (const H5::FileIException&) {
        return State::HDF5FileException;
    } catch (const H5::DataSetIException&) {
        return State::HDF5DatasetException;
    } catch (const H5::DataSpaceIException&) {
        return State::HDF5DataspaceException;
    } catch (const H5::DataTypeIException&) {
        return State::HDF5DatatypeException;
    }
    return State::Success;
}

template <typename T>
eptlib::io::State eptlib::io::IOh5::
WriteDataset(const eptlib::Image<T> &img, const std::string &url, const std::string &urn) const {
    H5::Exception::dontPrint();
    try {
        std::string uri = URI(url,urn);
        size_t snip = uri.find_last_of("/")+1;
        std::string group_names = uri.substr(0,snip);
        std::string dataset_name = uri.substr(snip);
        // open or create the group
        H5::Group group;
        try {
            group = file_.openGroup(group_names);
        } catch (const H5::Exception&) {
            group = CreateGroup(group_names);
        }
        // create the dataset
        std::array<hsize_t,N_DIM> dims(img.GetSize());
        std::reverse(dims.begin(),dims.end());
        H5::DataSpace dspace(N_DIM,dims.data());
        H5::DataType dtype(HDF5Types<T>::Type());
        H5::DataSet dset;
        try {
            dset = group.createDataSet(dataset_name,dtype,dspace);
        } catch (const H5::Exception&) {
            dset = group.openDataSet(dataset_name);
        }
        // write the data in the dataset
        dset.write(img.GetData().data(),dtype,dspace);
    } catch (const H5::FileIException&) {
        return State::HDF5FileException;
    } catch (const H5::GroupIException&) {
        return State::HDF5FileException;
    } catch (const H5::DataSetIException&) {
        return State::HDF5DatasetException;
    } catch (const H5::DataSpaceIException&) {
        return State::HDF5DataspaceException;
    } catch (const H5::DataTypeIException&) {
        return State::HDF5DatatypeException;
    }
    return State::Success;
}

template eptlib::io::State eptlib::io::IOh5::ReadDataset(Image<float>*,const std::string&,const std::string&) const;
template eptlib::io::State eptlib::io::IOh5::ReadDataset(Image<double>*,const std::string&,const std::string&) const;
template eptlib::io::State eptlib::io::IOh5::ReadDataset(Image<int>*,const std::string&,const std::string&) const;
template eptlib::io::State eptlib::io::IOh5::ReadDataset(Image<long>*,const std::string&,const std::string&) const;

template eptlib::io::State eptlib::io::IOh5::WriteDataset(const Image<float>&,const std::string&,const std::string&) const;
template eptlib::io::State eptlib::io::IOh5::WriteDataset(const Image<double>&,const std::string&,const std::string&) const;
template eptlib::io::State eptlib::io::IOh5::WriteDataset(const Image<int>&,const std::string&,const std::string&) const;
template eptlib::io::State eptlib::io::IOh5::WriteDataset(const Image<long>&,const std::string&,const std::string&) const;
