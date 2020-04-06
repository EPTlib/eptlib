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

#include "eptlib/io_hdf5.h"

#include "eptlib/util.h"

using namespace eptlib;
using namespace eptlib::io;

TEST(IOhdf5GTest,ReadDataset) {
    std::array<int,NDIM> nn;
    std::vector<double> data;

    std::string fname = "test/input/test_input.h5";
    IOh5 ifile(fname, Mode::In);

    std::string url = "/test/input/";
    std::string urn = "data";
    ifile.ReadDataset(data,nn, url,urn);

    std::array<int,NDIM> nn_expected{20,10,30};
    std::vector<double> data_expected(Prod(nn_expected));
    std::iota(data_expected.begin(),data_expected.end(),0.0);
    
    for (int d = 0; d<NDIM; ++d) {
        ASSERT_EQ(nn[d],nn_expected[d]);
    }
    for (int idx = 0; idx<Prod(nn); ++idx) {
        ASSERT_DOUBLE_EQ(data[idx],data_expected[idx]);
    }
}

TEST(IOhdf5GTest,WriteDataset) {
    std::array<int,NDIM> nn_expected{20,10,30};
    std::vector<double> data_expected(Prod(nn_expected));
    std::iota(data_expected.begin(),data_expected.end(),0.0);

    std::string fname = "test/input/test_output.h5";
    IOh5 ofile(fname, Mode::Out);

    std::string url = "/test/input/";
    std::string urn = "data";
    ofile.WriteDataset(data_expected,nn_expected, url,urn);

    std::array<int,NDIM> nn;
    std::vector<double> data;

    ofile.ReadDataset(data,nn, url,urn);
    for (int d = 0; d<NDIM; ++d) {
        ASSERT_EQ(nn[d],nn_expected[d]);
    }
    for (int idx = 0; idx<Prod(nn); ++idx) {
        ASSERT_DOUBLE_EQ(data[idx],data_expected[idx]);
    }
}
