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

#include "gtest/gtest.h"

#include "eptlib/io/io_hdf5.h"

#include <string>

TEST(IOHdf5GTest,GetFileName) {
    std::string ifname = "test/input/test_input.h5";
    eptlib::io::IOh5 ifile(ifname, eptlib::io::Mode::In);
    ASSERT_STREQ(ifile.GetFileName().c_str(), ifname.c_str());
    std::string ofname = "test/input/test_output.h5";
    eptlib::io::IOh5 ofile(ofname, eptlib::io::Mode::Out);
    ASSERT_STREQ(ofile.GetFileName().c_str(), ofname.c_str());
}

TEST(IOHdf5GTest,GetFileSize) {
    std::string ifname = "test/input/test_input.h5";
    eptlib::io::IOh5 ifile(ifname, eptlib::io::Mode::In);
    ASSERT_EQ(ifile.GetFileSize(), 51464);
}

TEST(IOHdf5GTest,GetFileDescription) {
    std::string ifname = "test/input/test_input.h5";
    eptlib::io::IOh5 ifile(ifname, eptlib::io::Mode::In);
    ASSERT_EQ(ifile.GetFileDescription(), (ifname + " [50 KiB]").c_str());
}

TEST(IOHdf5GTest,ReadDataset) {
    std::string fname = "test/input/test_input.h5";
    const eptlib::io::IOh5 ifile(fname, eptlib::io::Mode::In);
    std::string url = "/test/input/";
    std::string urn = "data";
    eptlib::Image<double> img;
    ifile.ReadDataset(&img, url,urn);
    ASSERT_EQ(img.GetSize(0),20);
    ASSERT_EQ(img.GetSize(1),10);
    ASSERT_EQ(img.GetSize(2),30);
    for (int idx = 0; idx<img.GetNVox(); ++idx) {
        ASSERT_DOUBLE_EQ(img(idx), idx);
    }
}

TEST(IOHdf5GTest,WriteDataset) {
    std::string fname = "test/input/test_output.h5";
    const eptlib::io::IOh5 ofile(fname, eptlib::io::Mode::Out);
    std::string url = "/test/output/";
    std::string urn = "data";
    eptlib::Image<double> img_out(20,10,30);
    std::iota(img_out.GetData().begin(),img_out.GetData().end(),0);
    ofile.WriteDataset(img_out, url,urn);
    eptlib::Image<double> img_in;
    ofile.ReadDataset(&img_in, url,urn);
    ASSERT_EQ(img_in.GetSize(0),20);
    ASSERT_EQ(img_in.GetSize(1),10);
    ASSERT_EQ(img_in.GetSize(2),30);
    for (int idx = 0; idx<img_in.GetNVox(); ++idx) {
        ASSERT_DOUBLE_EQ(img_in(idx), idx);
    }
}
