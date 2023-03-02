/*****************************************************************************
*
*     Program: EPTlib
*     Author: Alessandro Arduino <a.arduino@inrim.it>
*
*  MIT License
*
*  Copyright (c) 2023  Alessandro Arduino
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

#include "eptlib/io/io_util.h"

#include <string>

TEST(IOUtilGTest,GetAddress) {
    std::string address = "/test/input/test_input.xx:/group/dataset";
    std::string fname;
    std::string uri;
    eptlib::io::GetAddress(address, fname, uri);
    ASSERT_STREQ(fname.c_str(), "/test/input/test_input.xx");
    ASSERT_STREQ(uri.c_str(), "/group/dataset");
}

TEST(IOUtilGTest,BytesWithSuffix) {
    size_t     byte = ((size_t) 1)<<0;
    size_t kibibyte = ((size_t) 2)<<10;
    size_t mebibyte = ((size_t) 3)<<20;
    size_t gibibyte = ((size_t) 4)<<30;
    size_t tebibyte = ((size_t) 5)<<40;
    size_t pebibyte = ((size_t) 6)<<50;
    ASSERT_STREQ(eptlib::io::BytesWithSuffix(    byte).c_str(), "1 bytes");
    ASSERT_STREQ(eptlib::io::BytesWithSuffix(kibibyte).c_str(), "2 KiB");
    ASSERT_STREQ(eptlib::io::BytesWithSuffix(mebibyte).c_str(), "3 MiB");
    ASSERT_STREQ(eptlib::io::BytesWithSuffix(gibibyte).c_str(), "4 GiB");
    ASSERT_STREQ(eptlib::io::BytesWithSuffix(tebibyte).c_str(), "5 TiB");
    ASSERT_STREQ(eptlib::io::BytesWithSuffix(pebibyte).c_str(), "6144 TiB");
}
