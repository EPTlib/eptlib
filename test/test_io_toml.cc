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

#include "eptlib/io/io_toml.h"

#include "eptlib/util.h"

using namespace eptlib;
using namespace eptlib::io;

TEST(IOtomlGTest,ReadFile) {
    std::string string;
    int integer;
    double floating;
    std::array<int,NDIM> array;
    std::string var;

    std::string fname = "test/input/test_input.toml";
    IOtoml ifile(fname, Mode::In);

    ifile.GetValue<std::string>(string, "string");
    ifile.GetValue<int>(integer,"integer");
    ifile.GetValue<double>(floating,"floating");
    ifile.GetArrayOf<int>(array,"array");
    ifile.GetValue<std::string>(var,"group.var");

    ASSERT_STREQ(string.c_str(),"string");
    ASSERT_EQ(integer,5);
    ASSERT_DOUBLE_EQ(floating,0.5);
    for (int d = 0; d<NDIM; ++d) {
        ASSERT_EQ(array[d],d);
    }
    ASSERT_STREQ(var.c_str(),"var in group");
}
