/*****************************************************************************
*
*     Program: EPTlib
*     Author: Alessandro Arduino <a.arduino@inrim.it>
*
*  MIT License
*
*  Copyright (c) 2020-2022  Alessandro Arduino
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

#include "eptlib/io/io_toml.h"

#include <fstream>
#include <ios>

using namespace eptlib;
using namespace eptlib::io;

// IOtoml constructor
IOtoml::
IOtoml(const std::string &fname, const Mode mode) :
    fname_(fname), mode_(mode) {
    switch (mode_) {
        case Mode::In:
            file_.open(fname_);
            if (file_.is_open()) {
                toml::ParseResult parsed(toml::parse(file_));
                if (parsed.valid()) {
                    content_ = parsed.value;
                } else {
                    throw std::ios_base::failure(parsed.errorReason);
                }
            } else {
                throw std::ios_base::failure("Impossible to open file '"+fname_+"'");
            }
            break;
        case Mode::Out:
        case Mode::Append:
            throw std::ios_base::failure("TOML files are handled only in input");
            break;
    }
    return;
}

// IOtoml destructor
IOtoml::
~IOtoml() {
    file_.close();
    return;
}
