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

#ifndef IO_UTIL_H_
#define IO_UTIL_H_

#include "eptlib/util.h"

namespace eptlib {

namespace io {

    /**
     * File opening modes.
     */
    enum class Mode {
        /// Open the file in read-only mode.
        In = 0,
        /// Create a new file or overwrite an existing file.
        Out,
        /// Create a new file or append to an existing file.
        Append,
    };

    /**
     * IO states.
     */
    enum class State {
        /// Success.
        Success = 0,
        /// HDF5 file error.
        HDF5FileException,
        /// HDF5 dataset error.
        HDF5DatasetException,
        /// HDF5 dataspace error.
        HDF5DataspaceException,
        /// HDF5 datatype error.
        HDF5DatatypeException,
    };
    /**
     * Translates in a human-readable string the input IO state.
     * 
     * @param state is a State symbol.
     * @return the human-readable description of the IO state.
     */
    const std::string ToString(const State state);

    /**
     * Deduce the filename and the uri from a given address.
     * 
     * @param address complete address of filename and uri separated by the
     *      character ':'.
     * @param fname address of the file within the complete address.
     * @param uri address of the element within the complete address.
     */
    void GetAddress(const std::string &address, std::string &fname, std::string &uri);
    
    // ---------------------------------------------------------------------------
    // -------------------------  Implementation detail  -------------------------
    // ---------------------------------------------------------------------------

    // Translates in a human-readable string the input IO state
    inline const std::string ToString(const State state) {
        switch (state) {
            case State::Success:
                return "Success";
            case State::HDF5FileException:
                return "IO Error: HDF5, file exception";
            case State::HDF5DatasetException:
                return "IO Error: HDF5, dataset exception";
            case State::HDF5DataspaceException:
                return "IO Error: HDF5, dataspace exception";
            case State::HDF5DatatypeException:
                return "IO Error: HDF5, datatype exception";
        }
        return "";
    }

    // Read filename and uri from the address
    inline void GetAddress(const std::string &address, std::string &fname, std::string &uri) {
        size_t snip = 0;
        size_t snap = address.find_first_of(":",0);
        fname = address.substr(snip,snap);
        snip = ++snap;
        uri = address.substr(snip);
        return;
    }

}  // io

}  // eptlib

#endif  // IO_UTIL_H_
