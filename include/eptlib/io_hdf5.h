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

#ifndef IO_HDF5_H_
#define IO_HDF5_H_

#include <H5Cpp.h>

#include <algorithm>
#include <array>
#include <vector>
#include <string>

#include "util.h"

namespace eptlib {

namespace io {

    /**
     * File opening modes.
     */
    typedef enum Mode {
        /// Open the file in read-only mode.
        In = 0,
        /// Create a new file or overwrite an existing file.
        Out,
        /// Create a new file or append to an existing file.
        Append,
    } Mode_t;

    /**
     * IO states.
     */
    typedef enum State {
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
    } State_t;
    /**
     * Translates in a human-readable string the input IO state.
     * 
     * @param state is a State symbol.
     * @return the human-readable description of the IO state.
     */
    const std::string ToString(const State_t state);
    
    /**
     * Class for interacting with .h5 files.
     */
    class IOh5 {
        public:
            /**
             * Constructor.
             * 
             * @param fname address of the file to open.
             * @param mode file opening mode.
             */
            IOh5(const std::string &fname, const Mode_t mode);
            /**
             * Destructor.
             */
            ~IOh5();
            /**
             * Read a dataset from the .h5 file.
             * 
             * @tparam T scalar typename.
             * 
             * @param data pointer to the data begin.
             * @param nn number of elements in each direction of the data.
             * @param url url of the dataset.
             * @param urn urn of the dataset.
             * 
             * @return the IO state.
             */
            template <typename T>
            State_t ReadDataset(std::vector<T> &data, std::array<int,NDIM> &nn, const std::string &url, const std::string &urn);
            /**
             * Write a datates into the .h5 file.
             * 
             * @tparam T scalar typename.
             * 
             * @param data pointer to the data begin.
             * @param nn number of elements in each direction of the data.
             * @param url url of the dataset.
             * @param urn urn of the dataset.
             * 
             * @return the IO state.
             */
            template <typename T>
            State_t WriteDataset(std::vector<T> &data, std::array<int,NDIM> &nn, const std::string &url, const std::string &urn);
        private:
            /// Address of the file to open.
            std::string fname_;
            /// File opening mode.
            Mode_t mode_;
            /// HDF5 file.
            H5::H5File file_;
    };

    // ---------------------------------------------------------------------------
    // -------------------------  Implementation detail  -------------------------
    // ---------------------------------------------------------------------------

    // Translates in a human-readable string the input IO state
    inline const std::string ToString(const State_t state) {
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

}  // io

}  // eptlib

#endif  // IO_HDF5_H_
