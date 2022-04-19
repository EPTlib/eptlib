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

#ifndef IO_HDF5_H_
#define IO_HDF5_H_

#include <H5Cpp.h>

#include <string>

#include "eptlib/image.h"
#include "eptlib/io/io_util.h"

namespace eptlib {

namespace io {
    
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
            IOh5(const std::string &fname, const Mode mode);
            /**
             * Destructor.
             */
            ~IOh5();
            /**
             * Read a dataset from the .h5 file.
             * 
             * @tparam T scalar typename.
             * 
             * @param img pointer to the destination image.
             * @param url url of the dataset.
             * @param urn urn of the dataset.
             * 
             * @return the IO state.
             */
            template <typename T>
            State ReadDataset(Image<T> *img, const std::string &url, const std::string &urn);
            /**
             * Write a datates into the .h5 file.
             * 
             * @tparam T scalar typename.
             * 
             * @param img source image.
             * @param url url of the dataset.
             * @param urn urn of the dataset.
             * 
             * @return the IO state.
             */
            template <typename T>
            State WriteDataset(const Image<T> &img, const std::string &url, const std::string &urn) const;
        private:
            /// Address of the file to open.
            std::string fname_;
            /// File opening mode.
            Mode mode_;
            /// HDF5 file.
            H5::H5File file_;
    };

}  // namespace io

}  // namespace eptlib

#endif  // IO_HDF5_H_
