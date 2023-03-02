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

#ifndef IO_TOML_H_
#define IO_TOML_H_

#include <toml/toml.h>

#include <array>
#include <string>

#include "eptlib/util.h"
#include "eptlib/io/io_util.h"

namespace eptlib {

namespace io {

    /**
     * Class for interacting with .toml files.
     */
    class IOtoml {
        public:
            /**
             * Constructor.
             * 
             * @param fname address of the file to open.
             * @param mode file opening mode.
             */
            IOtoml(const std::string &fname, const Mode mode);

            /**
             * Destructor.
             */
            virtual ~IOtoml();

            /**
             * Extract a value from the TOML file content.
             * 
             * @tparam T typename of the value.
             * 
             * @param value pointer to the destination of the value.
             * @param uri address to the value in the toml file.
             * 
             * @return a Success, a MissingData if the value is not found, or
             *     a WrongDataFormat if the value is found but it is not of type `T'.
             */
            template <typename T>
            EPTlibError GetValue(T *value, const std::string &uri) const {
                const toml::Value* x = content_.find(uri);
                if (!x) {
                    return EPTlibError::MissingData;
                }
                if (!x->is<T>()) {
                    return EPTlibError::WrongDataFormat;
                }
                *value = x->as<T>();
                return EPTlibError::Success;
            }

            /**
             * Extract a `char' from the TOML file content.
             * 
             * @param value pointer to the destination of the value.
             * @param uri address to the value in the toml file.
             * 
             * @return a Success, a MissingData if the value is not found, or
             *     a WrongDataFormat if the value is found but it is not a `char'.
             */
            template <>
            EPTlibError GetValue<char>(char *value, const std::string &uri) const {
                std::string tmp;
                EPTlibError error = GetValue<std::string>(&tmp, uri);
                if (error==EPTlibError::Success) {
                    *value = tmp[0];
                }
                return error;
            }

            /**
             * Extract an array from the TOML file content.
             * 
             * @tparam T typename of the values in the array.
             * 
             * @param array pointer to the destination of the array.
             * @param uri address to the array in the toml file.
             * 
             * @return a Success, a MissingData if the array is not found, or
             *     a WrongDataFormat if the array is found but its values are not of type `T'.
             * 
             * @deprecated
             */
            template <typename T>
            EPTlibError GetArrayOf(std::array<T,N_DIM> *array, const std::string &uri) const {
                toml::Array a;
                EPTlibError error = GetValue<toml::Array>(&a, uri);
                if (error==EPTlibError::Success) {
                    if (!a[0].is<T>()) {
                        return EPTlibError::WrongDataFormat;
                    }
                    for (int d = 0; d<N_DIM; ++d) {
                        (*array)[d] = a[d].as<T>();
                    }
                }
                return error;
            }
        private:
            /// TOML file.
            std::ifstream file_;
            /// TOML file content.
            toml::Value content_;
    };

}  // io

}  // eptlib

#endif  // IO_TOML_H_
