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

#ifndef EPTLIB_UTIL_H_
#define EPTLIB_UTIL_H_

#include <algorithm>
#include <complex>
#include <functional>
#include <numeric>
#include <string>
#include <type_traits>

namespace eptlib {

    /// Number of spatial dimensions.
    constexpr int N_DIM = 3;
    /// Pi.
    constexpr double PI = 3.14159265358979323846;
    /// Speed of light [m/s].
    constexpr double C0 = 299792458.0;
    /// Vacuum permeability [H/m].
    constexpr double MU0 = 4.0e-7*PI;
    /// Vacuum permittivity [F/m].
    constexpr double EPS0 = 1.0/MU0/C0/C0;

    /**
     * Error codes that can be provided by EPTlib functions and methods.
     */
    enum class EPTlibError {
        /// Success.
        Success = 0,
        /// Missing data error.
        MissingData,
        /// Out of range error.
        OutOfRange,
        /// Wrong data format error.
        WrongDataFormat,
        /// Unknown error.
        Unknown,
    };
    /**
     * Translates in a human-readable string the input EPTlibError symbol.
     *
     * @param error is an EPTlibError symbol.
     * @return the human-readable description of the input EPTlibError symbol.
     */
    inline const std::string ToString(const EPTlibError error) {
        switch (error) {
            case EPTlibError::Success:
                return "Success";
            case EPTlibError::MissingData:
                return "Missing data";
            case EPTlibError::OutOfRange:
                return "Out of range";
            case EPTlibError::WrongDataFormat:
                return "Wrong data format";
            case EPTlibError::Unknown:
                return "Unknown error";
        }
        return "";
    };

    /**
     * Return the license boilerplate as a string.
     * 
     * @return the license boilerplate.
     */
    inline const std::string LicenseBoilerplate() {
        const std::string boilerplate = "MIT License\n"
            "Copyright (c) 2020-2023  Alessandro Arduino\n"
            "Istituto Nazionale di Ricerca Metrologica (INRiM)\n"
            "\n"
            "THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR\n"
            "IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,\n"
            "FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE\n"
            "AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER\n"
            "LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,\n"
            "OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE\n"
            "SOFTWARE.\n";
        return boilerplate;
    };

    /**
     * Define at compile-time the return type of an unambiguous function.
     * 
     * @tparam Function the function typename.
     */
    template <typename Function>
    using FunctionReturnType = typename decltype(std::function{std::declval<Function>()})::result_type;

    /**
     * Compute the sum of all the elements in a container.
     * 
     * @tparam T container typename.
     * 
     * @param v container of elements.
     * 
     * @return the sum of all the elements in `v'.
     */
    template <typename T>
    inline typename T::value_type Sum(const T &v) {
        using type = typename T::value_type;
        return std::accumulate(v.begin(),v.end(),static_cast<type>(0),std::plus<type>());
    }

    /**
     * Compute the products of all the elements in a container.
     * 
     * @tparam T container typename.
     * 
     * @param v container of elements.
     * 
     * @return the products of all the elements in `v'.
     */
    template <typename T>
    inline typename T::value_type Prod(const T &v) {
        using type = typename T::value_type;
        return std::accumulate(v.begin(),v.end(),static_cast<type>(1),std::multiplies<type>());
    }

    /**
     * Compute the arithmetic mean of all the elements in a container.
     * 
     * @tparam T container typename.
     * 
     * @param v container of elements.
     * 
     * @return the arithmetic mean of all the elements in `v'.
     * 
     * Works only if the container elements are floating points.
     */
    template <typename T, std::enable_if_t<std::is_floating_point_v<typename T::value_type>, bool> = true>
    inline typename T::value_type ArithmeticMean(const T &v) {
        using type = typename T::value_type;
        return std::accumulate(v.begin(),v.end(),static_cast<type>(0),
            [&] (const type &partial_sum, const type &new_element) -> type {
                return partial_sum + new_element/static_cast<type>(v.size());
            }
        );
    }

    /**
     * Compute the algebraic maximum of all the elements in a containes.
     * 
     * @tparam T container typename.
     * 
     * @param v container of elements.
     * 
     * @return the algebraic maximum of all the elements in `v'.
     */
    template <typename T>
    inline typename T::value_type Max(const T &v) {
        return *std::max_element(v.begin(),v.end());
    }

    /**
     * Compute the algebraic minimum of all the elements in a containes.
     * 
     * @tparam T container typename.
     * 
     * @param v container of elements.
     * 
     * @return the algebraic minimum of all the elements in `v'.
     */
    template <typename T>
    inline typename T::value_type Min(const T &v) {
        return *std::min_element(v.begin(),v.end());
    }

    /**
     * Compute the maximum absolute value of all the elements in a containes.
     * 
     * @tparam T container typename.
     * 
     * @param v container of elements.
     * 
     * @return the maximum absolute value of all the elements in `v'.
     */
    template <typename T>
    inline auto MaxAbs(const T &v) {
        using type = typename T::value_type;
        return std::abs(
            *std::max_element(v.begin(),v.end(),
                [] (const type &a, const type &b) -> bool {
                    return std::abs(a) < std::abs(b);
                }
            )
        );
    }

    /**
     * Compute the minimum absolute value of all the elements in a containes.
     * 
     * @tparam T container typename.
     * 
     * @param v container of elements.
     * 
     * @return the minimum absolute value of all the elements in `v'.
     */
    template <typename T>
    inline auto MinAbs(const T &v) {
        using type = typename T::value_type;
        return std::abs(
            *std::min_element(v.begin(),v.end(),
                [] (const type &a, const type &b) -> bool {
                    return std::abs(a) < std::abs(b);
                }
            )
        );
    }

    /**
     * Move from indices i,j,k to index idx.
     * 
     * @param i index along direction x.
     * @param j index along direction y.
     * @param k index along direction z.
     * @param n0 number of voxels along direction x.
     * @param n1 number of voxels along direction y.
     * 
     * @return the index idx.
     */
    inline size_t IJKToIdx(const size_t i, const size_t j, const size_t k,
    const size_t n0, const size_t n1) {
        return i + n0*(j + n1*k);
    }

    /**
     * Move from index idx to indices i,j,k.
     * 
     * @param[out] i index along direction x.
     * @param[out] j index along direction y.
     * @param[out] k index along direction z.
     * @param[in] idx index of the voxel.
     * @param[in] n0 number of voxels along direction x.
     * @param[in] n1 number of voxels along direction y.
     */
    inline void IdxToIJK(size_t &i, size_t &j, size_t &k, size_t idx,
    const size_t n0, const size_t n1) {
        i = idx%n0;
        idx /= n0;
        j = idx%n1;
        k = idx/n1;
    }

}  // eptlib

#endif  // EPTLIB_UTIL_H_
