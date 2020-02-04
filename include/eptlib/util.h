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

#ifndef EPTLIB_UTIL_H_
#define EPTLIB_UTIL_H_

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
typedef enum EPTlibError {
  /// Success.
  Success = 0,
  /// Missing data error.
  MissingData,
  /// Out of range error.
  OutOfRange,
  /// Unknown error.
  Unknown,
} EPTlibError_t;
/**
 * Translates in a human-readable string the input EPTlibError symbol.
 *
 * @param error is an EPTlibError symbol.
 * @return the human-readable description of the input EPTlibError symbol.
 */
const std::string ToString(const EPTlibError_t error);

/**
 * Translates from multi-index to index assuming the first index the fastest.
 * 
 * @tparam T,U iterator typenames.
 * 
 * @param ii multi-index to the voxel.
 * @param nn number of voxels in each direction.
 * 
 * @return the single index to the voxel.
 * 
 * Arguments `ii' and `nn' must have the methods `begin', `end', `size' and
 * `operator[]' (any sequence containers from STL works fine).
 */
template <typename T, typename U>
int MultiIdxToIdx(const T &ii, const U &nn);
/**
 * Translates from index to multi-index assuming the first index the fastest.
 * 
 * @tparam T,U iterator typenames.
 * 
 * @param[out] ii multi-index to the voxel.
 * @param[in] idx single index to the voxel
 * @param[in] nn number of voxels in each direction.
 * 
 * Arguments `ii' and `nn' must have the methods `begin', `end', `size' and
 * `operator[]' (any sequence containers from STL works fine).
 */
template <typename T, typename U>
void IdxToMultiIdx(T &ii, int idx, const U &nn);

// ---------------------------------------------------------------------------
// -------------------------  Implementation detail  -------------------------
// ---------------------------------------------------------------------------

// Translates in a human-readable string the input EPTlibError symbol.
inline const std::string ToString(const EPTlibError_t error) {
    switch (error) {
        case EPTlibError::Success:
            return "Success";
        case EPTlibError::MissingData:
            return "Missing data";
        case EPTlibError::OutOfRange:
            return "Out of range";
        case EPTlibError::Unknown:
            return "Unknown error";
    }
    return "";
}

// Multi-index to index
template <typename T, typename U>
int MultiIdxToIdx(const T &ii, const U &nn) {
    assert(ii.size()==nn.size());
    auto n_dim = ii.size();
    int idx = 0;
    for (decltype(n_dim) d = n_dim; d>0; --d) {
        idx = ii[d-1] + idx*nn[d-1];
    }
    return idx;
}
// Index to multi-index
template <typename T, typename U>
void IdxToMultiIdx(T &ii, int idx, const U &nn) {
    assert(ii.size()==nn.size());
    auto n_dim = ii.size();
    for (decltype(n_dim) d = 0; d<n_dim; ++d) {
        ii[d] = idx%nn[d];
        idx = idx/nn[d];
    }
    return;
}

#endif  // EPTLIB_UTIL_H_
