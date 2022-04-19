/*****************************************************************************
*
*     Program: EPTlib
*     Author: Alessandro Arduino <a.arduino@inrim.it>
*
*  MIT License
*
*  Copyright (c) 2020-2021  Alessandro Arduino
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

#ifndef EPTLIB_EPT_HELMHOLTZ_CHI2_H_
#define EPTLIB_EPT_HELMHOLTZ_CHI2_H_

#include "eptlib/ept_interface.h"

#include <array>
#include <vector>

#include "eptlib/ept_helmholtz.h"
#include "eptlib/finite_difference.h"
#include "eptlib/shape.h"
#include "eptlib/util.h"

namespace eptlib {

/**
 * @brief Implementation of the Helmholtz-based EPT method with pixel-wise
 * optimised kernel shape.
 * 
 * The method is phase-based and provides only an estimation of the electric
 * conductivity. Because of the involved non-linearities, it can operate only
 * on unwrapped phase maps.
 */
class EPTHelmholtzChi2 : public EPTInterface {
public:
    /**
     * @brief Constructor
     * 
     * @param freq operative frequency of the MRI scanner.
     * @param nn number of voxels in each direction.
     * @param dd voxels sizes in each direction.
     * @param shapes list of masks over which apply the finite difference scheme.
     * @param degree degree of the interpolating polynomial for the finite
     *     difference scheme (default: 2).
     * 
     * The number of Tx and Rx channels is fixed equal to one.
     */
    EPTHelmholtzChi2(const double freq,const std::array<int,NDIM> &nn,
        const std::array<double,NDIM> &dd,const std::vector<Shape> &shapes,
        const int degree = 2);
    /**
     * @brief Virtual destructor.
     */
    virtual ~EPTHelmholtzChi2();
    /**
     * @brief Perform the Helmholtz-based EPT with pixel-wise optimised kernel shape.
     * 
     * @return an error index about the state of the tomography.
     */
    virtual EPTlibError Run() override;
    /**
     * @brief Get the result variance.
     * 
     * @param[out] var pointer to the variance destination.
     * 
     * @return a Success or MissingData error.
     */
    EPTlibError GetVar(Image<double> *var);
    /**
     * @brief Get the pixel-wise selected kernel shape.
     * 
     * @param[out] shape_index pointer to the selected shape destination.
     * 
     * @return a Success or MissingData error.
     */
    EPTlibError GetShapeIndex(Image<int> *shape_index);
    /**
     * Set/unset unphysical values as admittable.
     * 
     * @return if unphysical values are admittable.
     */
    bool ToggleUnphysicalValues();
private:
    /// Frequency of the MRI scanner.
    double freq_;
    /// List of masks over which apply the finite difference scheme.
    std::vector<Shape> shapes_;
    /// Degree of the interpolating polynomial for the finite difference scheme.
    int degree_;
    /// Quality map.
    Image<double> var_;
    /// Pixel-wise selected kernel shape.
    Image<int> shape_index_;
    /// Unphysical values flag.
    bool unphysical_values_;
};

}  // namespace eptlib

#endif  // EPTLIB_EPT_HELMHOLTZ_CHI2_H_
