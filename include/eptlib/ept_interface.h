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

#ifndef EPTLIB_EPT_INTERFACE_H_
#define EPTLIB_EPT_INTERFACE_H_

#include <numeric>
#include <array>

#include <boost/dynamic_bitset.hpp>

#include "eptlib/image.h"
#include "eptlib/median_filter.h"
#include "eptlib/shape.h"
#include "eptlib/util.h"

/// Namespace containing all symbols from the EPTlib library.
namespace eptlib {

/**
 * Abstract interface representing any EPT method.
 */
class EPTInterface {
    public:
        /**
         * Constructor.
         * 
         * @param freq operative frequency of the MRI scanner.
         * @param nn number of voxels in each direction.
         * @param dd voxel sizes in each direction.
         * @param tx_ch number of transmit channels.
         * @param rx_ch number of receive channels.
         */
        EPTInterface(const double freq, const std::array<int,NDIM> &nn,
            const std::array<double,NDIM> &dd, const int tx_ch=1,
            const int rx_ch=1);
        /**
         * Virtual destructor.
         */
        virtual ~EPTInterface() = 0;
        /**
         * Abstract method performing the tomography.
         * 
         * @return an error index specified by each method implementation.
         */
        virtual EPTlibError Run() = 0;
        /**
         * Set the transmit sensitivity of j-th Tx channel.
         * 
         * @param tx_sens pointer to transmit sensitivity image.
         * @param j transmit channel id.
         * 
         * @return a Success or OutOfRange error.
         * 
         * Note that just the _reference_ to tx_sens is stored.
         */
        EPTlibError SetTxSensitivity(const Image<double> *tx_sens,
            const int j = 0);
        /**
         * Set the transceive phase of j-th Tx channel w.r.t k-th Rx channel.
         * 
         * @param trx_phase pointer to transceive phase image.
         * @param j transmit channel id.
         * @param k receive channel id.
         * 
         * @return a Success or OutOfRange error.
         * 
         * Note that just the _reference_ to trx_phase is stored.
         */
        EPTlibError SetTRxPhase(const Image<double> *trx_phase,
            const int j = 0, const int k = 0);
        /**
         * Get the electric conductivity.
         * 
         * @param[out] sigma pointer to electric conductivity destination.
         * 
         * @return a Success or MissingData error.
         */
        EPTlibError GetElectricConductivity(Image<double> *sigma);
        /**
         * Get the relative permittivity.
         * 
         * @param[out] epsr pointer to relative permittivity destination.
         * 
         * @return a Success or MissingData error.
         */
        EPTlibError GetRelativePermittivity(Image<double> *epsr);
        /**
         * Set a post-processing median filter.
         * 
         * @param shape the shape of the mask over which apply the filter.
         * 
         * @return a Success or Unknown error.
         */
        EPTlibError SetPostPro(const Shape &shape);
        /**
         * Unset the post-processing filter.
         * 
         * @return a Success or Unknown error.
         */
        EPTlibError UnsetPostPro();
        /**
         * Apply the post-processing median filter to the output electric
         * conductivity and/or relative permittivity.
         * 
         * @param img pointer to the reference image. If it is not
         *     nullptr, then it is used as a reference.
         * 
         * @return a Success or MissingData error.
         */
        EPTlibError ApplyPostPro(const double *img = nullptr);
        /**
         * Set or unset the wrapped phase flag.
         * 
         * @return the updated wrapped phase flag.
         */
        bool TogglePhaseIsWrapped();
        /**
         * Get the wrapped phase flag.
         * 
         * @return the wrapped phase flag.
         */
        bool PhaseIsWrapped() const;
    protected:
        /// Operative harmonic frequency of the MRI scanner.
        double omega_;
        /// Number of transmit channels.
        int tx_ch_;
        /// Number of receive channels.
        int rx_ch_;
        /// Number of voxels in each direction.
        std::array<int,NDIM> nn_;
        /// Voxel sizes in each direction.
        std::array<double,NDIM> dd_;
        /// Total number of voxels.
        int n_vox_;
        /// Collection of pointers to transmit sensitivity distributions.
        std::vector<const Image<double>*> tx_sens_;
        /// Collection of pointers to transceive phase distributions (tx_ch faster).
        std::vector<const Image<double>*> trx_phase_;
        /// Transmit sensitivity data flag.
        boost::dynamic_bitset<> thereis_tx_sens_;
        /// Transceive phase data flag.
        boost::dynamic_bitset<> thereis_trx_phase_;
        /// Electric conductivity distribution.
        Image<double> sigma_;
        /// Relative permittivity distribution.
        Image<double> epsr_;
        /// Electric conductivity data flag.
        bool thereis_sigma_;
        /// Relative permittivity data flag.
        bool thereis_epsr_;
        /// Pointer to a post-processing median filter.
        MedianFilter *postpro_;
        /// Post-processing filter flag.
        bool thereis_postpro_;
        /// Wrapped phase flag.
        bool phase_is_wrapped_;
};

}  // namespace eptlib

#endif  // EPTLIB_EPT_INTERFACE_H_
