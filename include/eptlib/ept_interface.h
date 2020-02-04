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

#ifndef EPTLIB_EPT_INTERFACE_H_
#define EPTLIB_EPT_INTERFACE_H_

#include <algorithm>
#include <functional>
#include <numeric>
#include <string>
#include <vector>

#include <boost/dynamic_bitset.hpp>

#include "eptlib/config.h"
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
         * @tparam T,U iterator typenames.
         * 
         * @param freq operative frequency of the MRI scanner.
         * @param nn number of voxels in each direction.
         * @param dd voxel sizes in each direction.
         * @param tx_ch number of transmit channels.
         * @param rx_ch number of receive channels.
         * 
         * Arguments `nn' and `dd' must have the methods `begin', `end',
         * `size' and `operator[]' (any sequence container from STL works
         * fine).
         */
        template <typename T, typename U>
        EPTInterface(const real_t freq, const T &nn, const U &dd, const int tx_ch=1, const int rx_ch=1);
        /**
         * Virtual destructor.
         */
        virtual ~EPTInterface();
        /**
         * Abstract method performing the tomography.
         * 
         * @return an error index specified by each method implementation.
         */
        virtual EPTlibError_t Run() = 0;
        /**
         * Set the transmit sensitivity of j-th Tx channel.
         * 
         * @param tx_sens pointer to transmit sensitivity distribution.
         * @param j transmit channel id.
         * 
         * @return a Success or OutOfRange error.
         * 
         * Note that just the _reference_ to tx_sens is stored.
         */
        EPTlibError_t SetTxSensitivity(const real_t *tx_sens, const int j = 0);
        /**
         * Set the transceive phase of j-th Tx channel w.r.t k-th Rx channel.
         * 
         * @param trx_phase pointer to transceive phase distribution.
         * @param j transmit channel id.
         * @param k receive channel id.
         * 
         * @return a Success or OutOfRange error.
         * 
         * Note that just the _reference_ to trx_phase is stored.
         */
        EPTlibError_t SetTRxPhase(const real_t *trx_phase, const int j = 0, const int k = 0);
        /**
         * Get the electric conductivity.
         * 
         * @param[out] sigma pointer to electric conductivity destination.
         * 
         * @return a Success or MissingData error.
         */
        EPTlibError_t GetElectricConductivity(real_t *sigma);
        /**
         * Get the relative permittivity.
         * 
         * @param[out] epsr pointer to relative permittivity destination.
         * 
         * @return a Success or MissingData error.
         */
        EPTlibError_t GetRelativePermittivity(real_t *epsr);
    protected:
        /// Operative harmonic frequency of the MRI scanner.
        real_t omega_;
        /// Number of spatial dimensions.
        int n_dim_;
        /// Number of transmit channels.
        int tx_ch_;
        /// Number of receive channels.
        int rx_ch_;
        /// Number of voxels in each direction.
        std::vector<int> nn_;
        /// Voxel sizes in each direction.
        std::vector<real_t> dd_;
        /// Total number of voxels.
        int n_vox_;
        /// Collection of pointers to transmit sensitivity distributions.
        std::vector<const real_t*> tx_sens_;
        /// Collection of pointers to transceive phase distributions (tx_ch faster).
        std::vector<const real_t*> trx_phase_;
        /// Transmit sensitivity data flag.
        boost::dynamic_bitset<> thereis_tx_sens_;
        /// Transceive phase data flag.
        boost::dynamic_bitset<> thereis_trx_phase_;
        /// Electric conductivity distribution.
        std::vector<real_t> sigma_;
        /// Relative permittivity distribution.
        std::vector<real_t> epsr_;
        /// Electric conductivity data flag.
        bool thereis_sigma_;
        /// Relative permittivity data flag.
        bool thereis_epsr_;
};

// ---------------------------------------------------------------------------
// -------------------------  Implementation detail  -------------------------
// ---------------------------------------------------------------------------

// EPTInterface constructor
template <typename T, typename U>
EPTInterface::
EPTInterface(const real_t freq, const T &nn, const U &dd, const int tx_ch, const int rx_ch) :
    omega_(2.0*PI*freq), n_dim_(static_cast<int>(nn.size())), tx_ch_(tx_ch), rx_ch_(rx_ch) {
    // check for coherence between nn and dd sizes
    assert(n_dim_==dd.size());
    // copy the content of nn and dd
    nn_.resize(n_dim_);
    dd_.resize(n_dim_);
    std::copy(nn.begin(),nn.end(),nn_.begin());
    std::copy(dd.begin(),dd.end(),dd_.begin());
    n_vox_ = std::accumulate(nn_.begin(),nn_.end(),1,std::multiplies<int>());
    // initialise tx_sens_ and trx_phase_
    tx_sens_.resize(tx_ch_,nullptr);
    trx_phase_.resize(rx_ch_*tx_ch_,nullptr);
    // initialise the data flags
    thereis_tx_sens_.resize(tx_ch_,false);
    thereis_trx_phase_.resize(rx_ch_*tx_ch_,false);
    thereis_sigma_ = false;
    thereis_epsr_ = false;
    return;
}

}  // namespace eptlib

#endif  // EPTLIB_EPT_INTERFACE_H_
