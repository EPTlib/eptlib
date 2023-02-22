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

#ifndef EPTLIB_EPT_INTERFACE_H_
#define EPTLIB_EPT_INTERFACE_H_

#include <array>
#include <memory>
#include <vector>

#include "eptlib/image.h"
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
         * @param n0 number of voxels along direction x.
         * @param n1 number of voxels along direction y.
         * @param n2 number of voxels along direction z.
         * @param d0 resolution in meter along direction x.
         * @param d1 resolution in meter along direction y.
         * @param d2 resolution in meter along direction z.
         * @param freq operative frequency of the MRI scanner.
         * @param n_tx_ch number of transmit channels (default: 1).
         * @param n_rx_ch number of receive channels (default: 1).
         * @param trx_phase_is_wrapped wrapped phase flag (default: false).
         */
        EPTInterface(const size_t n0, const size_t n1, const size_t n2,
            const double d0, const double d1, const double d2,
            const double freq, const size_t n_tx_ch = 1, const size_t n_rx_ch = 1,
            const bool trx_phase_is_wrapped = false);
        /**
         * Virtual destructor.
         */
        virtual ~EPTInterface();

        /**
         * Abstract method performing the tomography.
         * 
         * @return an error index specified by each method implementation.
         */
        virtual EPTlibError Run() = 0;

        /**
         * Set the transmit sensitivity of Tx channel tx_ch.
         * 
         * @param tx_sens transmit sensitivity image.
         * @param tx_ch transmit channel id.
         * 
         * @return Success, or OutOfRange or WrongDataFormat error.
         */
        EPTlibError SetTxSensitivity(const Image<double> &tx_sens,
            const size_t tx_ch = 0);

        /**
         * Set the transceive phase of Tx channel tx_ch w.r.t Rx channel rx_ch.
         * 
         * @param trx_phase transceive phase image.
         * @param tx_ch transmit channel id.
         * @param rx_ch receive channel id.
         * 
         * @return Success, or OutOfRange or WrongDataFormat error.
         */
        EPTlibError SetTRxPhase(const Image<double> &trx_phase,
            const size_t tx_ch = 0, const size_t rx_ch = 0);

        /**
         * @brief Get the pointer to the Tx sensitivity of the Tx channel tx_ch.
         * 
         * @param tx_ch index of the Tx channel.
         * 
         * @return the pointer to the Tx sensitivity of the Tx channel tx_ch.
         */
        inline const Image<double> * GetTxSens(const size_t tx_ch) const {
            return tx_sens_[tx_ch];
        }

        /**
         * @brief Get the pointer to the TRx phase of the Tx channel tx_ch and the Rx
         * channel rx_ch.
         * 
         * @param tx_ch index of the Tx channel.
         * @param rx_ch index of the Rx channel.
         * 
         * @return the pointer to the TRx phase of the Tx channel tx_ch and the Rx channel rx_ch.
         */
        inline const Image<double> * GetTRxPhase(const size_t tx_ch, const size_t rx_ch) const {
            return trx_phase_[tx_ch+n_tx_ch_*rx_ch];
        }

        /**
         * Get the electric conductivity.
         * 
         * @return pointer to the electric conductivity.
         */
        inline std::unique_ptr<Image<double> >& GetElectricConductivity() {
            return sigma_;
        };

        /**
         * Get the relative permittivity.
         * 
         * @param[out] epsr pointer to relative permittivity destination.
         * 
         * @return a Success or MissingData error.
         */
        inline std::unique_ptr<Image<double> >& GetRelativePermittivity() {
            return epsr_;
        };

        /**
         * Set or unset the wrapped phase flag.
         * 
         * @return the updated wrapped phase flag.
         */
        inline bool TogglePhaseIsWrapped() {
            trx_phase_is_wrapped_ = !trx_phase_is_wrapped_;
            return trx_phase_is_wrapped_;
        };

        /**
         * Get the wrapped phase flag.
         * 
         * @return the wrapped phase flag.
         */
        inline bool PhaseIsWrapped() const {
            return trx_phase_is_wrapped_;
        };

        /**
         * @brief Check if the Tx sensitivity of the Tx channel tx_ch is set.
         * 
         * @param tx_ch index of the Tx channel.
         * 
         * @return true if the Tx sensitivity is set.
         * @return false if the Tx sensitivity is not set.
         */
        inline bool ThereIsTxSens(const size_t tx_ch) const {
            return tx_sens_[tx_ch]!=nullptr;
        }

        /**
         * @brief Check if all the Tx sensitivities are set.
         * 
         * @return true if all the Tx sensitivities are set.
         * @return false if any Tx sensitivity is not set.
         */
        inline bool ThereAreAllTxSens() const {
            for (size_t tx_ch = 0; tx_ch<n_tx_ch_; ++tx_ch) {
                if (!ThereIsTxSens(tx_ch)) {
                    return false;
                }
            }
            return true;
        }

        /**
         * @brief Check if the TRx phase of the Tx channel tx_ch and the Rx
         * channel rx_ch is set.
         * 
         * @param tx_ch index of the Tx channel.
         * @param rx_ch index of the Rx channel.
         * 
         * @return true if the TRx phase is set.
         * @return false if the TRx phase is not set.
         */
        inline bool ThereIsTRxPhase(const size_t tx_ch, const size_t rx_ch) const {
            return trx_phase_[tx_ch+n_tx_ch_*rx_ch]!=nullptr;
        }

        /**
         * @brief Check if all the TRx phases are set.
         * 
         * @return true if all the TRx phases are set.
         * @return false if any TRx phase is not set.
         */
        inline bool ThereAreAllTRxPhase() const {
            for (size_t rx_ch = 0; rx_ch<n_rx_ch_; ++rx_ch) {
                for (size_t tx_ch = 0; tx_ch<n_tx_ch_; ++tx_ch) {
                    if (!ThereIsTRxPhase(tx_ch,rx_ch)) {
                        return false;
                    }
                }
            }
            return true;
        }

        /**
         * @brief Check if the electric conductivity is set.
         * 
         * @return true if the electric conductivity is set.
         * @return false if the electric conductivity is not set.
         */
        inline bool ThereIsSigma() const {
            return sigma_!=nullptr;
        }

        /**
         * @brief Check if the relative permittivity is set.
         * 
         * @return true if the relative permittivity is set.
         * @return false if the relative permittivity is not set.
         */
        inline bool ThereIsEpsr() const {
            return epsr_!=nullptr;
        }
    protected:
        /// Number of voxels in each direction.
        std::array<size_t,N_DIM> nn_;
        /// Voxel sizes in each direction.
        std::array<double,N_DIM> dd_;
        /// Operative harmonic frequency of the MRI scanner.
        double omega_;
        /// Number of transmit channels.
        int n_tx_ch_;
        /// Number of receive channels.
        int n_rx_ch_;
        /// Collection of pointers to transmit sensitivity distributions.
        std::vector<const Image<double> *> tx_sens_;
        /// Collection of pointers to transceive phase distributions (tx_ch faster).
        std::vector<const Image<double> *> trx_phase_;
        /// Electric conductivity distribution.
        std::unique_ptr<Image<double> > sigma_;
        /// Relative permittivity distribution.
        std::unique_ptr<Image<double> > epsr_;
        /// If true, the TRx phase is wrapped.
        bool trx_phase_is_wrapped_;

        /**
         * @brief Check if two arrays of sizes are equal component by component.
         * 
         * @param nn1 array 1
         * @param nn2 array 2
         * @return true if the arrays are equal component by component.
         * @return false if at least one component of the arrays is different.
         */
        bool CheckSizes(const std::array<size_t,N_DIM> &nn1,
            const std::array<size_t,N_DIM> &nn2) const {
            return nn1[0]==nn2[0] && nn1[1]==nn2[1] && nn1[2]==nn2[2];
        }
};

}  // namespace eptlib

#endif  // EPTLIB_EPT_INTERFACE_H_
