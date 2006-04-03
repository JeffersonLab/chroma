// -*- C++ -*-
// $Id: stoch_var.h,v 3.0 2006-04-03 04:59:00 edwards Exp $
/*! \file
 *  \brief Stochastic variable construction
 *
 */

#ifndef __stoch_var_h__
#define __stoch_var_h__

namespace Chroma {

//! Stochastic variable construction
/*!
 * \ingroup hadron
 *
 * This routine averages timeslice sums of fermion disconnected loop 
 * operators over the number of stochastic sources.
 * It also calculates the standard deviation on the mean of the real
 * and imaginary parts of these operators.
 *
 * \param ferm_loop_sum    sum over stochastice samples of timeslice 
 *                         disconnected fermion loop operators
 * \param ferm_loop        The timeslice operators for EACH stochastic sample
 * \param sigma            standard deviation on the mean of the real part
 *                         of the operator
 * \param im_sigma         same as above for imaginary part
 * \param t_length         length of lattice in time dir
 * \param Nsamp            Number of stochastic samples.
 */

void 
stoch_var(multi1d<DComplex>& ferm_loop_sum, multi2d<DComplex>& ferm_loop, 
          multi1d<Real64>& sigma, multi1d<Real64>& imsigma, 
          int t_length, int Nsamp);

}  // end namespace Chroma

#endif
