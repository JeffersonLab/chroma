/*! File: stoch_var.cc
 *
 * This routine averages timeslice sums of fermion disconnected loop 
 * operators over the number of stochastic sources.
 * It also calculates the standard deviation on the mean of the real
 * and imaginary parts of these operators.
 * 
 * Params:
 *
 *      ferm_loop_sum -- sum over stochastice samples of timeslice 
 *                       disconnected fermion loop operators
 *      ferm_loop     -- The timeslice operators for EACH stochastic sample
 *      sigma         -- standard deviation on the mean of the real part
 *                       of the operator
 *      im_sigma      -- same as above for imaginary part
 *      t_length      -- length of lattice in time dir
 *      Nsamp         -- Number of stochastic samples.
 */

#include "chroma.h"
#ifndef __stoch_var_h__
#define __stoch_var_h__

void 
stoch_var(multi1d<DComplex>& ferm_loop_sum, multi2d<DComplex>& ferm_loop, 
          multi1d<Real64>& sigma, multi1d<Real64>& imsigma, 
          int t_length, int Nsamp);

#endif
