// -*- C++ -*-
// $Id: stoch_var.cc,v 3.0 2006-04-03 04:59:00 edwards Exp $
/*! \file
 *  \brief Stochastic variable construction
 *
 */

#include "chromabase.h"
#include "meas/hadron/stoch_var.h"

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
          multi1d<Real64>& sigma, multi1d<Real64>& im_sigma, 
          int t_length, int Nsamp)
{
  multi1d<Real64> mean_sq(t_length), im_mean_sq(t_length), s(t_length),
                  im_s(t_length), std_dev(t_length), im_std_dev(t_length),
                  re_loop(t_length), im_loop(t_length);

  s = im_s = zero;

  for (int t = 0; t < t_length; ++t){

    // Average over number of stochastic samples
    re_loop[t] = real(ferm_loop_sum[t])/Nsamp;
    im_loop[t] = imag(ferm_loop_sum[t])/Nsamp;
    ferm_loop_sum[t] = cmplx(re_loop[t], im_loop[t]);

    // Calculate the standard deviation on the mean of the
    // fermion loop operators.
    // Do this for real and imaginary parts
 
    // Square of the mean
    mean_sq[t] = pow(re_loop[t], 2);
    im_mean_sq[t] = pow(im_loop[t], 2);

    // s = sum( (x_i)^2)
    for (int i = 0; i < Nsamp; ++i){
      s[t] += pow(real(ferm_loop[i][t]) , 2);
      im_s[t] += pow(imag(ferm_loop[i][t]) , 2);
    }

    // std_dev = sqrt(< (x_i)^2 > - mean^2)
    std_dev[t] = sqrt((s[t]/Nsamp) - mean_sq[t]);
    im_std_dev[t] = sqrt((im_s[t]/Nsamp) - im_mean_sq[t]);

    // std_dev on mean = 1/(sqrt(Nsamp-1) * std_dev
    sigma[t] = (1/sqrt(Nsamp-1.0))*std_dev[t];
    im_sigma[t] = (1/sqrt(Nsamp-1.0))*im_std_dev[t];
  }

}    

}  // end namespace Chroma
