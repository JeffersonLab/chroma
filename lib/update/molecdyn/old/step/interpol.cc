// $Id: interpol.cc,v 3.0 2006-04-03 04:59:10 edwards Exp $

#error "NOT FULLY CONVERTED - NEED TO MOVE IntrplOrd into params of Integ. functor"

#include "chromabase.h"


//! Linearly interpolate psi and old_psi
/*!
 * \ingroup molecdyn
 *
 */
/* psi     -- present psi ( Modify ) */
/* old_psi -- previous psi ( Modify ) */
/* fract   -- fractional scale factor for interpolation ( Read ) */
/* Ncb     -- number of checkerboards  ( Read ) */
/* Npf     -- number of pseudofermions  ( Read ) */

void Interpol(multi1d<LatticeFermion>& psi,
	      multi1d<LatticeFermion>& old_psi,
	      const Real& fract,
	      int Npf)
{
  START_CODE();
  
  if ( FermiP == NO ) 
  {
    END_CODE();
    return;
  }
  
  switch(IntrplOrd)       // Uses global  int IntrplOrd
  {
  case -1:
    /* initial guess for solution = 0 */
    psi = 0;
    break;

  case 0:
    /* initial guess for solution = previous solution */
    break;

  case 1:
    /* initial guess for solution = linear interpolation */
    /* psi = (1-fract)*psi + fract*old_psi */
    for(int i = 0; i < Npf; ++i)
    {
      LatticeFermion tmp = (1-fract)*psi[i]  +  old_psi[i]*fract;
      old_psi[i] = psi[i];
      psi[i] = tmp;
    }
    break;

  default:
    QDP_error_exit("interpolation order not implemented", IntrplOrd);
  }
  
  END_CODE();
}
