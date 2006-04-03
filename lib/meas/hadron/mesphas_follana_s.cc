// $Id: mesphas_follana_s.cc,v 3.0 2006-04-03 04:59:00 edwards Exp $


/* This routine is specific to staggered fermions! */

/* This routine will calculate the Kogus-Sussking Phases (alpha)
 * and the other phases (beta) needed to compute spectroscopy. 
 * -- GLASGOW/UKQCD Phase Conventions:
 * 
 * alpha(x,0) = 1
 * alpha(x,mu) = (-1)^{ x_0 + x_1 + ... x_{mu-1} } 
 * 
 * beta(x,mu) = (-1)^{ x_{mu+1} + ... + x_{Nd-1} }
 * beta(x,Nd-1) = 1; 
 *
 */

/* alpha -- KS phases (Write)
/* beta  -- the beta phase factors ( Write ) */


#include "chromabase.h"
#include "mesphas_follana_s.h"

namespace Chroma {

  void mesPhasFollana(multi1d<LatticeInteger>& alpha,
		      multi1d<LatticeInteger>& beta)
  {
    multi1d<LatticeInteger> x(Nd);
    int mu;

    alpha.resize(Nd);
    beta.resize(Nd);

   
    /* Start off by getting the coordinates of x(0), x(1), ..., x(Nd-3) */
    for( mu = 0; mu < Nd; mu++) { 
      x[ mu ] = Layout::latticeCoordinate(mu);
    }
  
  
    /* NOTE: in the comments below I mean x0 is the true *non*-checkerboard */
    /*   x coordinate. */
    switch(Nd)
    {
    case 4:
            
      alpha[0] = LatticeInteger(1);
      alpha[1] = where((x[0] % 2) == 0, LatticeInteger(1), LatticeInteger(-1));
      alpha[2] = where( ((x[0]+x[1])%2) == 0, LatticeInteger(1), LatticeInteger(-1));
      alpha[3] = where( ((x[0]+x[1]+x[2])%2) == 0, LatticeInteger(1), LatticeInteger(-1));

      beta[0] = where( ((x[1]+x[2]+x[3])%2) == 0, LatticeInteger(1), LatticeInteger(-1));
      beta[1] = where( ((x[2] + x[3])%2) == 0, LatticeInteger(1), LatticeInteger(-1));
      beta[2] = where( (x[3] % 2) == 0, LatticeInteger(1), LatticeInteger(-1) );

      beta[3] = LatticeInteger(1);

      break;

    default:
      QDP_error_exit("Can only handle d=4 dimensions for now", Nd);
    }
  
            
  }

}  // end namespace Chroma
