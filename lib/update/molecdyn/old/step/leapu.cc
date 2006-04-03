// $Id: leapu.cc,v 3.0 2006-04-03 04:59:10 edwards Exp $

#warning "NEED TO DEAL WITH SCHRFUN"

#include "chromabase.h"
#include "util/gauge/expmat.h"


//! LEAPU:  Make a step in U
/*!
 * \ingroup molecdyn
 *
 * p_mom -- conjugate momenta to u ( Read )
 * u     -- gauge field ( Modify )
 * eps   -- step size to use ( Read )
 */

void LeapU(multi1d<LatticeColorMatrix>& u,
	   const multi1d<LatticeColorMatrix>& p_mom,
	   const Real& eps)
{
  START_CODE();

  LatticeColorMatrix tmp_1;
  LatticeColorMatrix tmp_2;
  LatticeBoolean lbad;

  for(int mu=0; mu < Nd; ++mu)
  {
    /* tmp_1 = eps*p_mom */
    tmp_1 = p_mom[mu] * eps;

    /* tmp_1 = e^(tmp_1) */
    if ( Nc != 1 )
    {
      /* No need for exact exp. here. */
      expmat(tmp_1, TWELTH_ORDER);
    }
    else
    {
      /* No cost for exact exp. here. */
      expmat(tmp_1, EXACT);	
    }
    /* tmp_2 = u_old * tmp_1 = u_old * exp[ eps * p_mom ] */
    tmp_2 = tmp_1 * u[mu];
  
    /* Fill new u and reunitarize ---- */
    u[mu] = tmp_2;
    int numbad;
    reunit(u[mu], lbad, REUNITARIZE_ERROR, numbad); 
  }

#if 0
  /* Reset the boundary fields */
  /* This is not necessary, but useful for stability */
  if ( SchrFun > 0 )
  {
    copymask(u, lSFmask, SFBndFld, REPLACE);
  }
#endif

  END_CODE();
}

