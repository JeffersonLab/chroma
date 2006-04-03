// $Id: refrsh.cc,v 3.0 2006-04-03 04:59:11 edwards Exp $

#include "chromabase.h"


//! Refrsh: refreshes momenta conjugate to u
/*!
 * \ingroup molecdyn
 *
 * \param p_mom    momenta conjugate to u (Write ) 
 */

void Refrsh(multi1d<LatticeColorMatrix>& p_mom)
{
  START_CODE();

  /* Generate new momenta */
  for(int mu = 0; mu < Nd; ++mu)
  {
    gaussian(p_mom[mu]);
    taproj(p_mom[mu]);
  }

#if 0
  /* If using Schroedinger functional, zero out the boundaries */
  if (SchrFun > 0)
  {
    FILLMASK(p_mom, lSFmask, ZERO);
  }
#endif

  END_CODE();
}
