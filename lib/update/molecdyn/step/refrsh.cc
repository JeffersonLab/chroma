// $Id: refrsh.cc,v 1.1 2003-12-30 19:52:28 edwards Exp $

#include "chromabase.h"

using namespace QDP;

//! Refrsh: refreshes momenta conjugate to u
/*!
 * \ingroup molecdyn
 *
 * \param p_mom    momenta conjugate to u (Write ) 
 */

void Refrsh(multi1d<LatticeColorMatrix>& p_mom)
{
  START_CODE("Refrsh");

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

  END_CODE("Refrsh");
}
