// $Id: mesmom.cc,v 1.1 2003-12-30 19:50:25 edwards Exp $

#warning "SHOULD ADD INCLUDE OF PROTOTYPE"

#include "chromabase.h"

using namespace QDP;

//! MESMOM: calculate (p_mom,p_mom)
/*!
 * \ingroup molecdyn
 *
 */

/* p_mom -- momenta conjugate to u (Read ) */
/* p_mom_sq = | p_mom |^2  (Write) */

void MesMom(const multi1d<LatticeColorMatrix>& p_mom,
	    Double& p_mom_sq)
{
  START_CODE("MesMom");

  p_mom_sq = 0;
  for(int mu=0; mu < Nd; ++mu)
    p_mom_sq += norm2(p_mom[mu]);
  
  p_mom_sq /= Double(Layout::vol()*Nd);

  END_CODE("MesMom");
}
