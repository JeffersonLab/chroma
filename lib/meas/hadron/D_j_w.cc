// $Id: D_j_w.cc,v 1.3 2004-01-27 21:06:46 ikuro Exp $
/*
 *
 */
//

#include "chromabase.h"
#include "meas/hadron/D_j_w.h"

using namespace QDP;

void D_j(const multi1d<LatticeColorMatrix>& u,
	 const LatticePropagator& prop,
	 LatticePropagator& prop_j,
	 int direction)
{
  int d=direction;

  prop_j = 0.5*(u[d]*shift(prop, FORWARD, d) - shift(adj(u[d])*prop, BACKWARD, d));
}
