// $Id: D_j_w.cc,v 1.2 2003-08-20 11:33:51 bjoo Exp $
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

  prop_j = u[d]*shift(prop, FORWARD, d) - shift(adj(u[d])*prop, BACKWARD, d);
}
