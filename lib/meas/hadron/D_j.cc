// $Id: D_j.cc,v 1.1 2003-04-17 20:12:03 dgr Exp $
/*
 *
 */
//

#include "chromabase.h"
#include "meas/hadron/D_j.h"

using namespace QDP;

void D_j(const multi1d<LatticeColorMatrix>& u,
	 const LatticePropagator& prop,
	 LatticePropagator& prop_j,
	 int direction)
{
  int d=direction;

  prop_j = u[d]*shift(prop, FORWARD, d) - shift(adj(u[d])*prop, BACKWARD, d);
}
