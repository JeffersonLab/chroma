// $Id: D_j_w.cc,v 1.4 2005-01-14 18:42:35 edwards Exp $
/*
 *
 */
//

#include "chromabase.h"
#include "meas/hadron/D_j_w.h"

namespace Chroma {

void D_j(const multi1d<LatticeColorMatrix>& u,
	 const LatticePropagator& prop,
	 LatticePropagator& prop_j,
	 int direction)
{
  int d=direction;

  prop_j = 0.5*(u[d]*shift(prop, FORWARD, d) - shift(adj(u[d])*prop, BACKWARD, d));
}

}  // end namespace Chroma
