// $Id: DjDk.cc,v 1.1 2003-06-19 17:45:26 ikuro Exp $
/*
 *
 */
//

#include "chromabase.h"
#include "meas/hadron/D_j.h"

using namespace QDP;

void DjDk(const multi1d<LatticeColorMatrix>& u,
	  const LatticePropagator& prop,
	  LatticePropagator& prop_jk,
	  int direction)
  // DjDkq:  k=direction1, j=direction2
{
  if (direction == 12) {   // dydz
    prop_jk = 0.5 * (u[2]*shift(prop, FORWARD, 2) - shift(adj(u[2])*prop, BACKWARD, 2));

    LatticePropagator p_tmp = prop_jk;

    prop_jk = 0.5 * (u[1]*shift(p_tmp, FORWARD, 1) - shift(adj(u[1])*p_tmp, BACKWARD, 1));
  }

  else if (direction == 11) {   // dydy
    prop_jk = u[1]*shift(prop, FORWARD, 1) + shift(adj(u[1])*prop, BACKWARD, 1)
             - 2.0 * prop;
  }

  else if (direction == 22) {   // dzdz
    prop_jk = u[2]*shift(prop, FORWARD, 2) + shift(adj(u[2])*prop, BACKWARD, 2)
             - 2.0 * prop;
  }
}
