// $Id: DjDk_w.cc,v 2.0 2005-09-25 21:04:34 edwards Exp $
/*
 *
 */
//

#include "chromabase.h"
#include "meas/hadron/DjDk_w.h"

namespace Chroma {

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

  else if (direction == 00) {   // dxdx
    prop_jk = u[0]*shift(prop, FORWARD, 0) + shift(adj(u[0])*prop, BACKWARD, 0)
             - 2.0 * prop;
  }

  else if (direction == 11) {   // dydy
    prop_jk = u[1]*shift(prop, FORWARD, 1) + shift(adj(u[1])*prop, BACKWARD, 1)
             - 2.0 * prop;
  }

  else if (direction == 22) {   // dzdz
    prop_jk = u[2]*shift(prop, FORWARD, 2) + shift(adj(u[2])*prop, BACKWARD, 2)
             - 2.0 * prop;
  }

  else {
    cout << "Invalid direction\n"
	 << "Exit here\n";

    exit(1);
      }
}

}  // end namespace Chroma
