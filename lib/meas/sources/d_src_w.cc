// $Id: d_src_w.cc,v 1.4 2005-01-14 18:42:37 edwards Exp $

#include "chromabase.h"
#include "meas/sources/p_src_w.h"

namespace Chroma {

//
// Take second derivative to a source.
// direction: "12"  ->  DyDz chi
//            "22"  ->  DzDz chi
//

template<typename T>
void d_src(const multi1d<LatticeColorMatrix>& u, 
	   T& chi,
	   int direction)
{
  T psi = chi;
  
  if (direction == 12) {   // dydz
    chi = 0.5 * (u[2]*shift(psi, FORWARD, 2) - shift(adj(u[2])*psi, BACKWARD, 2));
    psi = chi;
    chi = 0.5 * (u[1]*shift(psi, FORWARD, 1) - shift(adj(u[1])*psi, BACKWARD, 1));
  }

  if (direction == 22) {   // dzdz
    chi = u[2]*shift(psi, FORWARD, 2) + shift(adj(u[2])*psi, BACKWARD, 2)
          - 2.0 * psi;
  }
}


void d_src(const multi1d<LatticeColorMatrix>& u,
	   LatticeColorVector& chi,
	   int direction)
{
  d_src<LatticeColorVector>(u, chi, direction);
}

void d_src(const multi1d<LatticeColorMatrix>& u,
	   LatticePropagator& chi,
	   int direction)
{
  d_src<LatticePropagator>(u, chi, direction);
}

void d_src(const multi1d<LatticeColorMatrix>& u,
	   LatticeFermion& chi,
	   int direction)
{
  d_src<LatticeFermion>(u, chi, direction);
}

}  // end namespace Chroma
