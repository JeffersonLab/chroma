// $Id: d_src_w.cc,v 1.2 2003-08-19 17:34:14 bjoo Exp $

#include "chromabase.h"
#include "meas/sources/p_src_w.h"

using namespace QDP;

void d_src(const multi1d<LatticeColorMatrix>& u, 
              LatticeFermion& chi,
	      int direction)
{
  LatticeFermion psi = chi;
  
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

