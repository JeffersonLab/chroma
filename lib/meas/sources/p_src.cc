// $Id: p_src.cc,v 1.1 2003-04-17 14:50:22 dgr Exp $

#include "chromabase.h"
#include "meas/sources/p_src.h"

using namespace QDP;

void p_src(const multi1d<LatticeColorMatrix>& u, 
              LatticeFermion& chi,
	      int direction)
{
  int d=direction;

  LatticeFermion psi = chi;

  chi = u[d]*shift(psi, FORWARD, d) - shift(adj(u[d])*psi, BACKWARD, d);
}

