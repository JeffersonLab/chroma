// $Id: p_src_w.cc,v 1.3 2004-01-27 21:07:32 ikuro Exp $

#include "chromabase.h"
#include "meas/sources/p_src_w.h"

using namespace QDP;

void p_src(const multi1d<LatticeColorMatrix>& u, 
              LatticeFermion& chi,
	      int direction)
{
  int d=direction;

  LatticeFermion psi = chi;

  chi = 0.5*(u[d]*shift(psi, FORWARD, d) - shift(adj(u[d])*psi, BACKWARD, d));
}

