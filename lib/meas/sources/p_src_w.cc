// $Id: p_src_w.cc,v 1.5 2005-01-14 18:42:37 edwards Exp $

#include "chromabase.h"
#include "meas/sources/p_src_w.h"

namespace Chroma {

//
// Take first derivative on source.
//

template<typename T>
void p_src(const multi1d<LatticeColorMatrix>& u, 
	   T& chi,
	   int direction)
{
  int d=direction;

  T psi = chi;

  chi = 0.5*(u[d]*shift(psi, FORWARD, d) - shift(adj(u[d])*psi, BACKWARD, d));
}


void p_src(const multi1d<LatticeColorMatrix>& u,
	   LatticeColorVector& chi,
	   int direction)
{
  p_src<LatticeColorVector>(u, chi, direction);
}

void p_src(const multi1d<LatticeColorMatrix>& u,
	   LatticePropagator& chi,
	   int direction)
{
  p_src<LatticePropagator>(u, chi, direction);
}

void p_src(const multi1d<LatticeColorMatrix>& u,
	   LatticeFermion& chi,
	   int direction)
{
  p_src<LatticeFermion>(u, chi, direction);
}

}  // end namespace Chroma
