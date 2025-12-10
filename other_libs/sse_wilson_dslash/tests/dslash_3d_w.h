#ifndef DSLASH_3D_W_H
#define DSLASH_3D_W_H

#ifndef QDP_INCLUDE
#include "qdp.h"
#endif

using namespace QDP;

void dslash_3d(LatticeFermion& chi, 
	       const multi1d<LatticeColorMatrix>& u,
	       const LatticeFermion& psi,
	       int isign, int cb3d);

#endif
