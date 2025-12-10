#ifndef DSLASH_M_W_H
#define DSLASH_M_W_H

#ifndef QDP_INCLUDE
#include "qdp.h"
#endif

using namespace QDP;

void dslash(LatticeFermion& chi, 
	    const multi1d<LatticeColorMatrix>& u,
	    const LatticeFermion& psi,
	    int isign, int cb);

#endif
