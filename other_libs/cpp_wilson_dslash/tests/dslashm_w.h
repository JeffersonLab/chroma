#ifndef DSLASH_M_W_H
#define DSLASH_M_W_H

#ifndef QDP_INCLUDE
#include "qdp.h"
#endif

using namespace QDP;


void dslash(LatticeFermionF& chi, 
	    const multi1d<LatticeColorMatrixF>& u, 
	    const LatticeFermionF& psi,
	    int isign, int cb);


void dslash(LatticeFermionD& chi, 
	    const multi1d<LatticeColorMatrixD>& u, 
	    const LatticeFermionD& psi,
	    int isign, int cb);

#endif
