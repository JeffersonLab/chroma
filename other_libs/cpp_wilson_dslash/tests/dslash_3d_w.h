#ifndef DSLASH_3D_W_H
#define DSLASH_3D_W_H

#ifndef QDP_INCLUDE
#include "qdp.h"
#endif

using namespace QDP;


void dslash_3d(LatticeFermionF& chi, 
	       const multi1d<LatticeColorMatrixF>& u,
	       const LatticeFermionF& psi,
	       int isign, int cb3d);

void dslash_3d(LatticeFermionD& chi, 
	       const multi1d<LatticeColorMatrixD>& u,
	       const LatticeFermionD& psi,
	       int isign, int cb3d);

#endif
