// -*- C++ -*-
// $Id: dslashm_w.h,v 1.2 2003-10-10 03:46:46 edwards Exp $

#ifndef DSLASHM_INCLUDE
#define DSLASHM_INCLUDE

void dslash(LatticeFermion& chi, const multi1d<LatticeColorMatrix>& u, const LatticeFermion& psi,
	    int isign, int cb);

#endif
