// -*- C++ -*-
// $Id: dslashm_w.h,v 1.1 2003-04-09 17:19:40 edwards Exp $

#ifndef DSLASHM_INCLUDE
#define DSLASHM_INCLUDE

void dslash(LatticeFermion& chi, const multi1d<LatticeColorMatrix>& u, const LatticeFermion& psi,
	    int isign, int cb);

#endif
