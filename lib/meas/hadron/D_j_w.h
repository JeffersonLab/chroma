// -*- C++ -*-
// $Id: D_j_w.h,v 1.2 2003-10-10 03:46:46 edwards Exp $

#ifndef __D_j_h__
#define __D_j_h__

void D_j(const multi1d<LatticeColorMatrix>& u,
	 const LatticePropagator& prop,
	 LatticePropagator& prop_j,
	 int direction);

#endif
