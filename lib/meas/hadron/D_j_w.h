// -*- C++ -*-
// $Id: D_j_w.h,v 1.1 2003-08-19 17:15:24 bjoo Exp $

#ifndef __D_j_h__
#define __D_j_h__

void D_j(const multi1d<LatticeColorMatrix>& u,
	 const LatticePropagator& prop,
	 LatticePropagator& prop_j,
	 int direction);

#endif
