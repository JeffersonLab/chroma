// -*- C++ -*-
// $Id: D_j.h,v 1.1 2003-04-17 20:12:03 dgr Exp $

#ifndef __D_j_h__
#define __D_j_h__

void D_j(const multi1d<LatticeColorMatrix>& u,
	 const LatticePropagator& prop,
	 LatticePropagator& prop_j,
	 int direction);

#endif
