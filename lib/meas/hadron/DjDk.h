// -*- C++ -*-
// $Id: DjDk.h,v 1.1 2003-06-19 17:45:26 ikuro Exp $

#ifndef __DjDk_h__
#define __DjDk_h__

void DjDk(const multi1d<LatticeColorMatrix>& u,
	  const LatticePropagator& prop,
	  LatticePropagator& prop_jk,
	  int direction);

#endif
