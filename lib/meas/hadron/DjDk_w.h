// -*- C++ -*-
// $Id: DjDk_w.h,v 1.2 2003-10-10 03:46:46 edwards Exp $

#ifndef __DjDk_h__
#define __DjDk_h__

void DjDk(const multi1d<LatticeColorMatrix>& u,
	  const LatticePropagator& prop,
	  LatticePropagator& prop_jk,
	  int direction);

#endif
