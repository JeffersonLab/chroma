// -*- C++ -*-
// $Id: DjDk_w.h,v 1.1 2003-08-19 17:15:24 bjoo Exp $

#ifndef __DjDk_h__
#define __DjDk_h__

void DjDk(const multi1d<LatticeColorMatrix>& u,
	  const LatticePropagator& prop,
	  LatticePropagator& prop_jk,
	  int direction);

#endif
