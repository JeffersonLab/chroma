// -*- C++ -*-
// $Id: DjDk_w.h,v 2.0 2005-09-25 21:04:34 edwards Exp $

#ifndef __DjDk_h__
#define __DjDk_h__

namespace Chroma {

void DjDk(const multi1d<LatticeColorMatrix>& u,
	  const LatticePropagator& prop,
	  LatticePropagator& prop_jk,
	  int direction);

}  // end namespace Chroma

#endif
