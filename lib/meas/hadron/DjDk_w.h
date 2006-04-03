// -*- C++ -*-
// $Id: DjDk_w.h,v 3.0 2006-04-03 04:58:59 edwards Exp $

#ifndef __DjDk_h__
#define __DjDk_h__

namespace Chroma {

void DjDk(const multi1d<LatticeColorMatrix>& u,
	  const LatticePropagator& prop,
	  LatticePropagator& prop_jk,
	  int direction);

}  // end namespace Chroma

#endif
