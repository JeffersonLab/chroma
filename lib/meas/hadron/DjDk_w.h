// -*- C++ -*-
// $Id: DjDk_w.h,v 1.3 2005-01-14 18:42:35 edwards Exp $

#ifndef __DjDk_h__
#define __DjDk_h__

namespace Chroma {

void DjDk(const multi1d<LatticeColorMatrix>& u,
	  const LatticePropagator& prop,
	  LatticePropagator& prop_jk,
	  int direction);

}  // end namespace Chroma

#endif
