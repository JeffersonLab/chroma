// -*- C++ -*-
// $Id: D_j_w.h,v 1.3 2005-01-14 18:42:35 edwards Exp $

#ifndef __D_j_h__
#define __D_j_h__

namespace Chroma {

void D_j(const multi1d<LatticeColorMatrix>& u,
	 const LatticePropagator& prop,
	 LatticePropagator& prop_j,
	 int direction);

}  // end namespace Chroma

#endif
