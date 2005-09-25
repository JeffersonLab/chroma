// -*- C++ -*-
// $Id: D_j_w.h,v 2.0 2005-09-25 21:04:34 edwards Exp $

#ifndef __D_j_h__
#define __D_j_h__

namespace Chroma {

void D_j(const multi1d<LatticeColorMatrix>& u,
	 const LatticePropagator& prop,
	 LatticePropagator& prop_j,
	 int direction);

}  // end namespace Chroma

#endif
