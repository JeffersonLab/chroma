// -*- C++ -*-
// $Id: D_j_w.h,v 3.0 2006-04-03 04:58:59 edwards Exp $

#ifndef __D_j_h__
#define __D_j_h__

namespace Chroma {

void D_j(const multi1d<LatticeColorMatrix>& u,
	 const LatticePropagator& prop,
	 LatticePropagator& prop_j,
	 int direction);

}  // end namespace Chroma

#endif
