// -*- C++ -*-
// $Id: sun_proj.h,v 1.1 2003-03-30 15:35:22 edwards Exp $

#ifndef __sun_proj__
#define __sun_proj__

//! Project a complex Nc x Nc matrix W onto SU(Nc) by maximizing Tr(VW)
/*!
 * Arguments:
 *
 *  \param w            complex Nc x Nc matrix (Read)
 *  \param v            the projected SU(Nc) Matrix (Write)
 *  \param BlkAccu      accuracy in SU(Nc) projection (Read)
 *  \param BlkMax       max number of iterations in SU(Nc) projection (Read)
 */

void sun_proj(const LatticeColorMatrix& w, LatticeColorMatrix& v,
	      const Real& BlkAccu, int BlkMax);

#endif
