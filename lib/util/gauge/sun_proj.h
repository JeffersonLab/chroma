// -*- C++ -*-
// $Id: sun_proj.h,v 1.2 2003-12-29 19:52:57 edwards Exp $
/*! \file
 *  \brief Project a complex Nc x Nc matrix W onto SU(Nc) by maximizing Tr(VW)
 */

#ifndef __sun_proj_h__
#define __sun_proj_h__

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
