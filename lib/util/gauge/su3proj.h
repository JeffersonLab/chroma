// -*- C++ -*-
// $Id: su3proj.h,v 3.1 2007-02-22 21:11:50 bjoo Exp $
/*! \file
 *  \ingroup gauge
 *  \author Subsetting added by A. Hart
 *  \param[in,out] u        the projected SU(3) Matrix (Modify)
 *  \param[in] w            matrix against which to maximize (Read)
 *  \param[in] su2_index    SU(2) subgroup index (Read)
 *  \param[in] mstag        An (un)ordered subset of sites (Read)
 *  \brief Project a GL(3,C) color matrix onto SU(3)
 *
 *  Project a complex Nc x Nc matrix W onto SU(Nc) by maximizing Tr(VW)
 */

#ifndef __su3proj_h__
#define __su3proj_h__

namespace Chroma { 
void su3proj(LatticeColorMatrix& u, 
	     const LatticeColorMatrix& w, 
	     int su2_index);

void su3proj(LatticeColorMatrix& u, 
	     const LatticeColorMatrix& w, 
	     int su2_index,
	     const Subset& mstag);

};
#endif
