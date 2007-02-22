// -*- C++ -*-
// $Id: sun_proj.h,v 3.2 2007-02-22 21:11:50 bjoo Exp $
/*! \file
 *  \ingroup gauge
 *  \author Subsetting added by A. Hart
 *  \param[in] w            complex Nc x Nc matrix
 *  \param[out] v           the projected SU(Nc) Matrix
 *  \param[in] BlkAccu      accuracy in SU(Nc) projection
 *  \param[in] BlkMax       max number of iterations in SU(Nc) projection
 *  \param[in] mstag        an (un)ordered subset of lattice sites
 *  \brief Project a complex Nc x Nc matrix W onto SU(Nc) by maximizing Tr(VW)
 *
 *  Project a complex Nc x Nc matrix W onto SU(Nc) by maximizing Tr(VW)
 */

#ifndef __sun_proj_h__
#define __sun_proj_h__

namespace Chroma 
{ 

  // No subsets 
  void sun_proj(const LatticeColorMatrix& w, 
		LatticeColorMatrix& v,
		const Real& BlkAccu, 
		int BlkMax);

  // Ordered subsets
  void sun_proj(const LatticeColorMatrix& w, 
		LatticeColorMatrix& v,
		const Real& BlkAccu, 
		int BlkMax,
		const Subset& mstag);


} // End namespace
#endif
