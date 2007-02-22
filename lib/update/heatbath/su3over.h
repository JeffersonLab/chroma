// -*- C++ -*-
// $Id: su3over.h,v 3.1 2007-02-22 21:11:49 bjoo Exp $
/*! \file
 *  \brief Do one SU(2) subgroup microcanonical overrelaxation update of SU(Nc)
 */

#ifndef __su3over_h__
#define __su3over_h__

namespace Chroma {

//! Do one SU(2) subgroup microcanonical overrelaxation update of SU(Nc) matrix
/*!
 * \ingroup heatbath
 *
 * Do one SU(2) subgroup microcanonical overrelaxation update of SU(Nc)
 * matrix U keeping action tr(U*W) constant.
 *
 * \param u            field to be updated ( Modify )
 * \param w            "staple" field in the action ( Read )
 * \param su2_index    SU(2) subgroup index ( Read ) 
 * \param sub          Subset for operations ( Read ) 
 */

void su3over(LatticeColorMatrix& u,
	     const LatticeColorMatrix& w,
	     int su2_index,
	     const Subset& sub);

}  // end namespace Chroma


#endif
