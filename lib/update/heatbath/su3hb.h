// -*- C++ -*-
// $Id: su3hb.h,v 3.1 2007-02-22 21:11:49 bjoo Exp $
/*! \file
 *  \brief Do one SU(2) subgroup heatbath update of SU(Nc) matrix U with action
 */

#ifndef __su3hb_h__
#define __su3hb_h__

namespace Chroma {

  //! Algorithm type
  enum HeatbathType 
    {
      HEATBATH_TYPE_KPHB = 201,
      HEATBATH_TYPE_CrHB = 202,
    };


//! Do one SU(2) subgroup heatbath update of SU(Nc) matrix U with action
/*!
 * \ingroup heatbath
 *
 * BetaMC * [tr(U*W) + hc] / (2Nc).
 * Do at most nheat update trials of the Kennedy-Pendleton or Creutz SU(2)
 * heatbath algorithm.
 *
 * \param u          field to be updated ( Modify )
 * \param w          "staple" field in the action ( Read )
 * \param su2_index  SU(2) subgroup index ( Read )
 * \param nheat      maximal number of update trials ( Read )
 * \param ntrials    total number of link trials ( Write )
 * \param nfails     total number of failed trials ( Write ) 
 * \param sub        Subset for updating ( Read )
 */

void su3hb(LatticeColorMatrix& u,
	   const LatticeColorMatrix& w,
	   int su2_index,
	   int nheat,
	   int& ntrials,
	   int& nfails,
	   HeatbathType algorithm,
	   const Subset& sub);

}  // end namespace Chroma

#endif
