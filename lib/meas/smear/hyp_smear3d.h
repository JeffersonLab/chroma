// -*- C++ -*-
// $Id: hyp_smear3d.h,v 3.0 2006-04-03 04:59:05 edwards Exp $

#ifndef __hyp_smear3d_h__
#define __hyp_smear3d_h__

namespace Chroma { 
//! Construct the "hyp-smeared" links of Anna Hasenfratz, with
//! staple coefficients alpha1, alpha2 and alpha3

/*!
 * Arguments:
 *
 *  \param u		gauge field (Read)
 *  \param u_hyp	"hyp-smeared" gauge field (Write)
 *  \param alpha1	staple coefficient "1" (Read, not used)
 *  \param alpha2	staple coefficient "2" (Read)
 *  \param alpha3	staple coefficient "3" (Read)
 *  \param BlkAccu	accuracy in SU(Nc) projection (Read)
 *  \param BlkMax	max number of iterations in SU(Nc) projection (Read)
 *  \param j_decay	direction of no staple(Read)
 */

void Hyp_Smear3d(const multi1d<LatticeColorMatrix>& u,
	         multi1d<LatticeColorMatrix>& u_hyp,
	         const Real& alpha1, const Real& alpha2, const Real& alpha3,
	         const Real& BlkAccu, int BlkMax,int j_decay);

};

#endif
