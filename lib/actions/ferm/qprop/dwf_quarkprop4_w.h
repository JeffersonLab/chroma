// $Id: dwf_quarkprop4_w.h,v 1.3 2004-01-30 21:35:49 kostas Exp $
// $Log: dwf_quarkprop4_w.h,v $
// Revision 1.3  2004-01-30 21:35:49  kostas
// added uprec_dwf support
//
// Revision 1.2  2004/01/30 20:21:32  kostas
// fixed the prototype
// 
/*! \file
 * \brief Full quark propagator solver for domain wall fermions
 *
 * Given a complete propagator as a source, this does all the inversions needed
 */

#ifndef __dwf_quarkprop4_w_h__
#define __dwf_quarkprop4_w_h__

#include "fermact.h"
#include "actions/ferm/fermacts/fermacts.h"

using namespace QDP;


//! Given a complete propagator as a source, this does all the inversions needed
/*! \ingroup qprop
 *
 * This routine is actually generic to Domain Wall fermions (Array) fermions
 *
 * \param q_sol    quark propagator ( Write )
 * \param q_src    source ( Read )
 * \param invType  inverter type ( Read (
 * \param RsdCG    CG (or MR) residual used here ( Read )
 * \param MaxCG    maximum number of CG iterations ( Read )
 * \param ncg_had  number of CG iterations ( Write )
 */

void dwf_quarkProp4(LatticePropagator& q_sol, 
		    XMLWriter& xml_out,
		    const LatticePropagator& q_src,
		    const int t_src,
		    const EvenOddPrecDWFermActBaseArray<LatticeFermion>& S_f,
		    Handle<const ConnectState> state,
		    enum InvType invType,
		    const Real& RsdCG, 
		    int MaxCG, int& ncg_had);

//! Given a complete propagator as a source, this does all the inversions needed
/*! \ingroup qprop
 *
 * This routine is actually generic to Domain Wall fermions (Array) fermions
 *
 * \param q_sol    quark propagator ( Write )
 * \param q_src    source ( Read )
 * \param invType  inverter type ( Read (
 * \param RsdCG    CG (or MR) residual used here ( Read )
 * \param MaxCG    maximum number of CG iterations ( Read )
 * \param ncg_had  number of CG iterations ( Write )
 */
void dwf_quarkProp4(LatticePropagator& q_sol, 
		    XMLWriter& xml_out,
		    const LatticePropagator& q_src,
		    const int t_src,
		    const UnprecDWFermActBaseArray<LatticeFermion>& S_f,
		    Handle<const ConnectState> state,
		    enum InvType invType,
		    const Real& RsdCG, 
		    int MaxCG, int& ncg_had) ;

#endif
