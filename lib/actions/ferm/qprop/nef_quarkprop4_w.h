// $Id: nef_quarkprop4_w.h,v 1.1 2004-09-02 22:39:25 kostas Exp $
// $Log: nef_quarkprop4_w.h,v $
// Revision 1.1  2004-09-02 22:39:25  kostas
// started adding the mres/Za measurement for NEF and (OvDWF)
//
// Revision 1.4  2004/02/23 03:05:11  edwards
// Pass in j_decay.
//
// Revision 1.3  2004/01/30 21:35:49  kostas
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
 * \param t_src    time slice of source ( Read )
 * \param j_decay  direction of decay ( Read )
 * \param invType  inverter type ( Read )
 * \param RsdCG    CG (or MR) residual used here ( Read )
 * \param MaxCG    maximum number of CG iterations ( Read )
 * \param ncg_had  number of CG iterations ( Write )
 */

void dwf_quarkProp4(LatticePropagator& q_sol, 
		    XMLWriter& xml_out,
		    const LatticePropagator& q_src,
		    int t_src, int j_decay,
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
 * \param t_src    time slice of source ( Read )
 * \param j_decay  direction of decay ( Read )
 * \param invType  inverter type ( Read )
 * \param RsdCG    CG (or MR) residual used here ( Read )
 * \param MaxCG    maximum number of CG iterations ( Read )
 * \param ncg_had  number of CG iterations ( Write )
 */

void dwf_quarkProp4(LatticePropagator& q_sol, 
		    XMLWriter& xml_out,
		    const LatticePropagator& q_src,
		    int t_src, int j_decay,
		    const UnprecDWFermActBaseArray<LatticeFermion>& S_f,
		    Handle<const ConnectState> state,
		    enum InvType invType,
		    const Real& RsdCG, 
		    int MaxCG, int& ncg_had);

#endif
