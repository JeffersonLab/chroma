// $Id: nef_quarkprop4_w.h,v 1.2 2004-09-03 14:24:36 kostas Exp $
// $Log: nef_quarkprop4_w.h,v $
// Revision 1.2  2004-09-03 14:24:36  kostas
// added mres measurement for NEF fermions
//
/*! \file
 * \brief Full quark propagator solver for NEF domain wall fermions
 *
 * Given a complete propagator as a source, this does all the inversions needed
 */

#ifndef __nef_quarkprop4_w_h__
#define __nef_quarkprop4_w_h__

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

void nef_quarkProp4(LatticePropagator& q_sol, 
		    XMLWriter& xml_out,
		    const LatticePropagator& q_src,
		    int t_src, int j_decay,
		    const EvenOddPrecNEFDWFermActArray<LatticeFermion>& S_f,
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

void nef_quarkProp4(LatticePropagator& q_sol, 
		    XMLWriter& xml_out,
		    const LatticePropagator& q_src,
		    int t_src, int j_decay,
		    const UnprecNEFDWFermActArray<LatticeFermion>& S_f,
		    Handle<const ConnectState> state,
		    enum InvType invType,
		    const Real& RsdCG, 
		    int MaxCG, int& ncg_had);

#endif
