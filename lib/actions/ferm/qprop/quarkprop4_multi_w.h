#ifndef __quarkprop4_multi_w_h__
#define __quarkprop4_multi_w_h__

#include "fermact.h"

/*! \ingroup qprop
 *
 * This routine is actually generic to all Wilson-like fermions
 *
 * \param q_sol    quark propagator ( Write )
 * \param q_src    source ( Read )
 * \param invType  inverter type ( Read (
 * \param RsdCG    CG (or MR) residual used here ( Read )
 * \param MaxCG    maximum number of CG iterations ( Read )
 * \param ncg_had  number of CG iterations ( Write )
 */


#include "actions/ferm/fermacts/overlap_fermact_base_w.h"

void multiQuarkProp4(multi1d<LatticePropagator>& q_sol, 
		     XMLWriter& xml_out,
		     const LatticePropagator& q_src,
		     const OverlapFermActBase& S_f,
		     Handle<const ConnectState> state,
		     enum InvType invType,
		     const multi1d<Real>& masses,
		     const multi1d<Real>& RsdCG, 
		     int MaxCG, int& ncg_had);
#endif
