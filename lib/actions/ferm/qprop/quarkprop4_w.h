// $Id: quarkprop4_w.h,v 1.8 2004-01-06 20:16:47 edwards Exp $
/*! \file
 *  \brief Full quark propagator solver
 *
 *  Given a complete propagator as a source, this does all the inversions needed
 */

#ifndef __quarkprop4_w_h__
#define __quarkprop4_w_h__

#include "fermact.h"

//! Given a complete propagator as a source, this does all the inversions needed
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

void quarkProp4(LatticePropagator& q_sol, 
		XMLWriter& xml_out,
		const LatticePropagator& q_src,
		const WilsonTypeFermAct<LatticeFermion>& S_f,
		Handle<const ConnectState> state,
		enum InvType invType,
		const Real& RsdCG, 
		int MaxCG, int& ncg_had);

//! Given a complete propagator as a source, this does all the inversions needed
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

void quarkProp4(LatticePropagator& q_sol, 
		XMLWriter& xml_out,
		const LatticePropagator& q_src,
		const WilsonTypeFermAct<LatticeDWFermion>& S_f,
		Handle<const ConnectState> state,
		enum InvType invType,
		const Real& RsdCG, 
		int MaxCG, int& ncg_had);


//! Given a complete propagator as a source, this does all the inversions needed
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

void quarkProp4(LatticePropagator& q_sol, 
		XMLWriter& xml_out,
		const LatticePropagator& q_src,
		const WilsonTypeFermAct< multi1d<LatticeFermion> >& S_f,
		Handle<const ConnectState> state,
		enum InvType invType,
		const Real& RsdCG, 
		int MaxCG, int& ncg_had);


#include "actions/ferm/fermacts/prec_dwf_fermact_base_array_w.h"

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

void quarkProp4(LatticePropagator& q_sol, 
		XMLWriter& xml_out,
		const LatticePropagator& q_src,
		const EvenOddPrecDWFermActBaseArray<LatticeFermion>& S_f,
		Handle<const ConnectState> state,
		enum InvType invType,
		const Real& RsdCG, 
		int MaxCG, int& ncg_had);

#endif
