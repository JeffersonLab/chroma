// $Id: unprec_ovext_fermact_array_w.cc,v 1.7 2004-03-29 21:32:28 edwards Exp $
/*! \file
 *  \brief Unpreconditioned extended-Overlap (5D) (Naryanan&Neuberger) action
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/unprec_ovext_fermact_array_w.h"
#include "actions/ferm/linop/unprec_ovext_linop_array_w.h"
#include "actions/ferm/linop/lmdagm.h"

#include "actions/ferm/invert/invcg2_array.h"

//! Produce a linear operator for this action
/*!
 * The operator acts on the entire lattice
 *
 * \param state	    gauge field     	       (Read)
 */
const LinearOperator<multi1d<LatticeFermion> >* 
UnprecOvExtFermActArray::linOp(Handle<const ConnectState> state) const
{
  return new UnprecOvExtLinOpArray(state->getLinks(),WilsonMass,m_q,N5);
}


//! Produce a M^dag.M linear operator for this action
/*!
 * The operator acts on the entire lattice
 *
 * \param state	    gauge field     	       (Read)
 */
const LinearOperator<multi1d<LatticeFermion> >* 
UnprecOvExtFermActArray::lMdagM(Handle<const ConnectState> state) const
{
  return new lmdagm<multi1d<LatticeFermion> >(linOp(state));
}

//! Propagator of an un-preconditioned Extended-Overlap linear operator
/*!
 * \param psi      quark propagator ( Modify )
 * \param state    gauge field ( Read )
 * \param chi      source ( Read )
 * \param invType  inverter type ( Read (
 * \param RsdCG    CG (or MR) residual used here ( Read )
 * \param MaxCG    maximum number of CG iterations ( Read )
 * \param ncg_had  number of CG iterations ( Write )
 */

void 
UnprecOvExtFermActArray::qprop(LatticeFermion& psi, 
			       Handle<const ConnectState> state, 
			       const LatticeFermion& chi, 
			       enum InvType invType,
			       const Real& RsdCG, 
			       int MaxCG, int& ncg_had) const
{
  START_CODE("UnprecOvExtFermActArray::qprop");

  const int  N5 = size();   // array size better match
  const Real m_q = quark_mass();
  int n_count;
  
  int G5 = Ns*Ns - 1;

  // Initialize the 5D fields
  multi1d<LatticeFermion> chi5(N5);
  multi1d<LatticeFermion> psi5(N5);
  psi5 = zero;
  chi5 = zero;

  psi5[0] = psi;
  chi5[0] = Gamma(G5) * chi;

  // Construct the linear operator
  Handle<const LinearOperator< multi1d<LatticeFermion> > > A(linOp(state));

  switch(invType)
  {
  case CG_INVERTER: 
    // psi5 = (H_o)^(-2) chi5
    InvCG2(*A, chi5, psi5, RsdCG, MaxCG, n_count);

    // chi5 = H_o * (H_o)^(-2) * gamma_5 * chi
    (*A)(chi5, psi5, MINUS);
    break;
  
  case MR_INVERTER:
  case BICG_INVERTER:
    QDP_error_exit("Unsupported inverter type", invType);
    break;
  
  default:
    QDP_error_exit("Unknown inverter type", invType);
  }
  
  if ( n_count == MaxCG )
    QDP_error_exit("no convergence in the inverter", n_count);
  
  ncg_had = n_count;
  
  // Overall normalization
  Real ftmp1 = Real(1) / Real(1 - m_q);

  // Normalize and remove contact term
  psi = ftmp1*(chi5[0] - chi);

  END_CODE("UnprecOvExtFermActArray::qprop");
}
