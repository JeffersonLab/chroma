// $Id: fermact_qprop_array.cc,v 1.6 2004-01-02 03:19:41 edwards Exp $
/*! \file
 *  \brief Propagator solver for a generic non-preconditioned fermion operator
 *
 *  Solve for the propagator of a generic non-preconditioned fermion operator
 */

#include "chromabase.h"
#include "fermact.h"
#include "actions/ferm/invert/invcg2_array.h"

using namespace QDP;

//! Propagator of a generic non-preconditioned fermion linear operator
/*! \ingroup qprop
 *
 * This routine is actually generic to all non-preconditioned (not red/black) fermions
 *
 * Compute the lattice fermion for a generic non-red/black fermion
 * using the source in "chi" - so, the source can
 * be of any desired form. The result will appear in "psi", which on input
 * contains an initial guess for the solution.

 * \param psi      quark propagator ( Modify )
 * \param state    gauge field ( Read )
 * \param chi      source ( Read )
 * \param invType  inverter type ( Read (
 * \param RsdCG    CG (or MR) residual used here ( Read )
 * \param MaxCG    maximum number of CG iterations ( Read )
 * \param ncg_had  number of CG iterations ( Write )
 */

template<typename T>
static 
void qprop_t(const FermionAction< multi1d<T> >& me,
	     multi1d<T>& psi, 
	     Handle<const ConnectState> state, 
	     const multi1d<T>& chi, 
	     enum InvType invType,
	     const Real& RsdCG, 
	     int MaxCG, int& ncg_had)
{
  START_CODE("FermionActionArray::qprop");

  int n_count;
  
  /* Construct the linear operator */
  /* This allocates field for the appropriate action */
  Handle<const LinearOperator< multi1d<T> > > A(me.linOp(state));

  switch(invType)
  {
  case CG_INVERTER: 
  {
    /* chi_1 = M_dag(u) * chi_1 */
    multi1d<T> tmp(me.size());
    (*A)(tmp, chi, MINUS);
    
    /* psi = (M^dag * M)^(-1) chi */
    InvCG2 (*A, tmp, psi, RsdCG, MaxCG, n_count);
  }
  break;
  
#if 0
  case MR_INVERTER:
    /* psi = M^(-1) chi */
    InvMR (*A, chi, psi, MRover, RsdCG, MaxCG, n_count);
    break;

  case BICG_INVERTER:
    /* psi = M^(-1) chi */
    InvBiCG (*A, chi, psi, RsdCG, MaxCG, n_count);
    break;
#endif
  
  default:
    QDP_error_exit("Unknown inverter type", invType);
  }
  
  if ( n_count == MaxCG )
    QDP_error_exit("no convergence in the inverter", n_count);
  
  ncg_had = n_count;
  
  END_CODE("FermionActionArray::qprop");
}


template<>
void 
FermionAction< multi1d<LatticeFermion> >::qpropT(multi1d<LatticeFermion>& psi, 
						 Handle<const ConnectState> state, 
						 const multi1d<LatticeFermion>& chi, 
						 enum InvType invType,
						 const Real& RsdCG, 
						 int MaxCG, int& ncg_had) const
{
  qprop_t<LatticeFermion>(*this, psi, state, chi, invType, RsdCG, MaxCG, ncg_had);
}





