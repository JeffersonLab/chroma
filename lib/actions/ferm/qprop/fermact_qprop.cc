// $Id: fermact_qprop.cc,v 1.10 2004-09-08 02:48:26 edwards Exp $
/*! \file
 *  \brief Propagator solver for a generic non-preconditioned fermion operator
 *
 *  Solve for the propagator of a generic non-preconditioned fermion operator
 */

#include "chromabase.h"
#include "fermact.h"
#include "invtype.h"
#include "actions/ferm/invert/invcg2.h"

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
void qprop_t(const FermionAction<T>& me,
	     T& psi, 
	     Handle<const ConnectState> state, 
	     const T& chi, 
	     const InvertParam_t& invParam,
	     int& ncg_had)
{
  START_CODE();

  int n_count;
  
  /* Construct the linear operator */
  /* This allocates field for the appropriate action */
  Handle<const LinearOperator<T> > A(me.linOp(state));

  switch(invParam.invType)
  {
  case CG_INVERTER: 
  {
    /* tmp = M_dag(u) * chi_1 */
    T  tmp;
    (*A)(tmp, chi, MINUS);
    
    /* psi = (M^dag * M)^(-1) chi */
    InvCG2 (*A, tmp, psi, invParam.RsdCG, invParam.MaxCG, n_count);
  }
  break;
  
#if 0
  case MR_INVERTER:
    /* psi = M^(-1) chi */
    InvMR (*A, chi, psi, invParam.MRover, invParam.RsdCG, invParam.MaxCG, n_count);
    break;

  case BICG_INVERTER:
    /* psi = M^(-1) chi */
    InvBiCG (*A, chi, psi, invParam.RsdCG, invParam.MaxCG, n_count);
    break;
#endif
  
  default:
    QDP_error_exit("Unknown inverter type", invParam.invType);
  }
  
  if ( n_count == invParam.MaxCG )
    QDP_error_exit("no convergence in the inverter", n_count);
  
  ncg_had = n_count;
  
  END_CODE();
}


template<>
void FermionAction<LatticeFermion>::qpropT(LatticeFermion& psi, 
					   Handle<const ConnectState> state, 
					   const LatticeFermion& chi, 
					   const InvertParam_t& invParam,
					   int& ncg_had) const
{
  qprop_t(*this, psi, state, chi, invParam, ncg_had);
}

template<>
void FermionAction<LatticeFermion>::qprop(LatticeFermion& psi, 
					  Handle<const ConnectState> state, 
					  const LatticeFermion& chi, 
					  const InvertParam_t& invParam,
					  int& ncg_had) const
{
  qprop_t(*this, psi, state, chi, invParam, ncg_had);
}



