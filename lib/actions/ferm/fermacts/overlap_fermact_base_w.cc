// $Id: overlap_fermact_base_w.cc,v 1.1 2003-12-04 02:56:53 edwards Exp $
/*! \file
 *  \brief Base class for unpreconditioned overlap-like fermion actions
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/overlap_fermact_base_w.h"
#include "actions/ferm/invert/invcg2.h"
#include "actions/ferm/fermacts/ev_state.h"

using namespace QDP;

//! Propagator for unpreconditioned overlap-like fermion actions
/*!
 * \param psi      quark propagator ( Modify )
 * \param state_   gauge field ( Read )
 * \param chi      source ( Read )
 * \param invType  inverter type ( Read (
 * \param RsdCG    CG (or MR) residual used here ( Read )
 * \param MaxCG    maximum number of CG iterations ( Read )
 * \param ncg_had  number of CG iterations ( Write )
 */

void 
OverlapFermActBase::qprop(LatticeFermion& psi, 
			  const ConnectState& state_, 
			  const LatticeFermion& chi, 
			  enum InvType invType,
			  const Real& RsdCG, 
			  int MaxCG, int& ncg_had) const
{
  START_CODE("OverlapFermActBase::qprop");

  const EVConnectState<LatticeFermion>& state = dynamic_cast<const EVConnectState<LatticeFermion>&>(state_);

  const Real m_q = quark_mass();

  int n_count;
  
  // Construct the linear operator
  const LinearOperatorProxy<LatticeFermion> A(linOp(state));

  switch(invType)
  {
  case CG_INVERTER: 
  {
    LatticeFermion tmp1;

#if 1
    // For the moment, use simple minded inversion
    A(tmp1, chi, MINUS);

    InvCG2(A, tmp1, psi, RsdCG, MaxCG, n_count);

#else
    /* psi = (D^dag * D)^(-1) * D^dag * chi */
    /*     = D^dag * (D * D^dag)^(-1) * chi */
    /*     = D^(-1) * chi */

    /* Check if source is chiral */
    int ichiral = ischiral(chi);

    if (ichiral == 0 || this->isChiral())  // also check if action is not chiral
    {
      // Source or action is not chiral: call the CG2 wrapper
      A(tmp1, chi, MINUS);

      InvCG2(A, tmp1, psi, RsdCG, MaxCG, n_count);
    }
    else
    {
      // Source is chiral: use the CG1
      // Construct the linear operator
      const LinearOperatorProxy<LatticeFermion> B(lMdagM(state));  // OOPS, NEED TO PASS ICHIRAL!!!!

      InvCG1(B, chi, tmp1, RsdCG, MaxCG, n_count);
      A(psi, tmp1, MINUS);
    }
#endif
  }
  break;
  
#if 0
  case MR_INVERTER:
    // psi = D^(-1)* chi
    InvMR (A, chi, psi, MRover, RsdCG, MaxCG, n_count);
    break;

  case BICG_INVERTER:
    // psi = D^(-1) chi
    InvBiCG (A, chi, psi, RsdCG, MaxCG, n_count);
    break;
#endif
  
  default:
    QDP_error_exit("Unknown inverter type", invType);
  }
  
  if ( n_count == MaxCG )
    QDP_error_exit("no convergence in the inverter", n_count);
  
  ncg_had = n_count;
  
  // Overall normalization
  Real ftmp1 = Real(1) / Real(1 - m_q);

  // Normalize and remove contact term
  psi -= chi;
  psi *= ftmp1;

  END_CODE("OverlapFermActBase::qprop");
}

