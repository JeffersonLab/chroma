// $Id: overlap_fermact_qprop_w.cc,v 1.2 2004-05-12 15:45:10 bjoo Exp $
/*! \file
 *  \brief Base class for unpreconditioned overlap-like fermion actions
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/overlap_fermact_base_w.h"
#include "actions/ferm/invert/invcg2.h"
#include "actions/ferm/invert/invsumr.h"

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


template<typename T>
static 
void qprop_t(const OverlapFermActBase<T>& me,
	     T& psi,
	     Handle<const ConnectState> state_, 
	     const T& chi, 
	     enum InvType invType,
	     const Real& RsdCG, 
	     int MaxCG, int& ncg_had) const
{
  START_CODE("OverlapFermActBase::qprop");

  const EVConnectState<T>& state = dynamic_cast<const EVConnectState<T>&>(*state_);

  const Real m_q = me.quark_mass();

  int n_count;
  

  switch(invType)
  {
  case CG_INVERTER: 
  {

    // Construct the linear operator
    Handle<const LinearOperator<T> > A(me.linOp(state_));

    T tmp1;

#if 1
    // For the moment, use simple minded inversion
    (*A)(tmp1, chi, MINUS);

    InvCG2(*A, tmp1, psi, RsdCG, MaxCG, n_count);

#else
    /* psi = (D^dag * D)^(-1) * D^dag * chi */
    /*     = D^dag * (D * D^dag)^(-1) * chi */
    /*     = D^(-1) * chi */

    /* Check if source is chiral */
    int ichiral = me.ischiral(chi);

    if (ichiral == 0 || me.isChiral())  // also check if action is not chiral
    {
      // Source or action is not chiral: call the CG2 wrapper
      (*A)(tmp1, chi, MINUS);

      InvCG2(*A, tmp1, psi, RsdCG, MaxCG, n_count);
    }
    else
    {
      // Source is chiral: use the CG1
      // Construct the linear operator
      Handle<const LinearOperator<T> > B(me.lMdagM(state_));  // OOPS, NEED TO PASS ICHIRAL!!!!

      InvCG1(*B, chi, tmp1, RsdCG, MaxCG, n_count);
      (*A)(psi, tmp1, MINUS);
    }
#endif
  }
  break;
  
#if 0
  case MR_INVERTER:
    {
      // Construct the linear operator
      Handle<const LinearOperator<T> > A(me.linOp(state_));
      
      // psi = D^(-1)* chi
      InvMR (*A, chi, psi, MRover, RsdCG, MaxCG, n_count);
    }
    break;
      
  case BICG_INVERTER:
    // psi = D^(-1) chi
    {
      // Construct the linear operator
      Handle<const LinearOperator<T> > A(me.linOp(state_));

      InvBiCG (*A, chi, psi, RsdCG, MaxCG, n_count);
    }
    break;
#endif

  case SUMR_INVERTER:
    {
      // Solve by SUMR solver -- for shifted unitary matrices
      //
      // Solve zeta I + rho gamma_5 eps(H)
      // where gamma_5 eps(H) is unitary
      //
      // zeta = (1 + mu)/(1-mu)
      // rho  = 1
      Real rho = Real(1);
      Real mu = me.quark_mass();
      Complex zeta = cmplx(( Real(1) + mu ) / (Real(1) - mu),0);
      {
	Handle<const LinearOperator<T> > U(me.lgamma5epsH(state_));
	
	// Now solve:
	InvSUMR(*U, chi, psi, zeta, rho, RsdCG, MaxCG, n_count);
	
	// Restore to normal scaling
	Real fact = Real(2)/(Real(1) - mu);
	psi *= fact;
	
	// Check back
      }

      // Get a proper operator 
      Handle<const LinearOperator<T> > A(me.linOp(state_));
      LatticeFermion Dpsi;
      (*A)(Dpsi, psi, PLUS);
      Dpsi = chi - Dpsi;
      Dpsi /= sqrt(norm2(chi));
      QDPIO::cout << "OvQprop || chi - D psi || = " << sqrt(norm2(Dpsi))
		  << "  n_count = " << n_count << " iters" << endl;
    }
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
  psi -= chi;
  psi *= ftmp1;

  END_CODE("OverlapFermActBase::qprop");
}


template<>
void OverlapFermActBase<LatticeFermion>::qprop(LatticeFermion& psi, 
					       Handle<const ConnectState> state, 
					       const LatticeFermion& chi, 
					       enum InvType invType,
					       const Real& RsdCG, 
					       int MaxCG, int& ncg_had) const
{
  qprop_t<LatticeFermion>(*this, psi, state, chi, invType, RsdCG, MaxCG, ncg_had);
}
