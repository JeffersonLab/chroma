// $Id: overlap_fermact_base_w.cc,v 1.3 2004-01-08 11:53:08 bjoo Exp $
/*! \file
 *  \brief Base class for unpreconditioned overlap-like fermion actions
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/overlap_state.h"
#include "actions/ferm/fermacts/overlap_fermact_base_w.h"
#include "actions/ferm/invert/invcg1.h"
#include "actions/ferm/invert/invcg2.h"


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
			  Handle<const ConnectState> state, 
			  const LatticeFermion& chi, 
			  enum InvType invType,
			  const Real& RsdCG, 
			  int MaxCG, int& ncg_had) const
{
  START_CODE("OverlapFermActBase::qprop");
  Handle< const LinearOperator<LatticeFermion> > M(linOp(state));
  int n_count;

  Real mass = quark_mass();

  switch( invType ) {
  case CG_INVERTER:
    {

      LatticeFermion tmp;
      
      // Check whether the source is chiral.
      Chirality ichiral = isChiralVector(chi);
      if( ichiral == CH_NONE || ( isChiral() == false )) { 
	
	
	(*M)(tmp, chi, MINUS);
      
	// Source is not chiral. In this case we should use,
	// InvCG2 with M
	InvCG2(*M, tmp, psi, RsdCG, MaxCG, n_count);
      }
      else {
	
	// Source is chiral. In this case we should use InvCG1
	// with the special MdagM
	Handle< const LinearOperator<LatticeFermion> > MM(lMdagM(state, ichiral));
	InvCG1(*MM, chi, tmp, RsdCG, MaxCG, n_count);
	(*M)(psi, tmp, MINUS);
      }
    }
    break;

#if 0
  case MR_INVERTER:
    // psi = D^(-1)* chi
    InvMR (*M, chi, psi, MRover, RsdCG, MaxCG, n_count);
    break;

  case BICG_INVERTER:
    // psi = D^(-1) chi
    InvBiCG (*M, chi, psi, RsdCG, MaxCG, n_count);
    break;
#endif
  
  default:
    QDP_error_exit("Zolotarev4DFermActBj::qprop Solver Type not implemented\n");
    break;
  };

  if ( n_count == MaxCG ) { 
    QDP_error_exit("Zolotarev4DFermAct::qprop: No convergence in solver: n_count = %d\n", n_count);
  }

  // Update the ncg_had counter
  ncg_had = n_count;
  
  // Normalize and remove contact term 
  Real ftmp = Real(1) / ( Real(1) - mass );
  
  psi -= chi;
  psi *= ftmp;

  END_CODE("OverlapFermActBase::qprop");
}
  


