// $Id: unprec_dwf_qprop_w.cc,v 1.3 2003-11-13 18:19:43 edwards Exp $
/*! \file
 *  \brief Unpreconditioned domain-wall fermion propagator solver
 *
 * Propagator of an un-preconditioned DWF linear operator
 * The conventions used here are specified in 
 * Phys.Rev.D63:094505,2001 (hep-lat/0005002). 
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/unprec_dwf_fermact_w.h"
#include "actions/ferm/invert/invcg2.h"
#include "actions/ferm/linop/dwffld_w.h"

using namespace QDP;

//! Propagator of an un-preconditioned DWF linear operator
/*! \ingroup qprop
 *
 * \param psi      quark propagator ( Modify )
 * \param u        gauge field ( Read )
 * \param chi      source ( Read )
 * \param invType  inverter type ( Read (
 * \param RsdCG    CG (or MR) residual used here ( Read )
 * \param MaxCG    maximum number of CG iterations ( Read )
 * \param ncg_had  number of CG iterations ( Write )
 */

void UnprecDWFermAct::qprop(LatticeFermion& psi, 
			    const multi1d<LatticeColorMatrix>& u, 
			    const LatticeFermion& chi, 
			    enum InvType invType,
			    const Real& RsdCG, 
			    int MaxCG, int& ncg_had) const
{
  START_CODE("UnprecDWTypeFermAct::qprop");

  int n_count;
  
  // Initialize the 5D fields
  //  tmp5 = (chi,0,0,0,..,0)^T
  LatticeDWFermion chi5, tmp5 = zero;
  pokeDW(tmp5, chi, 0);

  // chi5 = P . tmp5
  DwfFld(chi5, tmp5, PLUS);

  // tmp5 = D5(1) . chi5 =  D5(1) . P . (chi,0,0,..,0)^T 
  {
    // Create a Pauli-Villars linop and use it for just this part
    const LinearOperator<LatticeDWFermion>* B = linOpPV(u);

    tmp5 = (*B)(chi5, MINUS);

    delete B;
  }

  //  psi5 = (psi,0,0,0,...,0)^T
  LatticeDWFermion psi5 = zero;
  pokeDW(psi5, psi, 0);

  QDPIO::cout << "|psi5|^2 = " << norm2(psi5) << endl;
  QDPIO::cout << "|chi5|^2 = " << norm2(chi5) << endl;

  // Construct the linear operator
  const LinearOperator<LatticeDWFermion>* A = linOp(u);

  switch(invType)
  {
  case CG_INVERTER: 
    // chi5 = D5^\dagger(m) . tmp5 =  D5^dagger(m) . D5(1) . P . (chi,0,0,..,0)^T
    chi5 = (*A)(tmp5, MINUS);
    
    // psi5 = (D^dag * D)^(-1) chi5
    InvCG2 (*A, chi5, psi5, RsdCG, MaxCG, n_count);
    break;
  
#if 0
  case MR_INVERTER:
    // psi5 = D^(-1) * tmp5
    InvMR (*A, tmp5, psi5, MRover, RsdCG, MaxCG, n_count);
    break;

  case BICG_INVERTER:
    // psi5 = D^(-1) tmp5
    InvBiCG (*A, tmp5, psi5, RsdCG, MaxCG, n_count);
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

  // Project out first slice after  tmp5 <- P^(-1) . psi5
  DwfFld(tmp5, psi5, MINUS);

  // Normalize and remove contact term
  psi = ftmp1*(peekDW(tmp5, 0) - chi);

  // Call the virtual destructor of A
  delete A;

  END_CODE("UnprecDWTypeFermAct::qprop");
}
