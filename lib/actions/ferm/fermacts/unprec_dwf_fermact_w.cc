// $Id: unprec_dwf_fermact_w.cc,v 1.2 2003-11-08 04:21:47 edwards Exp $
/*! \file
 *  \brief Unpreconditioned domain-wall fermion action
 */

#include "chromabase.h"
#include "actions/ferm/linop/unprec_dwf_linop_w.h"
#include "actions/ferm/fermacts/unprec_dwf_fermact_w.h"
#include "actions/ferm/linop/lmdagm_w.h"

//! Creation routine
/*! \ingroup fermact
 *
 * \param WilsonMass   DWF height    (Read)
 * \param m_q          quark mass    (Read)
 */
void UnprecDWFermAct::create(const Real& WilsonMass_, const Real& m_q_)
{
  WilsonMass = WilsonMass_;
  m_q = m_q_;
  a5  = 1.0;
//    CoeffWilsr_s = (AnisoP) ? Wilsr_s / xiF_0 : 1;
}

//! Produce a linear operator for this action
/*!
 * \ingroup fermact
 *
 * The operator acts on the entire lattice
 *
 * \param u 	    gauge field     	       (Read)
 */
const LinearOperator<LatticeDWFermion>* 
UnprecDWFermAct::linOp(const multi1d<LatticeColorMatrix>& u) const
{
  return new UnprecDWLinOp(u,WilsonMass,m_q);
}

//! Produce a M^dag.M linear operator for this action
/*!
 * \ingroup fermact
 *
 * The operator acts on the entire lattice
 *
 * \param u 	    gauge field     	       (Read)
 */
const LinearOperator<LatticeDWFermion>* 
UnprecDWFermAct::lMdagM(const multi1d<LatticeColorMatrix>& u) const
{
  LinearOperator<LatticeDWFermion>* mdagm = new lmdagm<LatticeDWFermion>(UnprecDWLinOp(u,WilsonMass,m_q));
  return mdagm;
}



//-------------------------------------------------------------------------------------
#include "actions/ferm/invert/invcg2.h"

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
  
  QDPIO::cout << "|psi|^2 = " << norm2(psi) << endl;
  QDPIO::cout << "|chi|^2 = " << norm2(chi) << endl;

  /* Construct the linear operator */
  /* This allocates field for the appropriate action */
  const LinearOperator<LatticeDWFermion>* A = linOp(u);

  LatticeDWFermion tmp, psi_tmp, chi_tmp;

  psi_tmp = zero;
  pokeDW(psi_tmp, psi, 0);

  chi_tmp = zero;
  pokeDW(chi_tmp, chi, 0);     // WRONG!!!

  QDPIO::cout << "|psi_tmp|^2 = " << norm2(psi_tmp) << endl;
  QDPIO::cout << "|chi_tmp|^2 = " << norm2(chi_tmp) << endl;

  switch(invType)
  {
  case CG_INVERTER: 
    /* chi_tmp = M_dag(u) * chi_tmp */
    tmp = (*A)(chi_tmp, MINUS);
    
    /* psi_tmp = (M^dag * M)^(-1) chi_tmp */
    InvCG2 (*A, tmp, psi_tmp, RsdCG, MaxCG, n_count);
    break;
  
#if 0
  case MR_INVERTER:
    /* psi_tmp = M^(-1) chi_tmp */
    InvMR (*A, chi_tmp, psi_tmp, MRover, RsdCG, MaxCG, n_count);
    break;

  case BICG_INVERTER:
    /* psi_tmp = M^(-1) chi_tmp */
    InvBiCG (*A, chi_tmp, psi_tmp, RsdCG, MaxCG, n_count);
    break;
#endif
  
  default:
    QDP_error_exit("Unknown inverter type", invType);
  }
  
  if ( n_count == MaxCG )
    QDP_error_exit("no convergence in the inverter", n_count);
  
  ncg_had = n_count;
  
  psi = peekDW(psi_tmp, 0);    // WRONG!!!

  // Call the virtual destructor of A
  delete A;

  END_CODE("UnprecDWTypeFermAct::qprop");
}
