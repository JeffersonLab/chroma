// $Id: unprec_ovext_fermact_array_w.cc,v 1.2 2003-11-20 05:43:41 edwards Exp $
/*! \file
 *  \brief Unpreconditioned extended-Overlap (5D) (Naryanan&Neuberger) action
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/unprec_ovext_fermact_array_w.h"
#include "actions/ferm/linop/unprec_ovext_linop_array_w.h"
#include "actions/ferm/linop/lmdagm_w.h"

#include "actions/ferm/invert/invcg2_array.h"

//! Creation routine
/*!
 * \param WilsonMass_   DWF height    (Read)
 * \param m_q_          quark mass    (Read)
 * \param N5_           extent of flavor space   (Read)
 */
void UnprecOvExtFermActArray::create(const Real& WilsonMass_, const Real& m_q_, int N5_)
{
  WilsonMass = WilsonMass_;
  m_q = m_q_;
  N5  = N5_;

  a5  = 1.0;
}


//! Produce a linear operator for this action
/*!
 * The operator acts on the entire lattice
 *
 * \param u 	    gauge field     	       (Read)
 */
const LinearOperator<multi1d<LatticeFermion> >* 
UnprecOvExtFermActArray::linOp(const multi1d<LatticeColorMatrix>& u) const
{
  return new UnprecOvExtLinOpArray(u,WilsonMass,m_q,N5);
}


//! Produce a M^dag.M linear operator for this action
/*!
 * \ingroup fermact
 *
 * The operator acts on the entire lattice
 *
 * \param u 	    gauge field     	       (Read)
 */
const LinearOperator<multi1d<LatticeFermion> >* 
UnprecOvExtFermActArray::lMdagM(const multi1d<LatticeColorMatrix>& u) const
{
  LinearOperator<multi1d<LatticeFermion> >* mdagm = 
    new lmdagm<multi1d<LatticeFermion> >(UnprecOvExtLinOpArray(u,WilsonMass,m_q,N5));
  return mdagm;
}

//! Propagator of an un-preconditioned Extended-Overlap linear operator
/*!
 * \param psi      quark propagator ( Modify )
 * \param u        gauge field ( Read )
 * \param chi      source ( Read )
 * \param invType  inverter type ( Read (
 * \param RsdCG    CG (or MR) residual used here ( Read )
 * \param MaxCG    maximum number of CG iterations ( Read )
 * \param ncg_had  number of CG iterations ( Write )
 */

void 
UnprecOvExtFermActArray::qprop(LatticeFermion& psi, 
			       const multi1d<LatticeColorMatrix>& u, 
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
  const LinearOperator< multi1d<LatticeFermion> >* A = linOp(u);

  switch(invType)
  {
  case CG_INVERTER: 
    // psi5 = (H_o)^(-2) chi5
    InvCG2 (*A, chi5, psi5, RsdCG, MaxCG, n_count);

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

  // Call the virtual destructor of A
  delete A;

  END_CODE("UnprecOvExtFermActArray::qprop");
}


//! Computes the derivative of the fermionic action respect to the link field
/*!
 *         |  dS      dS_f
 * ds_u -- | ----   + -----   ( Write )
 *         |  dU       dU
 *
 * psi -- [1./(M_dag*M)]*chi_  ( read ) 
 *
 * \param ds_u     result      ( Write )
 * \param u        gauge field ( Read )
 * \param psi      solution to linear system ( Read )
 */

void
UnprecOvExtFermActArray::dsdu(multi1d<LatticeColorMatrix>& ds_u,
			      const multi1d<LatticeColorMatrix>& u, 
			      const multi1d<LatticeFermion>& psi) const
{
  START_CODE("UnprecWilsonFermAct::dsdu");
  
//  multi1d<LatticeColorMatrix> ds_u(Nd);

  ds_u = 0;

  QDPIO::cerr << "UnprecWilsonFermAct::dsdu not implemented" << endl;
  QDP_abort(1);

  END_CODE("UnprecWilsonFermAct::dsdu");
}
