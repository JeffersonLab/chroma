// $Id: staggered_qprop.cc,v 1.8 2004-10-14 14:19:32 mcneile Exp $
/*! \file
 *  \brief Propagator solver for a generic non-preconditioned fermion operator
 *
 *  Solve for the propagator of a generic non-preconditioned fermion operator
 */

#include "chromabase.h"
#include "fermact.h"
#include "linearop.h"
#include "actions/ferm/invert/invcg1.h"

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
 * \param u        gauge field ( Read )
 * \param chi      source ( Read )
 * \param RsdCG    CG (or MR) residual used here ( Read )
 * \param MaxCG    maximum number of CG iterations ( Read )
 * \param ncg_had  number of CG iterations ( Write )
 */

void 
EvenOddStaggeredTypeFermAct<LatticeFermion>::qprop(LatticeFermion& psi, 
						   Handle<const ConnectState> state,
						   const LatticeFermion& chi,
						   const InvertParam_t& invParam,
						   int& ncg_had)
{
  START_CODE();

  int n_count;
  
  /* Construct the linear operator */
  /* This allocates field for the appropriate action */
  Handle<const EvenOddLinearOperator<LatticeFermion> > M(linOp(state));
  Handle<const LinearOperator<LatticeFermion> > A(lMdagM(state));

  LatticeFermion tmp, tmp1, tmp2;
  tmp = tmp1 = tmp2 = zero;
  Real invm;

  //  switch(invType)
  //{
  //case CG_INVERTER: 

    // Make preconditioned source:  tmp_1_e = M_ee chi_e + M_eo^{dag} chi_o

    M->evenEvenLinOp(tmp, chi, PLUS);
    M->evenOddLinOp(tmp2, chi, MINUS);
    tmp[rb[0]] += tmp2;
    

    /* psi = (M^dag * M)^(-1) chi  = A^{-1} chi*/
    InvCG1(*A, tmp, psi, invParam.RsdCG, invParam.MaxCG, n_count);
   
    // psi[rb[0]] is returned, so reconstruct psi[rb[1]]
    invm = Real(1)/(2*getQuarkMass());
    
    // tmp_1_o = D_oe psi_e 
    M->oddEvenLinOp(tmp1, psi, PLUS);

    // tmp_1_o = (1/2m) D_oe psi_e
    tmp1[rb[1]] *= invm;

    // tmp_2_o = (1/2m) chi_o
    tmp2[rb[1]]  = invm * chi;

    // psi_o = (1/2m) chi_o - (1/2m) D_oe psi_e 
    psi[rb[1]] = tmp2 - tmp1;
    //    break;  

#if 0
  case MR_INVERTER:
    /* psi = M^(-1) chi */
    InvMR (M, chi, psi, MRover, RsdCG, MaxCG, n_count);
    break;

  case BICG_INVERTER:
    /* psi = M^(-1) chi */
    InvBiCG (M, chi, psi, RsdCG, MaxCG, n_count);
    break;
#endif
  
    //  default:
    // QDP_error_exit("Unknown inverter type", invType);
    //}

  if ( n_count == invParam.MaxCG )
    QDP_error_exit("no convergence in the inverter", n_count);

  ncg_had = n_count;
  
  // Call the virtual destructor of A
  // delete A;

  END_CODE();
}


//template<>
//void AsqtadFermAction<LatticeFermion>::qprop(LatticeFermion& psi,
//                                           const 
//multi1d<LatticeColorMatrix>& u_fat,
//					   const 
//multi1d<LatticeColorMatrix>& u_triple,
//                                           const LatticeFermion& chi,
//                                           int invType,
//                                           const Real& RsdCG,
//                                           int MaxCG, int& ncg_had) 
//const
//{   
//  qprop_t(psi, u_fat, u_triple, chi, invType, RsdCG, MaxCG, ncg_had);
//}   

