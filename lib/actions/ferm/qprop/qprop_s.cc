// $Id: qprop_s.cc,v 1.1 2003-12-10 12:38:14 bjoo Exp $
/*! \file
 *  \brief Propagator solver for a generic non-preconditioned fermion operator
 *
 *  Solve for the propagator of a generic non-preconditioned fermion operator
 */

#include "chromabase.h"
#include "fermact.h"
#include "primitives.h"
#include "common_declarations.h"
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
 * \param u        gauge field ( Read )
 * \param chi      source ( Read )
 * \param RsdCG    CG (or MR) residual used here ( Read )
 * \param MaxCG    maximum number of CG iterations ( Read )
 * \param ncg_had  number of CG iterations ( Write )
 */

void AsqtadFermTypeAction<LatticeFermion>::qprop(LatticeFermion& psi, 
		                                const multi1d<LatticeColorMatrix>& u_fat,
			  			const multi1d<LatticeColorMatrix>& u_triple, 
	        	  			const LatticeFermion& chi, 
			  			int invType,
		          			const Real& RsdCG, 
			  			int MaxCG, const Real& mass, int& ncg_had) const
{
  START_CODE("AsqtadFermTypeAction::qprop");

  int n_count;
  
  /* Construct the linear operator */
  /* This allocates field for the appropriate action */
  const LinearOperator<LatticeFermion>* M = linOp(u_fat,u_triple);
  const LinearOperator<LatticeFermion>* A = lMdagM(u_fat, u_triple);

  LatticeFermion tmp, tmp1, tmp2;
  Real invm;

  XMLFileWriter xml_out("output3.xml");
//  push(xml_out, "more_tests");
//  Write(xml_out, u_fat);
//  Write(xml_out, u_triple);
//  pop(xml_out);

  switch(invType)
  {
  case CG_INVERTER: 
    /* chi_1 = M_dag(u) * chi_1 */
    (*M)(tmp, chi, MINUS);
   
//    tmp1[rb[0]] = tmp;

    /* psi = (M^dag * M)^(-1) chi */
    InvCG2 (*A, tmp, psi, RsdCG, MaxCG, n_count);
   
// psi[rb[0]] is returned, so reconstruct psi[rb[1]]

    invm = Real(1)/(2*mass);

    (*M)(tmp1, psi, PLUS);
    tmp1 *= invm;
    tmp2 = invm * chi;

    psi[rb[1]] = tmp2 - tmp1;

    push(xml_out, "psi");
    Write(xml_out, psi);
    pop(xml_out);
    
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
  
  // Call the virtual destructor of A
  delete A;

  END_CODE("AsqtadFermTypeAction::qprop");
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

