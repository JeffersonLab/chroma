// $Id: unprec_qprop_w.cc,v 1.2 2003-04-09 21:14:01 edwards Exp $
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

void UnprecWilsonTypeFermAct::qprop(LatticeFermion& psi, 
				    const multi1d<LatticeColorMatrix>& u, 
				    const LatticeFermion& chi, 
				    const Real& RsdCG, 
				    int MaxCG, int& ncg_had) const
{
  START_CODE("UnprecWilsonTypeFermAct::qprop");

  int n_count;
  
  /* Construct the linear operator */
  /* This allocates field for the appropriate action */
  const LinearOperator* A = linOp(u);

  LatticeFermion tmp;

  switch(InvType)
  {
  case CG_INVERTER: 
    /* chi_1 = M_dag(u) * chi_1 */
    tmp = (*A)(chi, MINUS);
    
    /* psi = (M^dag * M)^(-1) chi */
    InvCG2 (*A, tmp, psi, RsdCG, MaxCG, n_count);
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
    QDP_error_exit("Unknown inverter type", InvType);
  }
  
  if ( n_count == MaxCG )
    QDP_error_exit("no convergence in the inverter", n_count);
  
  ncg_had += n_count;
  
  // Call the virtual destructor of A
  delete A;

  END_CODE("UnprecWilsonTypeFermAct::qprop");
}
