// $Id: wlfrmnse_w.cc,v 3.0 2006-04-03 04:59:10 edwards Exp $

#error "NOT FULLY CONVERTED"

#include "chromabase.h"


/* This routine is specific to Wilson fermions! */

//! Calculates the Pseudofermion fields and the pseudofermionic action
/*!
 * \ingroup molecdyn
 *
 */
/* u           -- gauge field ( Read ) */
/* chi         -- pseudofermion field. ( Write ) */
/* psi         -- pseudofermion field. ( Write ) */
/* n_congrd    -- Number of CG iteration ( Write ) */

/* eta = Gaussian random noise */
/* chi = M_dag * eta */
/* psi = (M_dag*M)^(-1) * chi */

void WlFrmNse(const multi1d<LatticeColorMatrix>& u,
	      LatticeFermion& chi,
	      LatticeFermion& psi,
	      int& n_congrd,
{
  START_CODE();
  
  LINEAR_OPERATOR(A);

  phfctr(u, FORWARD);              /* ON */

  /* Construct the linear operator */
  /* This allocates field for the appropriate action */
  ConsLinOp(A, u, KappaMC, FermAct);

  /* Initialize psi */
  psi = 0;
  
  /* Fill eta with random gaussian noise: < eta_dagger * eta > = 1 */
  LatticeFermion eta;
  gaussian(eta);
  
#if 0
  /* If using Schroedinger functional, zero out the boundaries */
  if ( SchrFun > 0 )
  {
    FILLMASK(eta, lSFmaskF(1), ZERO);
  }
#endif

  /* chi = M_dag*eta */
  A(A, eta, chi, 1, MINUS);

  /* psi = (1/(M_dag*M))*chi */
  InvCG2(A, chi, psi, RsdCGMC, 1, n_congrd);

  
  /* Clean up from the operator */
  DestLinOp (A, FermAct);
    
  phfctr (u, BACKWARD);              /* OFF */

  END_CODE();
}
