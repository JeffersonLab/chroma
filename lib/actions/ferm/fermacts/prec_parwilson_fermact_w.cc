// $Id: prec_parwilson_fermact_w.cc,v 1.2 2004-07-28 02:38:01 edwards Exp $
/*! \file
 *  \brief Even-odd preconditioned Wilson fermion action with parity breaking term
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/prec_parwilson_fermact_w.h"
#include "actions/ferm/linop/prec_parwilson_linop_w.h"
#include "actions/ferm/linop/lmdagm.h"



//! Produce a linear operator for this action
/*!
 * The operator acts on the odd subset
 *
 * \param state 	    gauge field     	       (Read)
 */
const EvenOddPrecLinearOperator<LatticeFermion>* 
EvenOddPrecParWilsonFermAct::linOp(Handle<const ConnectState> state) const
{
  return new EvenOddPrecParWilsonLinOp(state->getLinks(),Mass,H);
}

//! Produce a M^dag.M linear operator for this action
/*!
 * The operator acts on the odd subset
 *
 * \param state 	    gauge field     	       (Read)
 */
const LinearOperator<LatticeFermion>* 
EvenOddPrecParWilsonFermAct::lMdagM(Handle<const ConnectState> state) const
{
  return new lmdagm<LatticeFermion>(linOp(state));
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
 * \param state    gauge field ( Read )
 * \param psi      solution to linear system ( Read )
 */

void
EvenOddPrecParWilsonFermAct::dsdu(multi1d<LatticeColorMatrix>& ds_u,
				  Handle<const ConnectState> state,
				  const LatticeFermion& psi) const
{
  START_CODE();
  
#if 1
  QDPIO::cerr << "EvenOddPrecParWilsonFermAct::dsdu not implemented" << endl;
  QDP_abort(1);
#else

#error "Not quite correct implementation"

  const multi1d<LatticeColorMatrix>& u = state->getLinks();
				 
  LatticeColorMatrix utmp_1;
  LatticeFermion phi;
  LatticeFermion rho;
  LatticeFermion sigma;
  LatticeFermion ftmp_1;
  LatticeFermion ftmp_2;
  Double ddummy;
  Real dummy;
  int nu;
  int cb;
  
  /* Do the Wilson fermion dS_f/dU with parity breaking term */
    
//  CONSTRUCT_LINEAR_OPERATOR(A, lwlhmpsi, u, KappaMD);

  /*  phi = M(u)*psi  */
  A (A, psi, phi, 1, PLUS);
    
  /* rho = (1-i H gamma_5) * Dslash(0<-1) * psi */
  dslash (u, psi, ftmp_1, PLUS, 1);
  PARBREAK(ftmp_1, H_parity, rho, MINUS);
    
  /* phi = (KappaMD^2)*phi/(1+h^2) = -(KappaMD^2)*M*psi/(1+h^2) */
  dummy = -(KappaMD*KappaMD) / (WORD_VALUE(WORD_H_parity, ONE) + H_parity*H_parity);
  phi = phi * dummy;
    
  /* sigma = (1+i H gamma_5) * Dslash_dag(0<-1) * phi */
  dslash (u, phi, ftmp_1, MINUS, 1);
  PARBREAK(ftmp_1, H_parity, sigma, PLUS);
    
        
  for(int mu = 0; mu < Nd; ++mu)
  {
    cb = 0;

    /* ftmp_2 = (gamma(mu))*psi */
    SPIN_PRODUCT(psi,(INTEGER_LSHIFT_FUNCTION(1,mu)),ftmp_2);
    /* ftmp_2(x) = -(psi(x) - ftmp_2(x)) = -(1 - gamma(mu))*psi( x )  */
    ftmp_2 -= psi;
    utmp_1 = -(shift(ftmp_2, cb, FORWARD, mu) * adj(sigma));

    /* ftmp_2 = (gamma(mu))*phi */
    SPIN_PRODUCT(phi,(INTEGER_LSHIFT_FUNCTION(1,mu)),ftmp_2);
    /* ftmp_2 = phi + ftmp_2 = (1 + gamma(mu))*phi( x)  */
    ftmp_2 += phi;
    utmp_1 += shift(ftmp_2, cb, FORWARD, mu) * adj(rho);
    ds_u[mu][cb] += u[mu][cb] * utmp_1;
      
    cb = 1;

    /* ftmp_2 = (gamma(mu))*ftmp_1 */
    SPIN_PRODUCT(rho,(INTEGER_LSHIFT_FUNCTION(1,mu)),ftmp_2);
    /* ftmp_2 = -(rho - ftmp_2) = -(1 - gamma(mu))*rho( x )  */
    ftmp_2 -= rho;
    utmp_1 = -(shift(ftmp_2, cb, FORWARD, mu) * adj(phi));
      
    /* ftmp_2 = (gamma(mu))*sigma */
    SPIN_PRODUCT(sigma,(INTEGER_LSHIFT_FUNCTION(1,mu)),ftmp_2);
    /* ftmp_1 = ftmp_1 + ftmp_2 = (1 + gamma(mu))*sigma( x + mu)  */
    ftmp_2 += sigma;
    utmp_1 += shift(ftmp_2, cb, FORWARD, mu) * adj(psi);
    ds_u[mu][cb] += u[mu][cb] * utmp_1;
      
  }
        
#endif

  END_CODE();
}
