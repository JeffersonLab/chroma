// $Id: prec_wilson_fermact_w.cc,v 1.10 2004-08-05 15:56:25 bjoo Exp $
/*! \file
 *  \brief Even-odd preconditioned Wilson fermion action
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/prec_wilson_fermact_w.h"
#include "actions/ferm/linop/prec_wilson_linop_w.h"
#include "actions/ferm/linop/lmdagm.h"


//! Produce a linear operator for this action
/*!
 * The operator acts on the odd subset
 *
 * \param state 	    gauge field     	       (Read)
 */
const EvenOddPrecLinearOperator<LatticeFermion>* 
EvenOddPrecWilsonFermAct::linOp(Handle<const ConnectState> state) const
{
  const EvenOddPrecLinearOperator<LatticeFermion>* foo;

  if (aniso.anisoP)
    foo = new EvenOddPrecWilsonLinOp(state->getLinks(),Mass,aniso);
  else
    foo = new EvenOddPrecWilsonLinOp(state->getLinks(),Mass);

  return foo;
}

//! Produce a M^dag.M linear operator for this action
/*!
 * The operator acts on the odd subset
 *
 * \param state 	    gauge field     	       (Read)
 */
const LinearOperator<LatticeFermion>* 
EvenOddPrecWilsonFermAct::lMdagM(Handle<const ConnectState> state) const
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
 *
 * In this function I assume that ds_u may already have the gauge piece in there...
 */

void
EvenOddPrecWilsonFermAct::dsdu(multi1d<LatticeColorMatrix>& ds_u,
			       Handle<const ConnectState> state,
			       const LatticeFermion& psi) const
{
  START_CODE();
  
  if (aniso.anisoP)
  {
    QDPIO::cerr << "Currently do not support anisotropy" << endl;
    QDP_abort(1);
  }

  Real prefactor = -Real(1)/(Real(4)*(Real(Nd) + Mass));

				 
  LatticeColorMatrix utmp_1;
  LatticeFermion phi;
  LatticeFermion rho;
  LatticeFermion sigma;

  LatticeFermion ftmp_2;

  // Do the usual Wilson fermion dS_f/dU
  // const LinearOperatorProxy<LatticeFermion> A(linOp(u));
  const Handle< const LinearOperator<LatticeFermion> >&  M(linOp(state));

  // Need the wilson dslash
  // Use u from state with BC's on 
  const multi1d<LatticeColorMatrix>& u = state->getLinks();
  WilsonDslash  D(u);
  //  phi = M(u)*psi
  (*M)(phi, psi, PLUS);
    
  /* rho = Dslash(0<-1) * psi */
  D.apply(rho, psi, PLUS, 1);
 
  /* sigma = Dslash_dag(0 <- 1) * phi */
  D.apply(sigma, phi, MINUS, 1);
    
  for(int mu = 0; mu < Nd; ++mu)
  {

    // ftmp_2(x) = -(psi(x) - ftmp_2(x)) = -(1 - gamma(mu))*psi( x )
    ftmp_2[rb[1]] = Gamma(1<<mu) * psi;
    ftmp_2[rb[1]] -= psi;

    // utmp_1 = - Trace_spin [ ( 1 - gamma(mu) )*psi_{x+mu)*sigma^{dagger} ]
    //        = - Trace_spin [ sigma^{dagger} ( 1 - gamma_mu ) psi_{x+mu} ]

    utmp_1[rb[0]] = - trace(adj(sigma)*shift(ftmp_2, FORWARD, mu));
 
    // ftmp_2 = phi + ftmp_2 = (1 + gamma(mu))*phi( x) 
    ftmp_2[rb[1]] = Gamma(1<<mu) * phi;
    ftmp_2[rb[1]] += phi;

    // utmp_1 += ( 1 + gamma(mu) )*phi_{x+mu)*rho^{dagger}_x 
    utmp_1[rb[0]] += trace(adj(rho)*shift(ftmp_2, FORWARD, mu));

    // dsdu[mu][0] += u[mu][0] * utmp_1 
    //                = u[mu][0] [   ( 1 - gamma(mu) )*psi_{x+mu)*sigma^{dagger}_x
    //                             + ( 1 + gamma(mu) )*phi_{x+mu)*rho^{dagger}_x   ]
    ds_u[mu][rb[0]] += prefactor * u[mu] * utmp_1;
      
    // Checkerboard 1

    // ftmp_2 = -(rho - ftmp_2) = -(1 - gamma(mu))*rho( x ) 
    ftmp_2[rb[0]] = Gamma(1<<mu)*rho;
    ftmp_2[rb[0]] -= rho;

    // utmp_1 = ( 1 - gamma(mu) )*rho_{x+mu)*phi^{dagger}_x
    utmp_1[rb[1]] = -trace(adj(phi)*shift(ftmp_2, FORWARD, mu));
      
    // ftmp_2 = (gamma(mu))*sigma 
    ftmp_2[rb[0]] = Gamma(1<<mu)*sigma;
    ftmp_2[rb[0]] += sigma;


    utmp_1[rb[1]] += trace(adj(psi)*shift(ftmp_2, FORWARD, mu));
    ds_u[mu][rb[1]] += prefactor * u[mu] * utmp_1;

  }

  END_CODE();
}
