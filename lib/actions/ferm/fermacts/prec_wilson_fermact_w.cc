// $Id: prec_wilson_fermact_w.cc,v 1.12 2004-09-08 02:48:25 edwards Exp $
/*! \file
 *  \brief Even-odd preconditioned Wilson fermion action
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/prec_wilson_fermact_w.h"
#include "actions/ferm/linop/prec_wilson_linop_w.h"
#include "actions/ferm/linop/lmdagm.h"

#include "actions/ferm/fermacts/fermfactory_w.h"

namespace Chroma
{
  //! Hooks to register the class with the fermact factory
  namespace EvenOddPrecWilsonFermActEnv
  {
    //! Callback function
    WilsonTypeFermAct<LatticeFermion>* createFermAct(Handle< FermBC<LatticeFermion> > fbc,
						     XMLReader& xml_in,
						     const std::string& path)
    {
      return new EvenOddPrecWilsonFermAct(fbc, EvenOddPrecWilsonFermActParams(xml_in, path));
    }

     //! Name to be used
    const std::string name = "WILSON";

    //! Register the Wilson fermact
    const bool registered = Chroma::TheWilsonTypeFermActFactory::Instance().registerObject(name, createFermAct); 
  }


  //! Read parameters
  EvenOddPrecWilsonFermActParams::EvenOddPrecWilsonFermActParams(XMLReader& xml, const string& path)
  {
    XMLReader paramtop(xml, path);

    // Read the stuff for the action
    if (paramtop.count("Mass") != 0) 
    {
      read(paramtop, "Mass", Mass);
      if (paramtop.count("Kappa") != 0) 
      {
	QDPIO::cerr << "Error: found both a Kappa and a Mass tag" << endl;
	QDP_abort(1);
      }
    }
    else if (paramtop.count("Kappa") != 0)
    {
      Real Kappa;
      read(paramtop, "Kappa", Kappa);
      Mass = kappaToMass(Kappa);    // Convert Kappa to Mass
    }
    else
    {
      QDPIO::cerr << "Error: neither Mass or Kappa found" << endl;
      QDP_abort(1);
    }

    // There is always an aniso Param for wilson, so set it to default
    initHeader(anisoParam);

    //  Read optional anisoParam.
    if (paramtop.count("AnisoParam") != 0) 
      read(paramtop, "AnisoParam", anisoParam);
  }

  //! Read parameters
  void read(XMLReader& xml, const string& path, EvenOddPrecWilsonFermActParams& param)
  {
    EvenOddPrecWilsonFermActParams tmp(xml, path);
    param = tmp;
  }



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

    if (param.anisoParam.anisoP)
      foo = new EvenOddPrecWilsonLinOp(state->getLinks(),param.Mass,param.anisoParam);
    else
      foo = new EvenOddPrecWilsonLinOp(state->getLinks(),param.Mass);

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
  
    if (param.anisoParam.anisoP)
    {
      QDPIO::cerr << "Currently do not support anisotropy" << endl;
      QDP_abort(1);
    }

    Real prefactor = -Real(1)/(4*(Real(Nd) + param.Mass));
				 
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
    D.apply(rho, psi, PLUS, 0);

    /* sigma = Dslash_dag(0 <- 1) * phi */
    D.apply(sigma, phi, MINUS, 0);

    for(int mu = 0; mu < Nd; ++mu)
    {

      // ftmp_2(x) = -(psi(x) - ftmp_2(x)) = -(1 - gamma(mu))*psi( x )
      ftmp_2[rb[1]] = Gamma(1<<mu) * psi;
      ftmp_2[rb[1]]  -= psi;


      // utmp_1 = - Trace_spin [ ( 1 - gamma(mu) )*psi_{x+mu)*sigma^{dagger} ]
      //        = - Trace_spin [ sigma^{dagger} ( 1 - gamma_mu ) psi_{x+mu} ]
      utmp_1[rb[0]] = -traceSpin( outerProduct( shift(ftmp_2, FORWARD, mu), sigma) );

    
      // ftmp_2 = phi + ftmp_2 = (1 + gamma(mu))*phi( x) 
      ftmp_2[rb[1]] = Gamma(1<<mu) * phi;
      ftmp_2[rb[1]] += phi;

      // utmp_1 += ( 1 + gamma(mu) )*phi_{x+mu)*rho^{dagger}_x 
      utmp_1[rb[0]] += traceSpin( outerProduct( shift(ftmp_2, FORWARD, mu), rho) );

      // dsdu[mu][0] += u[mu][0] * utmp_1 
      //                = u[mu][0] [   ( 1 - gamma(mu) )*psi_{x+mu)*sigma^{dagger}_x
      //                             + ( 1 + gamma(mu) )*phi_{x+mu)*rho^{dagger}_x   ]
      ds_u[mu][rb[0]] += prefactor * u[mu] * utmp_1;
      
      // Checkerboard 1

      // ftmp_2 = -(rho - ftmp_2) = -(1 - gamma(mu))*rho( x ) 
      ftmp_2[rb[0]] = Gamma(1<<mu)*rho;
      ftmp_2[rb[0]] -= rho;

      // utmp_1 = ( 1 - gamma(mu) )*rho_{x+mu)*phi^{dagger}_x
      utmp_1[rb[1]] = -traceSpin( outerProduct( shift(ftmp_2, FORWARD, mu), phi) );
      
      // ftmp_2 = (gamma(mu))*sigma 
      ftmp_2[rb[0]] = Gamma(1<<mu)*sigma;
      ftmp_2[rb[0]] += sigma;


      utmp_1[rb[1]] += traceSpin( outerProduct( shift(ftmp_2, FORWARD, mu), psi) );
      ds_u[mu][rb[1]] += prefactor * u[mu] * utmp_1;

    }

    END_CODE();
  }

}
