// $Id: prec_wilson_linop_w.cc,v 1.8 2004-12-12 21:22:16 edwards Exp $
/*! \file
 *  \brief Even-odd preconditioned Wilson linear operator
 */

#include "chromabase.h"
#include "actions/ferm/linop/prec_wilson_linop_w.h"

namespace Chroma 
{ 
  //! Creation routine
  /*!
   * \param u_ 	  gauge field     	       (Read)
   * \param Mass_   fermion kappa   	       (Read)
   */
  void EvenOddPrecWilsonLinOp::create(const multi1d<LatticeColorMatrix>& u_, 
				      const Real& Mass_)
  {
    Mass = Mass_;
    D.create(u_);

    fact = Nd + Mass;
    invfact = 1/fact;
  }


  //! Creation routine with Anisotropy
  /*!
   * \param u_ 	  gauge field     	       (Read)
   * \param Mass_   fermion kappa   	       (Read)
   * \param aniso   anisotropy struct   	       (Read)
   */
  void EvenOddPrecWilsonLinOp::create(const multi1d<LatticeColorMatrix>& u_, 
				      const Real& Mass_,
				      const AnisoParam_t& aniso)
  {
    Mass = Mass_;

    multi1d<LatticeColorMatrix> u = u_;
    Real ff = where(aniso.anisoP, aniso.nu / aniso.xi_0, Real(1));
  
    if (aniso.anisoP)
    {
      // Rescale the u fields by the anisotropy
      for(int mu=0; mu < u.size(); ++mu)
      {
	if (mu != aniso.t_dir)
	  u[mu] *= ff;
      }
    }
    D.create(u);

    fact = 1 + (Nd-1)*ff + Mass;
    invfact = 1/fact;
  }


  //! Apply even-odd linop component
  /*!
   * The operator acts on the entire even sublattice
   *
   * \param chi 	  Pseudofermion field     	       (Write)
   * \param psi 	  Pseudofermion field     	       (Read)
   * \param isign   Flag ( PLUS | MINUS )   	       (Read)
   */
  void 
  EvenOddPrecWilsonLinOp::evenOddLinOp(LatticeFermion& chi, 
				       const LatticeFermion& psi, 
				       enum PlusMinus isign) const
  {
    START_CODE();

    Real mhalf = -0.5;

    D.apply(chi, psi, isign, 0);
    chi[rb[0]] *= mhalf;
  
    END_CODE();
  }

  //! Apply odd-even linop component
  /*!
   * The operator acts on the entire odd sublattice
   *
   * \param chi 	  Pseudofermion field     	       (Write)
   * \param psi 	  Pseudofermion field     	       (Read)
   * \param isign   Flag ( PLUS | MINUS )   	       (Read)
   */
  void 
  EvenOddPrecWilsonLinOp::oddEvenLinOp(LatticeFermion& chi, 
				       const LatticeFermion& psi, 
				       enum PlusMinus isign) const
  {
    START_CODE();

    Real mhalf = -0.5;

    D.apply(chi, psi, isign, 1);
    chi[rb[1]] *= mhalf;
  
    END_CODE();
  }



  //! Override inherited one with a few more funkies
  void EvenOddPrecWilsonLinOp::operator()(LatticeFermion & chi, 
					  const LatticeFermion& psi, 
					  enum PlusMinus isign) const
  {
    LatticeFermion tmp1, tmp2, tmp3;  // if an array is used here, 
    // the space is not reserved
  

    Real mquarterinvfact = -0.25*invfact;

    // tmp1[0] = D_eo psi[1]
    D.apply(tmp1, psi, isign, 0);

    // tmp2[1] = D_oe tmp1[0]
    D.apply(tmp2, tmp1, isign, 1);

    // Now we have tmp2[1] = D_oe D_eo psi[1]

    // now scale tmp2[1] with (-1/4)/fact = (-1/4)*(1/(Nd + m))
    // with a vscale -- using tmp2 on both sides should be OK, but
    // just to be safe use tmp3
    tmp3[rb[1]] = mquarterinvfact*tmp2;

    // now tmp3[1] should be = (-1/4)*(1/(Nd + m) D_oe D_eo psi[1]

    // Now get chi[1] = fact*psi[1] + tmp3[1]
    // with a vaxpy3 
    // chi[1] = (Nd + m) - (1/4)*(1/(Nd + m)) D_oe D_eo psi[1]
    //
    // in this order, this last job could be replaced with a 
    // vaxpy3_norm if we wanted the || M psi ||^2 over the subset.
    chi[rb[1]] = fact*psi + tmp3;
  }


  //! Derivative of even-odd linop component
  void 
  EvenOddPrecWilsonLinOp::derivEvenOddLinOp(multi1d<LatticeColorMatrix>& ds_u,
					    const LatticeFermion& chi, 
					    const LatticeFermion& psi, 
					    enum PlusMinus isign) const
  {
    START_CODE();

    QDPIO::cerr << "Prec wilson deriv: not implemented yet" << endl;
    QDP_abort(1);
  
    ds_u.resize(Nd);

    LatticeFermion tmp;
    LatticeFermion gamma_mu_psi;
    
    for(int mu = 0; mu < Nd; ++mu)
    {
      // Get gamma_mu * psi
      gamma_mu_psi = Gamma(1 << mu)*psi;

      // Construct the right derivative term...
      if( isign == PLUS ) 
      {
	// Undaggered:
	// ( 1 - gamma_mu ) psi
	tmp = psi - gamma_mu_psi;
      }
      else 
      { 
	// Daggered:
	// ( 1 + gamma_mu) psi
	tmp = psi + gamma_mu_psi;
      }
    
      // Do the shift here for delta_y, x+mu and form derivative
      ds_u[mu] = traceSpin(outerProduct(-Real(0.5)*shift(tmp, FORWARD, mu),chi));
    }
  
    END_CODE();
  }


  //! Derivative of odd-even linop component
  /*!
   * The operator acts on the entire odd sublattice
   *
   * \param chi 	  Pseudofermion field     	       (Write)
   * \param psi 	  Pseudofermion field     	       (Read)
   * \param isign   Flag ( PLUS | MINUS )   	       (Read)
   */
  void 
  EvenOddPrecWilsonLinOp::derivOddEvenLinOp(multi1d<LatticeColorMatrix>& ds_u,
					    const LatticeFermion& chi, 
					    const LatticeFermion& psi, 
					    enum PlusMinus isign) const
  {
    START_CODE();

    ds_u.resize(Nd);

    QDPIO::cerr << "Prec wilson deriv: not implemented yet" << endl;
    QDP_abort(1);
  
    END_CODE();
  }



#if 0

#error "This should be split apart into the even, odd critters to test the generic deriv() function. However, for optimization it can be done here since the diag pieces are trivial"

  //! Derivative of even-odd preconditioned linear operator
  /*!
   * In this function I assume that ds_u may already have the gauge piece in there...
   */
  void
  EvenOddPrecWilsonLinOp::deriv(multi1d<LatticeColorMatrix>& ds_u,
				const LatticeFermion& chi, const LatticeFermion& psi, 
				enum PlusMinus isign) const
  {
    START_CODE();
  
    ds_u.resize(Nd);

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

      // ds_u[mu][0] += u[mu][0] * utmp_1 
      //              = u[mu][0] [   ( 1 - gamma(mu) )*psi_{x+mu)*sigma^{dagger}_x
      //                           + ( 1 + gamma(mu) )*phi_{x+mu)*rho^{dagger}_x   ]
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

#endif

}; // End Namespace Chroma
