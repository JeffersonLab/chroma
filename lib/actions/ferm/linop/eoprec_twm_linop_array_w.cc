// $Id: eoprec_twm_linop_array_w.cc,v 1.1 2008-11-04 18:42:58 edwards Exp $
/*! \file
 *  \brief Even-odd preconditioned Twisted-mass linop where each flavor is one of two array elements
 */

#include "actions/ferm/linop/eoprec_parwilson_linop_w.h"

namespace Chroma 
{ 
  //! Creation routine
  /*!
   * \param fs 	       gauge field     	       (Read)
   * \param Mass_      fermion mass   	       (Read)
   * \param mu_sigma_  first parity mass term  (Read)
   */
  void EvenOddPrecParWilsonLinOp::create(Handle< FermState<T,P,Q> > fs,
					 const Real& Mass_, 
					 const Real& mu_sigma_, const Real& mu_delta_)
  {
    // This needs work
    Mass = Mass_;
    H = H_;
//    u = u_;
    D.create(fs);

    // This needs work
    fact = Nd + Mass;
    Real tmp = 1.0 / (fact*fact + mu_sigma_);
    invfact1 = fact * tmp;
    invfact2 = H * tmp;
  }


  //! Apply the the even-even block onto a source vector
  void 
  EvenOddPrecParWilsonLinOp::evenEvenLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
					   enum PlusMinus isign) const
  {
    switch (isign)
    {
    case PLUS:
      chi[rb[0]] = fact*psi + GammaConst<Ns,Ns*Ns-1>()*(H*timesI(psi));
      break;

    case MINUS:
      chi[rb[0]] = fact*psi - GammaConst<Ns,Ns*Ns-1>()*(H*timesI(psi));
      break;
    }
  }


  //! Return flops performed by the operator()
  unsigned long EvenOddPrecParWilsonLinOp::nFlops() const
  { 
    unsigned long cbsite_flops = 2*D.nFlops()+16*Nc*Ns;
    return cbsite_flops*(Layout::sitesOnNode()/2);
  }


  //! Apply the inverse of the even-even block onto a source vector
  void 
  EvenOddPrecParWilsonLinOp::evenEvenInvLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
					      enum PlusMinus isign) const
  {
    //  tmp   =  D'   (1-i isign H gamma_5)   D'    Psi
    //     O      O,E                          E,O     O
    switch (isign)
    {
    case PLUS:
      chi[rb[0]] = invfact1*psi - GammaConst<Ns,Ns*Ns-1>()*(invfact2*timesI(psi));
      break;

    case MINUS:
      chi[rb[0]] = invfact1*psi + GammaConst<Ns,Ns*Ns-1>()*(invfact2*timesI(psi));
      break;
    }
  }
  

  //! Apply the the odd-odd block onto a source vector
  void 
  EvenOddPrecParWilsonLinOp::oddOddLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
					 enum PlusMinus isign) const
  {
    switch (isign)
    {
    case PLUS:
      chi[rb[1]] = fact*psi + GammaConst<Ns,Ns*Ns-1>()*(H*timesI(psi));
      break;

    case MINUS:
      chi[rb[1]] = fact*psi - GammaConst<Ns,Ns*Ns-1>()*(H*timesI(psi));
      break;
    }
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
  EvenOddPrecParWilsonLinOp::evenOddLinOp(LatticeFermion& chi, 
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
  EvenOddPrecParWilsonLinOp::oddEvenLinOp(LatticeFermion& chi, 
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
  void 
  EvenOddPrecParWilsonLinOp::operator()(LatticeFermion & chi, 
					const LatticeFermion& psi, 
					enum PlusMinus isign) const
  {
    LatticeFermion tmp1, tmp2, tmp3;  // if an array is used here, 

    moveToFastMemoryHint(tmp1);
    moveToFastMemoryHint(tmp2);
    moveToFastMemoryHint(tmp3);

    Real mquarter = -0.25;

    //  tmp   =  D'   (1-i isign H gamma_5)   D'    Psi
    //     O      O,E                          E,O     O
    switch (isign)
    {
    case PLUS:
      // tmp1[0] = D_eo psi[1]
      D.apply(tmp1, psi, isign, 0);

      tmp2[rb[0]] = invfact1*tmp1 - invfact2*(GammaConst<Ns,Ns*Ns-1>()*timesI(tmp1));

      // tmp2[1] = D_oe tmp2[0]
      D.apply(tmp3, tmp2, isign, 1);

      chi[rb[1]] = fact*psi + mquarter*tmp3;
      chi[rb[1]] += H*(GammaConst<Ns,Ns*Ns-1>()*timesI(psi));
      break;

    case MINUS:
      // tmp1[0] = D_eo psi[1]
      D.apply(tmp1, psi, isign, 0);

      tmp2[rb[0]] = invfact1*tmp1 + invfact2*(GammaConst<Ns,Ns*Ns-1>()*timesI(tmp1));

      // tmp2[1] = D_oe tmp2[0]
      D.apply(tmp3, tmp2, isign, 1);

      chi[rb[1]] = fact*psi + mquarter*tmp3;
      chi[rb[1]] -= H*(GammaConst<Ns,Ns*Ns-1>()*timesI(psi));
      break;
    }

    getFermBC().modifyF(chi, rb[1]);
  }


  //! Derivative of even-odd linop component
  void 
  EvenOddPrecParWilsonLinOp::derivEvenOddLinOp(multi1d<LatticeColorMatrix>& ds_u,
					       const LatticeFermion& chi, 
					       const LatticeFermion& psi, 
					       enum PlusMinus isign) const
  {
    START_CODE();

    ds_u.resize(Nd);

    D.deriv(ds_u, chi, psi, isign, 0);
    for(int mu=0; mu < Nd; mu++) {
      ds_u[mu] *=  Real(-0.5);
    }

    END_CODE();
  }


  //! Derivative of odd-even linop component
  void 
  EvenOddPrecParWilsonLinOp::derivOddEvenLinOp(multi1d<LatticeColorMatrix>& ds_u,
					       const LatticeFermion& chi, 
					       const LatticeFermion& psi, 
					       enum PlusMinus isign) const
  {
    START_CODE();

    ds_u.resize(Nd);

    D.deriv(ds_u, chi, psi, isign, 1);
    for(int mu=0; mu < Nd; mu++) { 
      ds_u[mu]  *= Real(-0.5);
    }
    END_CODE();
  }



#if 0
// Code is here only as a reference. It should be deleted at some point.
// This is converted szin code.
// The above derivs will work, calling evenEvenInvLinOp to get that
// contribution. 
// However, for more performance, one should copy the regular Wilson
// deriv() routine and modify it.

#error "ANCIENT SZIN CODE"
#error "Not quite correct implementation"


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
  EvenOddPrecParWilsonFermAct::deriv(multi1d<LatticeColorMatrix>& ds_u,
				     const LatticeFermion& chi, 
				     const LatticeFermion& psi,
				     enum PlusMinus isign) const
  {
    START_CODE();
  
    QDPIO::cerr << "EvenOddPrecParWilsonFermAct::dsdu not implemented" << endl;
    QDP_abort(1);

    const multi1d<LatticeColorMatrix>& u = state->getLinks();
				 
    LatticeColorMatrix utmp_1;      moveToFastMemoryHint(utmp_1);
    LatticeFermion phi;             moveToFastMemoryHint(phi);
    LatticeFermion rho;             moveToFastMemoryHint(rho);
    LatticeFermion sigma;           moveToFastMemoryHint(sigma);
    LatticeFermion ftmp_1;          moveToFastMemoryHint(ftmp_1);
    LatticeFermion ftmp_2;          moveToFastMemoryHint(ftmp_2);
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

      /* THIS NEEDS TO BE CHANGED (removed) */
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

      /* THIS NEEDS TO BE CHANGED (removed) */
      ds_u[mu][cb] += u[mu][cb] * utmp_1;
      
    }
        
    getFermBC().zero(ds_u);

    END_CODE();
  }
#endif

} // End Namespace Chroma

