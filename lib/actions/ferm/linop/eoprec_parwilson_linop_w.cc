// $Id: eoprec_parwilson_linop_w.cc,v 3.1 2006-10-19 16:01:30 edwards Exp $
/*! \file
 *  \brief Even-odd preconditioned Wilson fermion linear operator with parity breaking term
 */

#include "actions/ferm/linop/eoprec_parwilson_linop_w.h"

namespace Chroma 
{ 
  //! Creation routine
  /*!
   * \param u_ 	   gauge field     	       (Read)
   * \param Mass_    fermion mass   	       (Read)
   * \param H_       parity breaking term	       (Read)
   */
  void EvenOddPrecParWilsonLinOp::create(Handle< FermState<T,P,Q> > fs,
					 const Real& Mass_, const Real& H_)
  {
    Mass = Mass_;
    H = H_;
//    u = u_;
    D.create(fs);

    fact = Nd + Mass;
    Real tmp = 1.0 / (fact*fact + H*H);
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

 void 
 EvenOddPrecParWilsonLinOp::evenEvenLinOp(Tower<T>& chi, const Tower<T>& psi,
  				         const P& p,
				         enum PlusMinus isign) const
  {
   switch (isign)
    {
    case PLUS:
    {
      for(int i=0; i < psi.size(); i++) {  
          chi[i][rb[0]] = fact*psi[i] + GammaConst<Ns,Ns*Ns-1>()*(H*timesI(psi[i]));
      }
    }
    break;

    case MINUS:
    {
      for(int i=0; i < psi.size(); i++)  { 
       chi[i][rb[0]] = fact*psi[i] - GammaConst<Ns,Ns*Ns-1>()*(H*timesI(psi[i]));
      }
    }
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
 
  void  
  EvenOddPrecParWilsonLinOp::evenEvenInvLinOp(Tower<T>& chi, const Tower<T>& psi,
				              const P& p, enum PlusMinus isign) const
  {
   //  tmp   =  D'   (1-i isign H gamma_5)   D'    Psi
    //     O      O,E                          E,O     O
    switch (isign)
    {
    case PLUS:
    { 
      for(int i=0; i < psi.size(); i++) {  
        chi[i][rb[0]] = invfact1*psi[i] - GammaConst<Ns,Ns*Ns-1>()*(invfact2*timesI(psi[i]));
      }
    }  
    break;

    case MINUS:
    {
      for(int i=0; i < psi.size(); i++) { 
         chi[i][rb[0]] = invfact1*psi[i] + GammaConst<Ns,Ns*Ns-1>()*(invfact2*timesI(psi[i]));
      }
    }  
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

  void 
  EvenOddPrecParWilsonLinOp::oddOddLinOp(Tower<T>& chi, const Tower<T>& psi,
				     const P& p,
				     enum PlusMinus isign) const
  {
 switch (isign)
    {
    case PLUS:
    {
      for(int i=0; i < psi.size(); i++) {  
          chi[i][rb[1]] = fact*psi[i] + GammaConst<Ns,Ns*Ns-1>()*(H*timesI(psi[i]));
      }
    }
    break;

    case MINUS:
    {
      for(int i=0; i < psi.size(); i++)  { 
       chi[i][rb[1]] = fact*psi[i] - GammaConst<Ns,Ns*Ns-1>()*(H*timesI(psi[i]));
      }
    }
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

 void 
  EvenOddPrecParWilsonLinOp::evenOddLinOp(Tower<T>& chi, const Tower<T>& psi,
				       const P& p,
				       enum PlusMinus isign) const
  { 
     START_CODE();
     Real mhalf = -0.5;
       
    D.applyTower(chi, psi, p, isign, 0);
    for(int i=0; i < chi.size(); i++) {
      chi[i][rb[0]] *= mhalf;
    }

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

 void 
  EvenOddPrecParWilsonLinOp::oddEvenLinOp(Tower<T>& chi, const Tower<T>& psi,
				       const P& p,
				       enum PlusMinus isign) const
  { 
     START_CODE();
     Real mhalf = -0.5;
       
    D.applyTower(chi, psi, p, isign, 1);
    for(int i=0; i < chi.size(); i++) {
      chi[i][rb[1]] *= mhalf;
    }

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

  void 
  EvenOddPrecParWilsonLinOp::operator()(Tower<T>& chi, const Tower<T>& psi,
				       const P& p,
				       enum PlusMinus isign) const
  {
    int N = psi.size();
    Tower<T> tmp1(N);
    Tower<T> tmp2(N);
    Tower<T> tmp3(N);
  
    Real mquarter = -0.25;

    //  tmp   =  D'   (1-i isign H gamma_5)   D'    Psi
    //     O      O,E                          E,O     O
    switch (isign)
    {
    case PLUS:
    {
      // tmp1[0] = D_eo psi[1]
      D.applyTower(tmp1, psi, p, isign, 0);

      for(int i=0; i < N; i++) { 
        tmp2[i][rb[0]] = invfact1*tmp1[i] - invfact2*(GammaConst<Ns,Ns*Ns-1>()*timesI(tmp1[i]));
      }

      // tmp2[1] = D_oe tmp2[0]
      D.applyTower(tmp3, tmp2, p, isign, 1);
      
      for(int i=0; i < N; i++) { 
      chi[i][rb[1]] = fact*psi[i] + mquarter*tmp3[i];
      chi[i][rb[1]] += H*(GammaConst<Ns,Ns*Ns-1>()*timesI(psi[i]));
      }
    }
    break;

    case MINUS:
    {
      // tmp1[0] = D_eo psi[1]
      D.applyTower(tmp1, psi, p, isign, 0);
   
      for(int i=0; i < N; i++) { 
        tmp2[i][rb[0]] = invfact1*tmp1[i] + invfact2*(GammaConst<Ns,Ns*Ns-1>()*timesI(tmp1[i]));
      }

      // tmp2[1] = D_oe tmp2[0]
      D.applyTower(tmp3, tmp2, p, isign, 1);

      for(int i=0; i < N; i++) { 
        chi[i][rb[1]] = fact*psi[i] + mquarter*tmp3[i];
        chi[i][rb[1]] -= H*(GammaConst<Ns,Ns*Ns-1>()*timesI(psi[i]));
      }
     }
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

  void 
  EvenOddPrecParWilsonLinOp::derivEvenOddLinOp(TowerArray< PQTraits<Q>::Base_t>& ds_u,
					       const Tower<T>& chi,
					       const Tower<T>& psi,
					       const P& p,
					       enum PlusMinus isign)
  {
    START_CODE();

    D.deriv(ds_u, chi, psi, p, isign, 0);
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

 void 
  EvenOddPrecParWilsonLinOp::derivOddEvenLinOp(TowerArray< PQTraits<Q>::Base_t>& ds_u,
					       const Tower<T>& chi,
					       const Tower<T>& psi,
					       const P& p,
					       enum PlusMinus isign)
  {
    START_CODE();

    D.deriv(ds_u, chi, psi, p, isign, 1);
    for(int mu=0; mu < Nd; mu++) {
      ds_u[mu] *=  Real(-0.5);
    }

    END_CODE();
 
  }


} // End Namespace Chroma

