// -*- C++ -*-
// $Id: prec_parwilson_linop_w.h,v 1.2 2004-12-12 21:22:16 edwards Exp $
/*! \file
 *  \brief Even-odd preconditioned Wilson fermion linear operator with parity breaking term
 */

#ifndef __prec_parwilson_linop_w_h__
#define __prec_parwilson_linop_w_h__

#include "linearop.h"
#include "actions/ferm/linop/dslash_w.h"

using namespace QDP;

namespace Chroma 
{ 
  //! Even-odd preconditioned Wilson fermion linear operator with parity breaking term
  /*!
   * \ingroup linop
   *
   * This routine is specific to Wilson fermions!
   *
   * The kernel for Wilson fermions with a parity breaking term is
   *
   *      M  =  (d+M) + i*H*gamma_5  - (1/2) D'
   */

  class EvenOddPrecParWilsonLinOp : public EvenOddPrecLinearOperator< LatticeFermion, multi1d<LatticeColorMatrix> >
  {
  public:
    //! Partial constructor
    EvenOddPrecParWilsonLinOp() {}

    //! Full constructor
    EvenOddPrecParWilsonLinOp(const multi1d<LatticeColorMatrix>& u_, 
			      const Real& Mass_, const Real& H_)
    {create(u_,Mass_,H_);}

    //! Destructor is automatic
    ~EvenOddPrecParWilsonLinOp() {}

    //! Creation routine
    void create(const multi1d<LatticeColorMatrix>& u_, 
		const Real& Mass_, const Real& H_);

    //! Apply the the even-even block onto a source vector
    inline
    void evenEvenLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
		       enum PlusMinus isign) const
    {
      chi[rb[0]] = fact*psi + Gamma(Ns*Ns-1)*(H*timesI(psi));
    }

    //! Apply the inverse of the even-even block onto a source vector
    inline 
    void evenEvenInvLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
			  enum PlusMinus isign) const
    {
      //  tmp   =  D'   (1-i isign H gamma_5)   D'    Psi
      //     O      O,E                          E,O     O
      chi[rb[0]] = invfact1*psi - Gamma(Ns*Ns-1)*(invfact2*timesI(psi));
    }
  
    //! Apply the the even-odd block onto a source vector
    void evenOddLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
		      enum PlusMinus isign) const;

    //! Apply the the odd-even block onto a source vector
    void oddEvenLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
		      enum PlusMinus isign) const;

    //! Apply the the odd-odd block onto a source vector
    inline 
    void oddOddLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
		     enum PlusMinus isign) const
    {
      chi[rb[0]] = fact*psi + Gamma(Ns*Ns-1)*(H*timesI(psi));
    }

  private:
    Real fact;  // tmp holding  Nd+Mass
    Real invfact1;  // tmp holding  [1/(Nd+Mass)]*[1/(1+H^2/(Nd+Mass)^2)]
    Real invfact2;  // tmp holding  H/[(Nd+Mass)^2 + H^2]

    Real Mass;
    Real H;
    multi1d<LatticeColorMatrix> u;
    WilsonDslash D;
  };

}; // End Namespace Chroma

using namespace Chroma;

#endif
