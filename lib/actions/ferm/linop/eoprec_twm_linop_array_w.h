// -*- C++ -*-
// $Id: eoprec_twm_linop_array_w.h,v 1.1 2008-11-04 18:42:58 edwards Exp $
/*! \file
 *  \brief Even-odd preconditioned Twisted-mass linop where each flavor is one of two array elements
 */

#ifndef __eoprec_twm_linop_array_w_h__
#define __eoprec_twm_linop_array_w_h__

#include "eoprec_constdet_linop.h"
#include "actions/ferm/linop/dslash_w.h"

namespace Chroma 
{ 
  //! Even-odd preconditioned Twisted-mass linop where each flavor is one of two array elements
  /*!
   * \ingroup linop
   *
   * This routine is specific to Wilson fermions!
   *
   * The kernel for Wilson fermions with a parity breaking term is
   *
   *      M  =  (d+M) + i*H*gamma_5  - (1/2) D'
   */

  class EvenOddPrecTwmLinOpArray : public EvenOddPrecConstDetLinearOperator<LatticeFermion, 
				    multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Partial constructor
    EvenOddPrecTwmLinOpArray() {}

    //! Full constructor
    EvenOddPrecTwmLinOpArray(Handle< FermState<T,P,Q> > fs,
			     const Real& Mass_, const Real& mu_sigma_, const Real& mu_delta_)
    {create(fs,Mass_,mu_sigma_,mu_delta_);}

    //! Destructor is automatic
    ~EvenOddPrecTwmLinOpArray() {}

    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const {return D.getFermBC();}

    //! Creation routine
    void create(Handle< FermState<T,P,Q> > fs,
		const Real& Mass_, const Real& mu_delta_, const Real& mu_sigma_);

    //! Apply the the even-even block onto a source vector
    void evenEvenLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
		       enum PlusMinus isign) const;

    //! Apply the inverse of the even-even block onto a source vector
    void evenEvenInvLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
			  enum PlusMinus isign) const;
  
    //! Apply the the even-odd block onto a source vector
    void evenOddLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
		      enum PlusMinus isign) const;

    //! Apply the the odd-even block onto a source vector
    void oddEvenLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
		      enum PlusMinus isign) const;

    //! Apply the the odd-odd block onto a source vector
    void oddOddLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
		     enum PlusMinus isign) const;

    //! Override inherited one with a few more funkies
    void operator()(LatticeFermion& chi, const LatticeFermion& psi, 
		    enum PlusMinus isign) const;

    //! Apply the even-even block onto a source vector
    inline
    void derivEvenEvenLinOp(multi1d<LatticeColorMatrix>& ds_u, 
			    const LatticeFermion& chi, const LatticeFermion& psi, 
			    enum PlusMinus isign) const
    {
      ds_u.resize(Nd);
      ds_u = zero;
    }

    //! Apply the the even-odd block onto a source vector
    void derivEvenOddLinOp(multi1d<LatticeColorMatrix>& ds_u, 
			   const LatticeFermion& chi, const LatticeFermion& psi, 
			   enum PlusMinus isign) const;
 
    //! Apply the the odd-even block onto a source vector
    void derivOddEvenLinOp(multi1d<LatticeColorMatrix>& ds_u, 
			   const LatticeFermion& chi, const LatticeFermion& psi, 
			   enum PlusMinus isign) const;

    //! Apply the the odd-odd block onto a source vector
    inline 
    void derivOddOddLinOp(multi1d<LatticeColorMatrix>& ds_u, 
			  const LatticeFermion& chi, const LatticeFermion& psi, 
			  enum PlusMinus isign) const
    {
      ds_u.resize(Nd);
      ds_u = zero;
    }

    //! Return flops performed by the operator()
    unsigned long nFlops() const;

  private:
    Real fact;  // tmp holding  Nd+Mass
    Real invfact1;  // tmp holding  [1/(Nd+Mass)]*[1/(1+H^2/(Nd+Mass)^2)]
    Real invfact2;  // tmp holding  H/[(Nd+Mass)^2 + H^2]

    Real Mass;
    Real mu_sigma;
    Real mu_delta;
    WilsonDslash D;
  };

} // End Namespace Chroma


#endif
