// -*- C++ -*-
// $Id: prec_wilson_linop_w.h,v 1.7 2004-12-12 21:22:16 edwards Exp $
/*! \file
 *  \brief Even-odd preconditioned Wilson fermion linear operator
 */

#ifndef __prec_wilson_linop_w_h__
#define __prec_wilson_linop_w_h__

#include "linearop.h"
#include "actions/ferm/linop/dslash_w.h"
#include "io/param_io.h"       // to get AnisoParam_t

using namespace QDP;

namespace Chroma 
{ 
  //! Even-odd preconditioned Wilson-Dirac operator
  /*!
   * \ingroup linop
   *
   * This routine is specific to Wilson fermions!
   *
   * The kernel for Wilson fermions is
   *
   *      M  =  (d+M) - (1/2) D'
   */
  class EvenOddPrecWilsonLinOp : public EvenOddPrecLinearOperator< LatticeFermion, multi1d<LatticeColorMatrix> >
  {
  public:
    //! Partial constructor
    EvenOddPrecWilsonLinOp() {}

    //! Full constructor
    EvenOddPrecWilsonLinOp(const multi1d<LatticeColorMatrix>& u_, const Real& Mass_)
    {create(u_,Mass_);}

    //! Full constructor with Anisotropy
    EvenOddPrecWilsonLinOp(const multi1d<LatticeColorMatrix>& u_, 
			   const Real& Mass_,
			   const AnisoParam_t& aniso)
    {create(u_,Mass_,aniso);}

    //! Destructor is automatic
    ~EvenOddPrecWilsonLinOp() {}

    //! Creation routine
    void create(const multi1d<LatticeColorMatrix>& u_, const Real& Mass_);

    //! Creation routine with Anisotropy
    void create(const multi1d<LatticeColorMatrix>& u_, 
		const Real& Mass_,
		const AnisoParam_t& aniso);

    //! Apply the the even-even block onto a source vector
    inline
    void evenEvenLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
		       enum PlusMinus isign) const
    {
      chi[rb[0]] = fact*psi;
    }

    //! Apply the inverse of the even-even block onto a source vector
    inline 
    void evenEvenInvLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
			  enum PlusMinus isign) const
    {
      chi[rb[0]] = invfact*psi;
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
      chi[rb[1]] = fact*psi;
    }

    //! Override inherited one with a few more funkies
    void operator()(LatticeFermion& chi, const LatticeFermion& psi, 
		    enum PlusMinus isign) const;


    //! Apply the even-even block onto a source vector
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
    void derivOddOddLinOp(multi1d<LatticeColorMatrix>& ds_u, 
			  const LatticeFermion& chi, const LatticeFermion& psi, 
			  enum PlusMinus isign) const
    {
      ds_u.resize(Nd);
      ds_u = zero;
    }

  private:
    Real fact;  // tmp holding  Nd+Mass
    Real invfact;  // tmp holding  1/(Nd+Mass)

    Real Mass;
    WilsonDslash D;
  };

}; // End Namespace Chroma

using namespace Chroma;

#endif
