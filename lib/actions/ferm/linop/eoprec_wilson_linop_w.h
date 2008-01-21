// -*- C++ -*-
// $Id: eoprec_wilson_linop_w.h,v 3.2 2008-01-21 20:18:50 edwards Exp $
/*! \file
 *  \brief Even-odd preconditioned Wilson fermion linear operator
 */

#ifndef __eoprec_wilson_linop_w_h__
#define __eoprec_wilson_linop_w_h__

#include "eoprec_constdet_linop.h"
#include "actions/ferm/linop/dslash_w.h"
#include "io/aniso_io.h"


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

  class EvenOddPrecWilsonLinOp : public EvenOddPrecConstDetLinearOperator<LatticeFermion, 
				 multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Partial constructor
    EvenOddPrecWilsonLinOp() {}

    //! Full constructor
    EvenOddPrecWilsonLinOp(Handle< FermState<T,P,Q> > fs, 
			   const Real& Mass_)
    {create(fs,Mass_);}

    //! Full constructor with Anisotropy
    EvenOddPrecWilsonLinOp(Handle< FermState<T,P,Q> > fs,
			   const Real& Mass_,
			   const AnisoParam_t& aniso)
    {create(fs,Mass_,aniso);}

    //! Full constructor with array of coefficients
    EvenOddPrecWilsonLinOp(Handle< FermState<T,P,Q> > fs,
			   const Real& Mass_,
			   const multi1d<Real>& coeffs_)
    {create(fs,Mass_,coeffs_);}

    //! Destructor is automatic
    ~EvenOddPrecWilsonLinOp() {}

    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const {return D.getFermBC();}

    //! Creation routine
    void create(Handle< FermState<T,P,Q> > fs, 
		const Real& Mass_);

    //! Creation routine with Anisotropy
    void create(Handle< FermState<T,P,Q> > fs,
		const Real& Mass_,
		const AnisoParam_t& aniso);

    //! Creation routine with general coefficients
    void create(Handle< FermState<T,P,Q> > fs,
		const Real& Mass_,
		const multi1d<Real>& coeffs_);

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
      for(int mu=0; mu < Nd; mu++) { 
	ds_u[mu] = zero;
      }
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
      for(int mu=0; mu < Nd; mu++) { 
	ds_u[mu] = zero;
      }
    }


    //! Return flops performed by the operator()
    unsigned long nFlops() const;

  private:
    Real          Mass;    /*!< yep, the mass */
    multi1d<Real> coeffs;  /*!< Nd array of coefficients of terms in the action */

    Real          fact;    /*<! tmp holding  Nd+Mass */
    Real          invfact; /*!< tmp holding  1/(Nd+Mass) */

    WilsonDslash  D;       /*!< Wilson dslash term */
  };

} // End Namespace Chroma


#endif
