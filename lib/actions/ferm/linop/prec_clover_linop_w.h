// -*- C++ -*-
// $Id: prec_clover_linop_w.h,v 1.3 2005-02-23 19:29:20 edwards Exp $
/*! \file
 *  \brief Even-odd preconditioned Clover fermion linear operator
 */

#error "THIS ROUTINE IS NOT CORRECTLY IMPLEMENTED"


#ifndef __prec_clover_linop_w_h__
#define __prec_clover_linop_w_h__

#include "linearop.h"
#include "actions/ferm/linop/dslash_w.h"
#include "io/aniso_io.h"


namespace Chroma 
{ 
  //! Even-odd preconditioned Clover-Dirac operator
  /*!
   * \ingroup linop
   *
   * This routine is specific to Wilson fermions!
   *
   * The kernel for Clover fermions is
   *
   *      M  =  (d+M) - (1/2) D'
   */
  class EvenOddPrecCloverLinOp : public EvenOddPrecLinearOperator< LatticeFermion, multi1d<LatticeColorMatrix> >
  {
  public:
    //! Partial constructor
    EvenOddPrecCloverLinOp() {}

    //! Full constructor
    EvenOddPrecCloverLinOp(const multi1d<LatticeColorMatrix>& u_, const Real& Mass_)
    {create(u_,Mass_);}

    //! Full constructor with Anisotropy
    EvenOddPrecCloverLinOp(const multi1d<LatticeColorMatrix>& u_, 
			   const Real& Mass_,
			   const AnisoParam_t& aniso)
    {create(u_,Mass_,aniso);}

    //! Destructor is automatic
    ~EvenOddPrecCloverLinOp() {}

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

    // Override inherited one with a few more funkies
    void operator()(LatticeFermion& chi, const LatticeFermion& psi, 
		    enum PlusMinus isign) const;


    //! Apply the even-even block onto a source vector
    void derivEvenEvenLinOp(multi1d<LatticeColorMatrix>& ds_u, 
			    const LatticeFermion& chi, const LatticeFermion& psi, 
			    enum PlusMinus isign) const
    {
      QDPIO::cerr << "Clover: not implemented" << endl;
      QDP_abort(1);
    }
  
    //! Apply the the even-odd block onto a source vector
    void derivEvenOddLinOp(multi1d<LatticeColorMatrix>& ds_u, 
			   const LatticeFermion& chi, const LatticeFermion& psi, 
			   enum PlusMinus isign) const
    {
      QDPIO::cerr << "Clover: not implemented" << endl;
      QDP_abort(1);
    }
 
    //! Apply the the odd-even block onto a source vector
    void derivOddEvenLinOp(multi1d<LatticeColorMatrix>& ds_u, 
			   const LatticeFermion& chi, const LatticeFermion& psi, 
			   enum PlusMinus isign) const
    {
      QDPIO::cerr << "Clover: not implemented" << endl;
      QDP_abort(1);
    }

    //! Apply the the odd-odd block onto a source vector
    void derivOddOddLinOp(multi1d<LatticeColorMatrix>& ds_u, 
			  const LatticeFermion& chi, const LatticeFermion& psi, 
			  enum PlusMinus isign) const
    {
      QDPIO::cerr << "Clover: not implemented" << endl;
      QDP_abort(1);
    }

  private:
    Real fact;  // tmp holding  Nd+Mass
    Real invfact;  // tmp holding  1/(Nd+Mass)

    Real Mass;
    WilsonDslash D;
  };

}; // End Namespace Chroma


#endif
