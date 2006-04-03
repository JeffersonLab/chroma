// -*- C++ -*-
// $Id: unprec_wilson_linop_w.h,v 3.0 2006-04-03 04:58:52 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Wilson fermion linear operator
 */

#ifndef __unprec_wilson_linop_w_h__
#define __unprec_wilson_linop_w_h__

#include "linearop.h"
#include "actions/ferm/linop/dslash_w.h"


namespace Chroma 
{ 
  //! Unpreconditioned Wilson-Dirac operator
  /*!
   * \ingroup linop
   *
   * This routine is specific to Wilson fermions!
   *
   *                                                      ~      ~+
   * This subroutine applies the unpreconditioned matrix  M  or  M   the vector
   * Psi,
   *
   *      	       	   {   ~
   *      	       	   {   M(U) . Psi      	       if  ISign = PLUS
   *      	   Chi  =  {
   *      	       	   {   ~   +
   *      	       	   {   M(U)  . Psi     	       if  ISign = MINUS
   *
   * Algorithm:

   * The kernel for Wilson fermions is

   *      M  =  (d+M) - (1/2) D'
   *
   */
  
  class UnprecWilsonLinOp : public UnprecLinearOperator<LatticeFermion, 
                    multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Partial constructor
    UnprecWilsonLinOp() {}

    //! Full constructor
    UnprecWilsonLinOp(Handle< FermState<T,P,Q> > fs, const Real& Mass_)
      {create(fs,Mass_);}

    //! Full constructor with Anisotropy
    UnprecWilsonLinOp(Handle< FermState<T,P,Q> > fs,
		      const Real& Mass_,
		      const AnisoParam_t& aniso)
      {create(fs,Mass_,aniso);}

    //! Destructor is automatic
    ~UnprecWilsonLinOp() {}

    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const {return D.getFermBC();}

    //! Creation routine
    void create(Handle< FermState<T,P,Q> > fs,
		const Real& Mass_);

    //! Creation routine with Anisotropy
    void create(Handle< FermState<T,P,Q> > fs,
		const Real& Mass_,
		const AnisoParam_t& aniso);

    //! Apply the operator onto a source vector
    void operator() (LatticeFermion& chi, const LatticeFermion& psi, enum PlusMinus isign) const;

    //! Derivative of unpreconditioned Wilson dM/dU
    void deriv(multi1d<LatticeColorMatrix>& ds_u, 
	       const LatticeFermion& chi, const LatticeFermion& psi, 
	       enum PlusMinus isign) const;

    //! Return flops performed by the operator()
    unsigned long nFlops() const;

  private:
    Real fact;  // tmp holding  Nd+Mass

    Real Mass;
    WilsonDslash D;
  };

} // End Namespace Chroma


#endif
