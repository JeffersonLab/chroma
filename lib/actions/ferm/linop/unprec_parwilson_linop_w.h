// -*- C++ -*-
// $Id: unprec_parwilson_linop_w.h,v 3.1 2007-02-22 21:11:47 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned Wilson fermion linear operator with parity breaking term
 */

#ifndef __unprec_parwilson_linop_w_h__
#define __unprec_parwilson_linop_w_h__

#include "linearop.h"
#include "actions/ferm/linop/dslash_w.h"


namespace Chroma 
{ 
  //! Unpreconditioned Wilson-Dirac operator with parity breaking term
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

   * Algorithm:
   *
   *
   * The kernel for Wilson fermions with a parity breaking term is
   *
   *      M  =  (d+M) + i*H*gamma_5  - (1/2) D'
   */

  class UnprecParWilsonLinOp : public UnprecLinearOperator<LatticeFermion, 
            multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Partial constructor
    UnprecParWilsonLinOp() {}

    //! Full constructor
    UnprecParWilsonLinOp(Handle< FermState<T,P,Q> > fs,
			 const Real& Mass_, const Real& H_)
    {create(fs,Mass_,H_);}

    //! Destructor is automatic
    ~UnprecParWilsonLinOp() {}

    //! Only defined on the odd subset
    const Subset& subset() const {return all;}

    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const {return D.getFermBC();}

    //! Creation routine
    void create(Handle< FermState<T,P,Q> > fs,
		const Real& Mass_, const Real& H_);

    //! Apply the operator onto a source vector
    void operator() (LatticeFermion& chi, const LatticeFermion& psi, enum PlusMinus isign) const;

    //! Derivative of unpreconditioned ParWilson dM/dU
    void deriv(multi1d<LatticeColorMatrix>& ds_u, 
	       const LatticeFermion& chi, const LatticeFermion& psi, 
	       enum PlusMinus isign) const;

  private:
    Real Mass;
    Real H;
//    multi1d<LatticeColorMatrix> u;
    WilsonDslash D;
  };

} // End Namespace Chroma


#endif
