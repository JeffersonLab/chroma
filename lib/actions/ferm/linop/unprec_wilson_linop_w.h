// -*- C++ -*-
// $Id: unprec_wilson_linop_w.h,v 1.10 2005-01-14 20:13:06 edwards Exp $
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
  
  class UnprecWilsonLinOp : public UnprecLinearOperator<LatticeFermion, multi1d<LatticeColorMatrix> >
  {
  public:
    //! Partial constructor
    UnprecWilsonLinOp() {}

    //! Full constructor
    UnprecWilsonLinOp(const multi1d<LatticeColorMatrix>& u_, const Real& Mass_)
      {create(u_,Mass_);}

    //! Destructor is automatic
    ~UnprecWilsonLinOp() {}

    //! Creation routine
    void create(const multi1d<LatticeColorMatrix>& u_, const Real& Mass_);

    //! Apply the operator onto a source vector
    void operator() (LatticeFermion& chi, const LatticeFermion& psi, enum PlusMinus isign) const;

    //! Derivative of unpreconditioned Wilson dM/dU
    void deriv(multi1d<LatticeColorMatrix>& ds_u, 
	       const LatticeFermion& chi, const LatticeFermion& psi, 
	       enum PlusMinus isign) const;

  private:
    Real Mass;
    multi1d<LatticeColorMatrix> u;
    WilsonDslash D;
  };

}; // End Namespace Chroma


#endif
