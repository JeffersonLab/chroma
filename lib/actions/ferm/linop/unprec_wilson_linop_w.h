// -*- C++ -*-
// $Id: unprec_wilson_linop_w.h,v 1.5 2003-10-20 20:31:50 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Wilson fermion linear operator
 */

#ifndef __unprec_wilson_linop_w_h__
#define __unprec_wilson_linop_w_h__

#include "linearop.h"
#include "actions/ferm/linop/dslash_w.h"

using namespace QDP;

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

 * Algorithm:

 * The kernel for Wilson fermions is

 *      M  =  1 - k D'
 *
 */

class UnprecWilsonLinOp : public LinearOperator<LatticeFermion>
{
public:
  //! Partial constructor
  UnprecWilsonLinOp() {}

  //! Full constructor
  UnprecWilsonLinOp(const multi1d<LatticeColorMatrix>& _u, const Real& _Kappa)
    {create(_u,_Kappa);}

  //! Destructor is automatic
  ~UnprecWilsonLinOp() {}

  //! Only defined on the odd subset
  const OrderedSubset& subset() const {return all;}

  //! Creation routine
  void create(const multi1d<LatticeColorMatrix>& _u, const Real& _Kappa);

  //! Apply the operator onto a source vector
  LatticeFermion operator() (const LatticeFermion& psi, enum LinOpSign isign) const;

private:
  Real Kappa;
  multi1d<LatticeColorMatrix> u;
  WilsonDslash D;
};

#endif
