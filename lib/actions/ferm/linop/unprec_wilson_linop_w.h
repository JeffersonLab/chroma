// -*- C++ -*-
// $Id: unprec_wilson_linop_w.h,v 1.8 2003-11-20 05:43:41 edwards Exp $
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

 *      M  =  (d+M) - (1/2) D'
 *
 */

class UnprecWilsonLinOp : public LinearOperator<LatticeFermion>
{
public:
  //! Partial constructor
  UnprecWilsonLinOp() {}

  //! Full constructor
  UnprecWilsonLinOp(const multi1d<LatticeColorMatrix>& u_, const Real& Mass_)
    {create(u_,Mass_);}

  //! Destructor is automatic
  ~UnprecWilsonLinOp() {}

  //! Only defined on the odd subset
  const OrderedSubset& subset() const {return all;}

  //! Creation routine
  void create(const multi1d<LatticeColorMatrix>& u_, const Real& Mass_);

  //! Apply the operator onto a source vector
  void operator() (LatticeFermion& chi, const LatticeFermion& psi, enum PlusMinus isign) const;

private:
  Real Mass;
  multi1d<LatticeColorMatrix> u;
  WilsonDslash D;
};

#endif
