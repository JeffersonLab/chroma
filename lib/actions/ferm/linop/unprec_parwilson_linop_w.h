// -*- C++ -*-
// $Id: unprec_parwilson_linop_w.h,v 1.1 2004-01-12 04:48:00 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Wilson fermion linear operator with parity breaking term
 */

#ifndef __unprec_parwilson_linop_w_h__
#define __unprec_parwilson_linop_w_h__

#include "linearop.h"
#include "actions/ferm/linop/dslash_w.h"

using namespace QDP;

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

class UnprecParWilsonLinOp : public LinearOperator<LatticeFermion>
{
public:
  //! Partial constructor
  UnprecParWilsonLinOp() {}

  //! Full constructor
  UnprecParWilsonLinOp(const multi1d<LatticeColorMatrix>& u_, 
		       const Real& Mass_, const Real& H_)
    {create(u_,Mass_,H_);}

  //! Destructor is automatic
  ~UnprecParWilsonLinOp() {}

  //! Only defined on the odd subset
  const OrderedSubset& subset() const {return all;}

  //! Creation routine
  void create(const multi1d<LatticeColorMatrix>& u_, 
	      const Real& Mass_, const Real& H_);

  //! Apply the operator onto a source vector
  void operator() (LatticeFermion& chi, const LatticeFermion& psi, enum PlusMinus isign) const;

private:
  Real Mass;
  Real H;
  multi1d<LatticeColorMatrix> u;
  WilsonDslash D;
};

#endif
