// -*- C++ -*-
// $Id: unprec_wilson_fermact_w.h,v 1.1 2003-04-09 05:57:15 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Wilson fermion action
 */

#ifndef __unprec_wilson_fermact_w_h__
#define __unprec_wilson_fermact_w_h__

#include "fermact.h"
#include "actions/ferm/linop/unprec_wilson_linop_w.h"

using namespace QDP;

// Forward declaration
// class UnprecWilsonLinOp : public LinearOperator;

//! Unpreconditioned Wilson fermion action
/*! \ingroup fermact
 *
 * Supports creation and application for fermion actions
 */

class UnprecWilsonFermAct : UnprecWilsonTypeFermAct 
{
public:
  //! Partial constructor
  UnprecWilsonFermAct() {}

  //! Full constructor
  UnprecWilsonFermAct(const Real& _Kappa)
    {create(_Kappa);}

  //! Creation routine
  void create(const Real& _Kappa);

  //! Produce a linear operator for this action
  const LinearOperator* linOp(const multi1d<LatticeColorMatrix>& u) const;

  //! Produce a linear operator M^dag.M for this action
  const LinearOperator* lMdagM(const multi1d<LatticeColorMatrix>& u) const;

  //! Compute dS_f/dU
  multi1d<LatticeColorMatrix> dsdu(const multi1d<LatticeColorMatrix>& u) const;

  //! Destructor is automatic
  ~UnprecWilsonFermAct() {}

private:
  Real Kappa;
};

#endif
