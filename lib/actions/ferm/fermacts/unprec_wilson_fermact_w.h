// -*- C++ -*-
// $Id: unprec_wilson_fermact_w.h,v 1.8 2003-10-10 03:46:46 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Wilson fermion action
 */

#ifndef __unprec_wilson_fermact_w_h__
#define __unprec_wilson_fermact_w_h__

#include "fermact.h"
#include "actions/ferm/linop/unprec_wilson_linop_w.h"

using namespace QDP;

//! Unpreconditioned Wilson fermion action
/*! \ingroup fermact
 *
 * Supports creation and application for fermion actions
 */

class UnprecWilsonFermAct : public UnprecWilsonTypeFermAct 
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

#if 1
// THIS SHOULD NOT BE NEEDED
  //! Compute quark propagator
  void qprop(LatticeFermion& psi, 
	     const multi1d<LatticeColorMatrix>& u, 
	     const LatticeFermion& chi, 
	     const Real& RsdCG, 
	     int MaxCG, int& ncg_had) const;
#endif

  //! Compute dS_f/dU
  multi1d<LatticeColorMatrix> dsdu(const multi1d<LatticeColorMatrix>& u,
				   const LatticeFermion& psi) const;

  //! Destructor is automatic
  ~UnprecWilsonFermAct() {}

private:
  Real Kappa;
};

#endif
