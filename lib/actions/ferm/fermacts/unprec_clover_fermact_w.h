// -*- C++ -*-
// $Id: unprec_clover_fermact_w.h,v 1.1 2003-11-22 21:33:24 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Clover fermion action
 */

#ifndef __unprec_clover_fermact_w_h__
#define __unprec_clover_fermact_w_h__

#include "fermact.h"
#include "actions/ferm/linop/unprec_clover_linop_w.h"

using namespace QDP;

//! Unpreconditioned Clover fermion action
/*! \ingroup fermact
 *
 * Unpreconditioned clover fermion action
 */

class UnprecCloverFermAct : public UnprecWilsonTypeFermAct<LatticeFermion>
{
public:
  //! Partial constructor
  UnprecCloverFermAct() {}

  //! Full constructor
  /*! Isotropic action */
  UnprecCloverFermAct(const Real& Mass_, const Real& ClovCoeff_, const Real& u0_)
    {create(Mass_, ClovCoeff_, u0_);}

  //! Creation routine
  void create(const Real& Mass_, const Real& ClovCoeff_, const Real& u0_);

  //! Produce a linear operator for this action
  const LinearOperator<LatticeFermion>* linOp(const multi1d<LatticeColorMatrix>& u) const;

  //! Produce a linear operator M^dag.M for this action
  const LinearOperator<LatticeFermion>* lMdagM(const multi1d<LatticeColorMatrix>& u) const;

  //! Override - compute dS_f/dU
  void dsdu(multi1d<LatticeColorMatrix>& result,
	    const multi1d<LatticeColorMatrix>& u,
	    const LatticeFermion& psi) const;

  //! Destructor is automatic
  ~UnprecCloverFermAct() {}

private:
  Real Mass;
  Real ClovCoeff;
  Real u0;
};

#endif
