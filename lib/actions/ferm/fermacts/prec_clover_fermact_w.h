// -*- C++ -*-
// $Id: prec_clover_fermact_w.h,v 1.2 2003-12-02 15:45:04 edwards Exp $
/*! \file
 *  \brief Even-odd preconditioned Clover fermion action
 */

#ifndef __prec_clover_fermact_w_h__
#define __prec_clover_fermact_w_h__

#include "fermact.h"
#include "actions/ferm/linop/prec_clover_linop_w.h"

using namespace QDP;

//! Even-odd preconditioned Clover fermion action
/*! \ingroup fermact
 *
 * Even-odd preconditioned clover fermion action. 
 * Only defined on odd subset.
 */

class EvenOddPrecCloverFermAct : public EvenOddPrecWilsonTypeFermAct<LatticeFermion>
{
public:
  //! Partial constructor
  EvenOddPrecCloverFermAct() {}

  //! Full constructor
  /*! Isotropic action */
  EvenOddPrecCloverFermAct(const Real& Mass_, const Real& ClovCoeff_, const Real& u0_)
    {create(Mass_, ClovCoeff_, u0_);}

  //! Creation routine
  void create(const Real& Mass_, const Real& ClovCoeff_, const Real& u0_);

  //! Produce a linear operator for this action
  const EvenOddPrecLinearOperator<LatticeFermion> linOp(const ConnectState& state) const;

  //! Produce a linear operator M^dag.M for this action
  const LinearOperator<LatticeFermion> lMdagM(const ConnectState& state) const;

  //! Override - compute dS_f/dU
  void dsdu(multi1d<LatticeColorMatrix>& result,
	    const ConnectState& state,
	    const LatticeFermion& psi) const;

  //! Destructor is automatic
  ~EvenOddPrecCloverFermAct() {}

private:
  Real Mass;
  Real ClovCoeff;
  Real u0;
};

#endif
