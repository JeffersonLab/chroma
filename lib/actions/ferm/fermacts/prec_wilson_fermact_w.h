// -*- C++ -*-
// $Id: prec_wilson_fermact_w.h,v 1.3 2003-12-02 21:30:35 edwards Exp $
/*! \file
 *  \brief Even-odd preconditioned Wilson fermion action
 */

#ifndef __prec_wilson_fermact_w_h__
#define __prec_wilson_fermact_w_h__

#include "fermact.h"

using namespace QDP;

//! Even-odd preconditioned Wilson fermion action
/*! \ingroup fermact
 *
 * Even-odd preconditioned wilson fermion action. 
 * Only defined on odd subset.
 */

class EvenOddPrecWilsonFermAct : public EvenOddPrecWilsonTypeFermAct<LatticeFermion>
{
public:
  //! Partial constructor
  EvenOddPrecWilsonFermAct() {}

  //! Full constructor
  /*! Isotropic action */
  EvenOddPrecWilsonFermAct(const Real& Mass_)
    {create(Mass_);}

  //! Creation routine
  void create(const Real& Mass_);

  //! Produce a linear operator for this action
  const EvenOddPrecLinearOperator<LatticeFermion>* linOp(const ConnectState& state) const;

  //! Produce a linear operator M^dag.M for this action
  const LinearOperator<LatticeFermion>* lMdagM(const ConnectState& state) const;

  //! Override - compute dS_f/dU
  void dsdu(multi1d<LatticeColorMatrix>& result,
	    const ConnectState& state,
	    const LatticeFermion& psi) const;

  //! Destructor is automatic
  ~EvenOddPrecWilsonFermAct() {}

private:
  Real Mass;
};

#endif
