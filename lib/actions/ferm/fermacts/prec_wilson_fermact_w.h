// -*- C++ -*-
// $Id: prec_wilson_fermact_w.h,v 1.1 2003-11-22 21:34:01 edwards Exp $
/*! \file
 *  \brief Even-odd preconditioned Wilson fermion action
 */

#ifndef __prec_wilson_fermact_w_h__
#define __prec_wilson_fermact_w_h__

#include "fermact.h"
#include "actions/ferm/linop/prec_wilson_linop_w.h"

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
  const EvenOddPrecLinearOperator<LatticeFermion>* linOp(const multi1d<LatticeColorMatrix>& u) const;

  //! Produce a linear operator M^dag.M for this action
  const LinearOperator<LatticeFermion>* lMdagM(const multi1d<LatticeColorMatrix>& u) const;

  //! Override - compute dS_f/dU
  void dsdu(multi1d<LatticeColorMatrix>& result,
	    const multi1d<LatticeColorMatrix>& u,
	    const LatticeFermion& psi) const;

  //! Destructor is automatic
  ~EvenOddPrecWilsonFermAct() {}

private:
  Real Mass;
};

#endif
