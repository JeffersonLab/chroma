// -*- C++ -*-
// $Id: unprec_wilson_fermact_w.h,v 1.15 2004-01-02 03:19:41 edwards Exp $
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

class UnprecWilsonFermAct : public UnprecWilsonTypeFermAct<LatticeFermion>
{
public:
  //! General FermBC
  UnprecWilsonFermAct(Handle< FermBC<LatticeFermion> > fbc_, 
		      const Real& Mass_) : 
    fbc(fbc_), Mass(Mass_) {}

  //! Copy constructor
  UnprecWilsonFermAct(const UnprecWilsonFermAct& a) : 
    fbc(a.fbc), Mass(a.Mass) {}

  //! Assignment
  UnprecWilsonFermAct& operator=(const UnprecWilsonFermAct& a)
    {fbc=a.fbc; Mass=a.Mass; return *this;}

  //! Return the fermion BC object for this action
  const FermBC<LatticeFermion>& getFermBC() const {return *fbc;}

  //! Produce a linear operator for this action
  const LinearOperator<LatticeFermion>* linOp(Handle<const ConnectState> state) const;

  //! Produce a linear operator M^dag.M for this action
  const LinearOperator<LatticeFermion>* lMdagM(Handle<const ConnectState> state) const;

  //! Compute dS_f/dU
  void dsdu(multi1d<LatticeColorMatrix>& result,
	    Handle<const ConnectState> state,
	    const LatticeFermion& psi) const;

  //! Destructor is automatic
  ~UnprecWilsonFermAct() {}

private:
  UnprecWilsonFermAct() {} //hide default constructor
  
private:
  Handle< FermBC<LatticeFermion> >  fbc;
  Real Mass;
};

#endif
