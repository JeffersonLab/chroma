// -*- C++ -*-
// $Id: prec_wilson_fermact_w.h,v 1.5 2004-01-23 10:35:36 bjoo Exp $
/*! \file
 *  \brief Even-odd preconditioned Wilson fermion action
 */

#ifndef __prec_wilson_fermact_w_h__
#define __prec_wilson_fermact_w_h__

#include "fermact.h"
#include "actions/ferm/linop/lgherm_w.h"

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
  //! General FermBC
  EvenOddPrecWilsonFermAct(Handle< FermBC<LatticeFermion> > fbc_, 
		      const Real& Mass_) : 
    fbc(fbc_), Mass(Mass_) {}

  //! Copy constructor
  EvenOddPrecWilsonFermAct(const EvenOddPrecWilsonFermAct& a) : 
    fbc(a.fbc), Mass(a.Mass) {}

  //! Assignment
  EvenOddPrecWilsonFermAct& operator=(const EvenOddPrecWilsonFermAct& a)
    {fbc=a.fbc; Mass=a.Mass; return *this;}

  //! Return the fermion BC object for this action
  const FermBC<LatticeFermion>& getFermBC() const {return *fbc;}

  //! Produce a linear operator for this action
  const EvenOddPrecLinearOperator<LatticeFermion>* linOp(Handle<const ConnectState> state) const;

  //! Produce a linear operator M^dag.M for this action
  const LinearOperator<LatticeFermion>* lMdagM(Handle<const ConnectState> state) const;

  //! Produce the gamma_5 hermitian operator H_w
  const LinearOperator<LatticeFermion>* gamma5HermLinOp(Handle< const ConnectState> state) const { 
    return new lgherm<LatticeFermion>(linOp(state));
  }

  //! Override - compute dS_f/dU
  void dsdu(multi1d<LatticeColorMatrix>& result,
	    Handle<const ConnectState> state,
	    const LatticeFermion& psi) const;

  //! Destructor is automatic
  ~EvenOddPrecWilsonFermAct() {}

private:
  Handle< FermBC<LatticeFermion> >  fbc;
  Real Mass;
};

#endif
