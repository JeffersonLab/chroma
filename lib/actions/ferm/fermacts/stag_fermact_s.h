// -*- C++ -*-
// $Id: stag_fermact_s.h,v 1.3 2004-11-06 13:00:24 mcneile Exp $
/*! \file
 *  \brief Staggered fermion action
 */

#ifndef __stag_fermact_w_h__
#define __stag_fermact_w_h__

#include "fermact.h"

using namespace QDP;

//! Staggered fermion action
/*! \ingroup fermact
 *
 */

class StagFermAct : public EvenOddStaggeredTypeFermAct<LatticeStaggeredFermion>
{
public:
  //! General FermBC
  StagFermAct(Handle< FermBC<LatticeStaggeredFermion> > fbc_, 
	      const Real& Mass_) : 
    fbc(fbc_), Mass(Mass_) {}

  //! Copy constructor
  StagFermAct(const StagFermAct& a) : 
    fbc(a.fbc), Mass(a.Mass) {}

  //! Assignment
  StagFermAct& operator=(const StagFermAct& a)
  {fbc=a.fbc; Mass=a.Mass; return *this;}

  //! Return the quark mass
  const Real getQuarkMass() const {return Mass;}

  //! Return the fermion BC object for this action
  const FermBC<LatticeStaggeredFermion>& getFermBC() const {return *fbc;}

  //! Produce a linear operator for this action
  const LinearOperator<LatticeStaggeredFermion>* linOp(Handle<const ConnectState> state) const;

  //! Produce a linear operator M^dag.M for this action
  const LinearOperator<LatticeStaggeredFermion>* lMdagM(Handle<const ConnectState> state) const;

  //! Compute dS_f/dU
  void dsdu(multi1d<LatticeColorMatrix>& result,
	    Handle<const ConnectState> state,
	    const LatticeStaggeredFermion& psi) const;

  //! Destructor is automatic
  ~StagFermAct() {}

private:
  StagFermAct() {} //hide default constructor
  
private:
  Handle< FermBC<LatticeStaggeredFermion> >  fbc;
  Real Mass;
};

#endif
