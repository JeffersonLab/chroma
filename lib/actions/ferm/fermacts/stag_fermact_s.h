// -*- C++ -*-
// $Id: stag_fermact_s.h,v 1.2 2004-03-29 21:32:28 edwards Exp $
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

class StagFermAct : public EvenOddStaggeredTypeFermAct<LatticeFermion>
{
public:
  //! General FermBC
  StagFermAct(Handle< FermBC<LatticeFermion> > fbc_, 
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
  ~StagFermAct() {}

private:
  StagFermAct() {} //hide default constructor
  
private:
  Handle< FermBC<LatticeFermion> >  fbc;
  Real Mass;
};

#endif
