// -*- C++ -*-
// $Id: unprec_dwf_fermact_w.h,v 1.8 2004-01-02 03:19:40 edwards Exp $
/*! \file
 *  \brief Unpreconditioned domain-wall fermion action
 */

#ifndef __unprec_dwf_fermact_w_h__
#define __unprec_dwf_fermact_w_h__

#include "fermact.h"
#include "actions/ferm/fermacts/unprec_dwf_fermact_base_w.h"

using namespace QDP;

//! Unpreconditioned domain-wall fermion action
/*! \ingroup fermact
 *
 * Unprecondition domain-wall fermion action. The conventions used here
 * are specified in Phys.Rev.D63:094505,2001 (hep-lat/0005002).
 */

class UnprecDWFermAct : public UnprecDWFermActBase
{
public:
  //! General FermBC
  UnprecDWFermAct(Handle< FermBC<LatticeDWFermion> > fbc_, 
		  const Real& WilsonMass_, const Real& m_q_) : 
    fbc(fbc_), WilsonMass(WilsonMass_), m_q(m_q_) {a5=1;}

  //! Copy constructor
  UnprecDWFermAct(const UnprecDWFermAct& a) : 
    fbc(a.fbc), WilsonMass(a.WilsonMass), m_q(a.m_q), a5(a.a5) {}

  //! Assignment
  UnprecDWFermAct& operator=(const UnprecDWFermAct& a)
    {fbc=a.fbc; WilsonMass=a.WilsonMass; m_q=a.m_q; a5=a.a5; return *this;}

  //! Return the fermion BC object for this action
  const FermBC<LatticeDWFermion>& getFermBC() const {return *fbc;}

  //! Return the quark mass
  Real quark_mass() const {return m_q;}

  //! Produce a linear operator for this action
  const LinearOperator<LatticeDWFermion>* linOp(Handle<const ConnectState> state) const;

  //! Produce a linear operator M^dag.M for this action
  const LinearOperator<LatticeDWFermion>* lMdagM(Handle<const ConnectState> state) const;

  //! Produce a linear operator for this action but with quark mass 1
  const LinearOperator<LatticeDWFermion>* linOpPV(Handle<const ConnectState> state) const;

  //! Destructor is automatic
  ~UnprecDWFermAct() {}

private:
  Handle< FermBC<LatticeDWFermion> >  fbc;
  Real WilsonMass;
  Real m_q;
  Real a5;
};

#endif
