// -*- C++ -*-
// $Id: unprec_dwf_fermact_w.h,v 1.6 2003-11-15 04:25:42 edwards Exp $
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
  //! Partial constructor
  UnprecDWFermAct() {}

  //! Full constructor
  UnprecDWFermAct(const Real& WilsonMass, const Real& m_q)
    {create(WilsonMass, m_q);}

  //! Creation routine
  void create(const Real& WilsonMass, const Real& m_q);

  //! Return the quark mass
  Real quark_mass() const {return m_q;}

  //! Produce a linear operator for this action
  const LinearOperator<LatticeDWFermion>* linOp(const multi1d<LatticeColorMatrix>& u) const;

  //! Produce a linear operator M^dag.M for this action
  const LinearOperator<LatticeDWFermion>* lMdagM(const multi1d<LatticeColorMatrix>& u) const;

  //! Produce a linear operator for this action but with quark mass 1
  const LinearOperator<LatticeDWFermion>* linOpPV(const multi1d<LatticeColorMatrix>& u) const;

  //! Destructor is automatic
  ~UnprecDWFermAct() {}

private:
  Real WilsonMass;
  Real m_q;
  Real a5;
};

#endif
