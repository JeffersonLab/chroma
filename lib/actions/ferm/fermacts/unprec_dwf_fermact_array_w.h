// -*- C++ -*-
// $Id: unprec_dwf_fermact_array_w.h,v 1.5 2003-12-02 15:45:04 edwards Exp $
/*! \file
 *  \brief Unpreconditioned domain-wall fermion action
 */

#ifndef __unprec_dwf_fermact_array_w_h__
#define __unprec_dwf_fermact_array_w_h__

#include "fermact.h"
#include "actions/ferm/fermacts/unprec_dwf_fermact_base_array_w.h"

using namespace QDP;

//! Unpreconditioned domain-wall fermion action
/*! \ingroup fermact
 *
 * Unprecondition domain-wall fermion action. The conventions used here
 * are specified in Phys.Rev.D63:094505,2001 (hep-lat/0005002).
 */

class UnprecDWFermActArray : public UnprecDWFermActBaseArray
{
public:
  //! Partial constructor
  UnprecDWFermActArray() {}

  //! Full constructor
  UnprecDWFermActArray(const Real& WilsonMass_, const Real& m_q_, int N5_)
    {create(WilsonMass_, m_q_, N5_);}

  //! Creation routine
  void create(const Real& WilsonMass_, const Real& m_q_, int N5_);

  //! Length of DW flavor index/space
  int size() const {return N5;}

  //! Return the quark mass
  Real quark_mass() const {return m_q;}

  //! Produce a linear operator for this action
  const LinearOperator< multi1d<LatticeFermion> >* linOp(const ConnectState& state) const;

  //! Produce a linear operator M^dag.M for this action
  const LinearOperator< multi1d<LatticeFermion> >* lMdagM(const ConnectState& state) const;

  //! Produce a linear operator for this action but with quark mass 1
  const LinearOperator< multi1d<LatticeFermion> >* linOpPV(const ConnectState& state) const;

  //! Destructor is automatic
  ~UnprecDWFermActArray() {}

private:
  Real WilsonMass;
  Real m_q;
  Real a5;
  int  N5;
};

#endif
