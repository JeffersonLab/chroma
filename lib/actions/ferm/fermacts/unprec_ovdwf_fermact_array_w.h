// -*- C++ -*-
// $Id: unprec_ovdwf_fermact_array_w.h,v 1.1 2003-11-15 03:23:06 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Overlap-DWF (Borici) action
 */

#ifndef __unprec_ovdwf_fermact_array_w_h__
#define __unprec_ovdwf_fermact_array_w_h__

#include "fermact_array.h"
#include "actions/ferm/linop/unprec_ovdwf_fermact_array_w.h"

using namespace QDP;

//! Unpreconditioned Overlap-style (Borica) DWF fermion action
/*! \ingroup fermact
 *
 * Unprecondition domain-wall fermion action. The conventions used here
 * are specified in Phys.Rev.D63:094505,2001 (hep-lat/0005002).
 */

class UnprecOvDWFermActArray : public UnprecDWFermActBaseArray
{
public:
  //! Partial constructor
  UnprecOvDWFermActArray() {}

  //! Full constructor
  UnprecOvDWFermActArray(const Real& WilsonMass_, const Real& m_q_, int N5_)
    {create(WilsonMass_, m_q_, N5_);}

  //! Creation routine
  void create(const Real& WilsonMass_, const Real& m_q_, int N5_);

  //! Length of DW flavor index/space
  int size() const {return N5;}

  //! Produce a linear operator for this action
  const LinearOperator< multi1d<LatticeFermion> >* linOp(const multi1d<LatticeColorMatrix>& u) const;

  //! Produce a linear operator M^dag.M for this action
  const LinearOperator< multi1d<LatticeFermion> >* lMdagM(const multi1d<LatticeColorMatrix>& u) const;

  //! Produce a linear operator for this action but with quark mass 1
  const LinearOperator< multi1d<LatticeFermion> >* linOpPV(const multi1d<LatticeColorMatrix>& u) const;

  //! Destructor is automatic
  ~UnprecOvDWFermActArray() {}

private:
  Real WilsonMass;
  Real m_q;
  Real a5;
  int  N5;
};

#endif
