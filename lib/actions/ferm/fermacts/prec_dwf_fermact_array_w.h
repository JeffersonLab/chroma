// -*- C++ -*-
// $Id: prec_dwf_fermact_array_w.h,v 1.1 2003-11-23 05:54:14 edwards Exp $
/*! \file
 *  \brief 4D style even-odd preconditioned domain-wall fermion action
 */

#ifndef __prec_dwf_fermact_array_w_h__
#define __prec_dwf_fermact_array_w_h__

#include "fermact.h"
#include "actions/ferm/fermacts/prec_dwf_fermact_base_array_w.h"

using namespace QDP;

//! 4D style even-odd preconditioned domain-wall fermion action
/*! \ingroup fermact
 *
 * 4D style even-odd preconditioned domain-wall fermion action. 
 * Follows notes of Orginos (10/2003)
 *
 * Hopefully, the conventions used here
 * are specified in Phys.Rev.D63:094505,2001 (hep-lat/0005002).
 */

class EvenOddPrecDWFermActArray : public EvenOddPrecDWFermActBaseArray
{
public:
  //! Partial constructor
  EvenOddPrecDWFermActArray() {}

  //! Full constructor
  EvenOddPrecDWFermActArray(const Real& WilsonMass_, const Real& m_q_, int N5_)
    {create(WilsonMass_, m_q_, N5_);}

  //! Creation routine
  void create(const Real& WilsonMass_, const Real& m_q_, int N5_);

  //! Length of DW flavor index/space
  int size() const {return N5;}

  //! Return the quark mass
  Real quark_mass() const {return m_q;}

  //! Produce a linear operator for this action
  const EvenOddPrecLinearOperator< multi1d<LatticeFermion> >* linOp(const multi1d<LatticeColorMatrix>& u) const;

  //! Produce a linear operator M^dag.M for this action
  const LinearOperator< multi1d<LatticeFermion> >* lMdagM(const multi1d<LatticeColorMatrix>& u) const;

  //! Produce a linear operator for this action but with quark mass 1
  const LinearOperator< multi1d<LatticeFermion> >* linOpPV(const multi1d<LatticeColorMatrix>& u) const;

  //! Destructor is automatic
  ~EvenOddPrecDWFermActArray() {}

private:
  Real WilsonMass;
  Real m_q;
  Real a5;
  int  N5;
};

#endif
