// -*- C++ -*-
// $Id: unprec_ovdwf_fermact_array_w.h,v 1.5 2004-01-02 03:19:40 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Overlap-DWF (Borici) action
 */

#ifndef __unprec_ovdwf_fermact_array_w_h__
#define __unprec_ovdwf_fermact_array_w_h__

#include "fermact.h"
#include "actions/ferm/fermacts/unprec_dwf_fermact_base_array_w.h"

using namespace QDP;

//! Unpreconditioned Overlap-style (Borici) DWF fermion action
/*! \ingroup fermact
 *
 * Unprecondition domain-wall fermion action. The conventions used here
 * are specified in Phys.Rev.D63:094505,2001 (hep-lat/0005002).
 */

class UnprecOvDWFermActArray : public UnprecDWFermActBaseArray
{
public:
  //! General FermBC
  UnprecOvDWFermActArray(Handle< FermBC< multi1d<LatticeFermion> > > fbc_, 
			 const Real& WilsonMass_, const Real& m_q_, int N5_) : 
    fbc(fbc_), WilsonMass(WilsonMass_), m_q(m_q_), N5(N5_) {a5=1;}

  //! Copy constructor
  UnprecOvDWFermActArray(const UnprecOvDWFermActArray& a) : 
    fbc(a.fbc), WilsonMass(a.WilsonMass), m_q(a.m_q), a5(a.a5), N5(a.N5) {}

  //! Assignment
  UnprecOvDWFermActArray& operator=(const UnprecOvDWFermActArray& a)
    {fbc=a.fbc; WilsonMass=a.WilsonMass; m_q=a.m_q; a5=a.a5; N5=a.N5; return *this;}

  //! Return the fermion BC object for this action
  const FermBC< multi1d<LatticeFermion> >& getFermBC() const {return *fbc;}

  //! Length of DW flavor index/space
  int size() const {return N5;}

  //! Return the quark mass
  Real quark_mass() const {return m_q;}

  //! Produce a linear operator for this action
  const LinearOperator< multi1d<LatticeFermion> >* linOp(Handle<const ConnectState> state) const;

  //! Produce a linear operator M^dag.M for this action
  const LinearOperator< multi1d<LatticeFermion> >* lMdagM(Handle<const ConnectState> state) const;

  //! Produce a linear operator for this action but with quark mass 1
  const LinearOperator< multi1d<LatticeFermion> >* linOpPV(Handle<const ConnectState> state) const;

  //! Destructor is automatic
  ~UnprecOvDWFermActArray() {}

private:
  Handle< FermBC< multi1d<LatticeFermion> > >  fbc;
  Real WilsonMass;
  Real m_q;
  Real a5;
  int  N5;
};

#endif
