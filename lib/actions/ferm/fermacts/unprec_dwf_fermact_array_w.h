// -*- C++ -*-
// $Id: unprec_dwf_fermact_array_w.h,v 1.9 2004-03-17 03:37:22 edwards Exp $
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

class UnprecDWFermActArray : public UnprecDWFermActBaseArray<LatticeFermion>
{
public:
  //! General FermBC
  UnprecDWFermActArray(Handle< FermBC< multi1d<LatticeFermion> > > fbc_, 
		       const Real& WilsonMass_, const Real& m_q_, int N5_) : 
    fbc(fbc_), WilsonMass(WilsonMass_), m_q(m_q_), N5(N5_) {a5=1;}

  //! Copy constructor
  UnprecDWFermActArray(const UnprecDWFermActArray& a) : 
    fbc(a.fbc), WilsonMass(a.WilsonMass), m_q(a.m_q), a5(a.a5), N5(a.N5) {}

  //! Assignment
  UnprecDWFermActArray& operator=(const UnprecDWFermActArray& a)
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

  //! Produce a hermitian version of the linear operator
  /*! This code is generic */
  const LinearOperator< multi1d<LatticeFermion> >* gamma5HermLinOp(Handle<const ConnectState> state) const
    {
      // Have not implemented this yet, but it is generic
      QDPIO::cerr << "UnprecDWFermActBaseArray::gamma5HermLinOp not implemented" << endl;
      QDP_abort(1);
      return 0;
    }

  //! Produce a linear operator for this action but with quark mass 1
  const LinearOperator< multi1d<LatticeFermion> >* linOpPV(Handle<const ConnectState> state) const;

  //! Destructor is automatic
  ~UnprecDWFermActArray() {}

private:
  Handle< FermBC< multi1d<LatticeFermion> > >  fbc;
  Real WilsonMass;
  Real m_q;
  Real a5;
  int  N5;
};

#endif
