// -*- C++ -*-
// $Id: prec_dwf_fermact_array_w.h,v 1.5 2004-01-29 17:09:14 edwards Exp $
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

class EvenOddPrecDWFermActArray : public EvenOddPrecDWFermActBaseArray<LatticeFermion>
{
public:
  //! General FermBC
  EvenOddPrecDWFermActArray(Handle< FermBC< multi1d<LatticeFermion> > > fbc_, 
			    const Real& WilsonMass_, const Real& m_q_, int N5_) : 
    fbc(fbc_), WilsonMass(WilsonMass_), m_q(m_q_), N5(N5_) 
    {
      a5=1;
      QDPIO::cout << "EvenOddPrecDWFermActArray: WilsonMass = " << WilsonMass 
		  << "m_q = " << m_q 
		  << "N5 = " << N5 
		  << endl;
    }

  //! Copy constructor
  EvenOddPrecDWFermActArray(const EvenOddPrecDWFermActArray& a) : 
    fbc(a.fbc), WilsonMass(a.WilsonMass), m_q(a.m_q), a5(a.a5), N5(a.N5) {}

  //! Assignment
  EvenOddPrecDWFermActArray& operator=(const EvenOddPrecDWFermActArray& a)
    {fbc=a.fbc; WilsonMass=a.WilsonMass; m_q=a.m_q; a5=a.a5; N5=a.N5; return *this;}

  //! Return the fermion BC object for this action
  const FermBC< multi1d<LatticeFermion> >& getFermBC() const {return *fbc;}

  //! Length of DW flavor index/space
  int size() const {return N5;}

  //! Return the quark mass
  Real quark_mass() const {return m_q;}

  //! Produce a linear operator for this action
  const EvenOddPrecLinearOperator< multi1d<LatticeFermion> >* linOp(Handle<const ConnectState> state) const;

  //! Produce a linear operator M^dag.M for this action
  const LinearOperator< multi1d<LatticeFermion> >* lMdagM(Handle<const ConnectState> state) const;

  //! Produce a linear operator for this action but with quark mass 1
  const LinearOperator< multi1d<LatticeFermion> >* linOpPV(Handle<const ConnectState> state) const;

  //! Destructor is automatic
  ~EvenOddPrecDWFermActArray() {}

private:
  Handle< FermBC< multi1d<LatticeFermion> > >  fbc;
  Real WilsonMass;
  Real m_q;
  Real a5;
  int  N5;
};

#endif
