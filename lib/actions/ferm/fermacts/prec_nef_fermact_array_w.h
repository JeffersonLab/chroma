// -*- C++ -*-
// $Id: prec_nef_fermact_array_w.h,v 1.1 2004-08-23 21:07:47 kostas Exp $
/*! \file
 *  \brief 4D style even-odd preconditioned NEF fermion action
 */

#ifndef __prec_nef_fermact_array_w_h__
#define __prec_nef_fermact_array_w_h__

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

class EvenOddPrecNEFFermActArray : public EvenOddPrecDWFermActBaseArray<LatticeFermion>
{
public:
  //! General FermBC
  EvenOddPrecNEFFermActArray(Handle< FermBC< multi1d<LatticeFermion> > > fbc_, 
			    const Real& WilsonMass_, 
			     const Real& b5_, const Real& c5_,
			     const Real& m_q_, int N5_) : 
    fbc(fbc_), WilsonMass(WilsonMass_),  b5(b5_),c5(c5_),  m_q(m_q_), N5(N5_) 
    {
      QDPIO::cout << "Construct EvenOddPrecNEFFermActArray: WilsonMass = " << WilsonMass 
		  << "  m_q = " << m_q 
		  << "  N5 = " << N5 
		  << "  b5 = " << b5
		  << "  c5 = " << c5
		  << endl;
    }

  //! Copy constructor
  EvenOddPrecNEFFermActArray(const EvenOddPrecNEFFermActArray& a) : 
    fbc(a.fbc), WilsonMass(a.WilsonMass),  b5(a.b5), c5(a.c5), 
    m_q(a.m_q), N5(a.N5) {}

  //! Assignment
  EvenOddPrecNEFFermActArray& operator=(const EvenOddPrecNEFFermActArray& a)
     {
       fbc=a.fbc; 
       WilsonMass=a.WilsonMass; 
       b5=a.b5; 
       c5=a.c5; 
       m_q=a.m_q; 
       N5=a.N5; 
       return *this;
     }

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
  ~EvenOddPrecNEFFermActArray() {}

private:
  Handle< FermBC< multi1d<LatticeFermion> > >  fbc;
  Real WilsonMass;
  Real b5;
  Real c5;
  Real m_q;
  int  N5;
};

#endif
