// -*- C++ -*-
// $Id: unprec_ovext_fermact_array_w.h,v 1.6 2004-01-23 20:38:17 edwards Exp $
/*! \file
 *  \brief Unpreconditioned extended-Overlap (5D) (Naryanan&Neuberger) action
 */

#ifndef __unprec_ovext_fermact_array_w_h__
#define __unprec_ovext_fermact_array_w_h__

#include "fermact.h"

using namespace QDP;

//! Unpreconditioned Extended-Overlap (N&N) linear operator
/*!
 * \ingroup fermact
 *
 * This operator applies the extended version of the hermitian overlap operator
 *   Chi  =   ((1+m_q)/(1-m_q)*gamma_5 + B) . Psi
 *  where  B  is the continued fraction of the pole approx. to eps(H(m))
 */

class UnprecOvExtFermActArray : public UnprecWilsonTypeFermAct< multi1d<LatticeFermion> >
{
public:
  //! General FermBC
  UnprecOvExtFermActArray(Handle< FermBC< multi1d<LatticeFermion> > > fbc_, 
		       const Real& WilsonMass_, const Real& m_q_, int N5_) : 
    fbc(fbc_), WilsonMass(WilsonMass_), m_q(m_q_), N5(N5_) {a5=1;}

  //! Copy constructor
  UnprecOvExtFermActArray(const UnprecOvExtFermActArray& a) : 
    fbc(a.fbc), WilsonMass(a.WilsonMass), m_q(a.m_q), a5(a.a5), N5(a.N5) {}

  //! Assignment
  UnprecOvExtFermActArray& operator=(const UnprecOvExtFermActArray& a)
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
  const LinearOperator< multi1d<LatticeFermion> >* gamma5HermLinOp(Handle<const ConnectState> state) const
    {
      QDPIO::cerr << "UnprecOvExtFermActArray::gamma5HermLinOp not implemented" << endl;
      QDP_abort(1);
    }

  //! Compute quark propagator over base type
  /*! 
   * Solves  M.psi = chi
   *
   * \param psi      quark propagator ( Modify )
   * \param u        gauge field ( Read )
   * \param chi      source ( Modify )
   * \param invType  inverter type ( Read (
   * \param RsdCG    CG (or MR) residual used here ( Read )
   * \param MaxCG    maximum number of CG iterations ( Read )
   * \param ncg_had  number of CG iterations ( Write )
   *
   * NOTE: maybe this should produce a quark prop foundry class object 
   */
  void qprop(LatticeFermion& psi, 
	     Handle<const ConnectState> state, 
	     const LatticeFermion& chi, 
	     enum InvType invType,
	     const Real& RsdCG, 
	     int MaxCG, int& ncg_had) const;

  //! Destructor is automatic
  ~UnprecOvExtFermActArray() {}

private:
  // Hide partial constructor
  UnprecOvExtFermActArray() {}

private:
  Handle< FermBC< multi1d<LatticeFermion> > >  fbc;
  Real WilsonMass;
  Real m_q;
  Real a5;
  int  N5;
};

#endif
