// -*- C++ -*-
// $Id: unprec_ovext_fermact_array_w.h,v 1.4 2003-12-02 15:45:04 edwards Exp $
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
  //! Partial constructor
  UnprecOvExtFermActArray() {}

  //! Full constructor
  UnprecOvExtFermActArray(const Real& WilsonMass_, const Real& m_q_, int N5_)
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
	     const ConnectState& state, 
	     const LatticeFermion& chi, 
	     enum InvType invType,
	     const Real& RsdCG, 
	     int MaxCG, int& ncg_had) const;

  //! Destructor is automatic
  ~UnprecOvExtFermActArray() {}

private:
  Real WilsonMass;
  Real m_q;
  Real a5;
  int  N5;
};

#endif
