// -*- C++ -*-
// $Id: unprec_clover_fermact_w.h,v 1.5 2004-01-02 03:19:40 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Clover fermion action
 */

#ifndef __unprec_clover_fermact_w_h__
#define __unprec_clover_fermact_w_h__

#include "fermact.h"

using namespace QDP;

//! Unpreconditioned Clover fermion action
/*! \ingroup fermact
 *
 * Unpreconditioned clover fermion action
 */

class UnprecCloverFermAct : public UnprecWilsonTypeFermAct<LatticeFermion>
{
public:
  //! General FermBC
  /*! Isotropic action */
  UnprecCloverFermAct(Handle< FermBC<LatticeFermion> > fbc_, 
			   const Real& Mass_, const Real& ClovCoeff_, const Real& u0_) : 
    fbc(fbc_), Mass(Mass_), ClovCoeff(ClovCoeff_), u0(u0_) {}

  //! Copy constructor
  UnprecCloverFermAct(const UnprecCloverFermAct& a) : 
    fbc(a.fbc), Mass(a.Mass), ClovCoeff(a.ClovCoeff), u0(a.u0) {}

  //! Assignment
  UnprecCloverFermAct& operator=(const UnprecCloverFermAct& a)
    {fbc=a.fbc; Mass=a.Mass; ClovCoeff=a.ClovCoeff; u0=a.u0; return *this;}

  //! Return the fermion BC object for this action
  const FermBC<LatticeFermion>& getFermBC() const {return *fbc;}

  //! Produce a linear operator for this action
  const LinearOperator<LatticeFermion> linOp(Handle<const ConnectState> state) const;

  //! Produce a linear operator M^dag.M for this action
  const LinearOperator<LatticeFermion> lMdagM(Handle<const ConnectState> state) const;

  //! Destructor is automatic
  ~UnprecCloverFermAct() {}

private:
  // Hide partial constructor
  UnprecCloverFermAct() {}

private:
  Handle< FermBC<LatticeFermion> >  fbc;
  Real Mass;
  Real ClovCoeff;
  Real u0;
};

#endif
