// -*- C++ -*-
// $Id: prec_clover_fermact_w.h,v 1.5 2004-01-23 10:35:36 bjoo Exp $
/*! \file
 *  \brief Even-odd preconditioned Clover fermion action
 */

#ifndef __prec_clover_fermact_w_h__
#define __prec_clover_fermact_w_h__

#include "fermact.h"
#include "actions/ferm/linop/lgherm_w.h"

using namespace QDP;

//! Even-odd preconditioned Clover fermion action
/*! \ingroup fermact
 *
 * Even-odd preconditioned clover fermion action. 
 * Only defined on odd subset.
 */

class EvenOddPrecCloverFermAct : public EvenOddPrecWilsonTypeFermAct<LatticeFermion>
{
public:
  //! General FermBC
  EvenOddPrecCloverFermAct(Handle< FermBC<LatticeFermion> > fbc_, 
			   const Real& Mass_, const Real& ClovCoeff_, const Real& u0_) : 
    fbc(fbc_), Mass(Mass_), ClovCoeff(ClovCoeff_), u0(u0_) {}

  //! Copy constructor
  EvenOddPrecCloverFermAct(const EvenOddPrecCloverFermAct& a) : 
    fbc(a.fbc), Mass(a.Mass), ClovCoeff(a.ClovCoeff), u0(a.u0) {}

  //! Assignment
  EvenOddPrecCloverFermAct& operator=(const EvenOddPrecCloverFermAct& a)
    {fbc=a.fbc; Mass=a.Mass; ClovCoeff=a.ClovCoeff; u0=a.u0; return *this;}

  //! Return the fermion BC object for this action
  const FermBC<LatticeFermion>& getFermBC() const {return *fbc;}

  //! Produce a linear operator for this action
  const EvenOddPrecLinearOperator<LatticeFermion> linOp(Handle<const ConnectState> state) const;

  //! Produce a linear operator M^dag.M for this action
  const LinearOperator<LatticeFermion> lMdagM(Handle<const ConnectState> state) const;

  const LinearOperator<LatticeFermion>* gamma5HermLinOp(Handle< const ConnectState> state) const { 
    return new lgherm<LatticeFermion>(linOp(state));
  }

  

  //! Override - compute dS_f/dU
  void dsdu(multi1d<LatticeColorMatrix>& result,
	    Handle<const ConnectState> state,
	    const LatticeFermion& psi) const;

  //! Destructor is automatic
  ~EvenOddPrecCloverFermAct() {}

private:
  Handle< FermBC<LatticeFermion> >  fbc;
  Real Mass;
  Real ClovCoeff;
  Real u0;
};

#endif
