// -*- C++ -*-
// $Id: unprec_parwilson_fermact_w.h,v 1.2 2004-01-23 10:35:36 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned Wilson fermion action with parity breaking term
 */

#ifndef __unprec_parwilson_fermact_w_h__
#define __unprec_parwilson_fermact_w_h__

#include "fermact.h"

using namespace QDP;

//! Unpreconditioned Wilson fermion action with parity breaking term
/*! \ingroup fermact
 *
 * Supports creation and application for fermion actions
 *
 * The kernel for Wilson fermions with a parity breaking term is
 *
 *      M  =  (d+M) + i*H*gamma_5  - (1/2) D'
 */

class UnprecParWilsonFermAct : public UnprecWilsonTypeFermAct<LatticeFermion>
{
public:
  //! General FermBC
  UnprecParWilsonFermAct(Handle< FermBC<LatticeFermion> > fbc_, 
			 const Real& Mass_, const Real& H_) : 
    fbc(fbc_), Mass(Mass_), H(H_) {}

  //! Copy constructor
  UnprecParWilsonFermAct(const UnprecParWilsonFermAct& a) : 
    fbc(a.fbc), Mass(a.Mass), H(a.H) {}

  //! Assignment
  UnprecParWilsonFermAct& operator=(const UnprecParWilsonFermAct& a)
    {fbc=a.fbc; Mass=a.Mass; H=a.H; return *this;}

  //! Return the fermion BC object for this action
  const FermBC<LatticeFermion>& getFermBC() const {return *fbc;}

  //! Produce a linear operator for this action
  const LinearOperator<LatticeFermion>* linOp(Handle<const ConnectState> state) const;

  //! Produce a linear operator M^dag.M for this action
  const LinearOperator<LatticeFermion>* lMdagM(Handle<const ConnectState> state) const;

 //! Produce the gamma_5 hermitin op gamma_5 M
  const LinearOperator<LatticeFermion>* gamma5HermLinOp(Handle<const ConnectState> state) const {

    QDP_error_exit("gamma5HermLinOp not implemented yet for this action\n");
    return 0;
  }


  //! Compute dS_f/dU
  void dsdu(multi1d<LatticeColorMatrix>& result,
	    Handle<const ConnectState> state,
	    const LatticeFermion& psi) const;

  //! Destructor is automatic
  ~UnprecParWilsonFermAct() {}

private:
  UnprecParWilsonFermAct() {} //hide default constructor
  
private:
  Handle< FermBC<LatticeFermion> >  fbc;
  Real Mass;
  Real H;
};

#endif
