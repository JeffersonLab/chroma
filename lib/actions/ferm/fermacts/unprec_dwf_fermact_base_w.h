// -*- C++ -*-
// $Id: unprec_dwf_fermact_base_w.h,v 1.6 2004-01-23 17:59:07 edwards Exp $
/*! \file
 *  \brief Base class for unpreconditioned domain-wall fermion action
 */

#ifndef __unprec_dwf_fermact_base_w_h__
#define __unprec_dwf_fermact_base_w_h__

#include "fermact.h"

using namespace QDP;

//! Base class for unpreconditioned domain-wall fermion action
/*! \ingroup fermact
 *
 * Unprecondition domain-wall fermion action. The conventions used here
 * are specified in Phys.Rev.D63:094505,2001 (hep-lat/0005002).
 */

class UnprecDWFermActBase : public UnprecWilsonTypeFermAct<LatticeDWFermion>
{
public:
  //! Return the quark mass
  virtual Real quark_mass() const = 0;

  //! Produce a hermitian version of the linear operator
  /*! This code is generic */
  virtual const LinearOperator<LatticeDWFermion>* gamma5HermLinOp(Handle<const ConnectState> state) const
    {
      // Have not implemented this yet, but it is generic
      QDPIO::cerr << "UnprecDWFermActBase::gamma5HermLinOp not implemented" << endl;
      QDP_abort(1);
    }

  //! Produce a linear operator for this action but with quark mass 1
  virtual const LinearOperator<LatticeDWFermion>* linOpPV(Handle<const ConnectState> state) const = 0;

  //! Redefine quark propagator routine for 4D fermions
  void qprop(LatticeFermion& psi, 
	     Handle<const ConnectState> state, 
	     const LatticeFermion& chi, 
	     enum InvType invType,
	     const Real& RsdCG, 
	     int MaxCG, int& ncg_had) const;
};

#endif
