// -*- C++ -*-
// $Id: prec_dwf_fermact_base_array_w.h,v 1.3 2003-12-02 15:45:04 edwards Exp $
/*! \file
 *  \brief Base class for even-odd preconditioned domain-wall-like fermion actions
 */

#ifndef __prec_dwf_fermact_base_array_w_h__
#define __prec_dwf_fermact_base_array_w_h__

#include "fermact.h"

using namespace QDP;

//! Base class for unpreconditioned domain-wall-like fermion actions
/*! \ingroup fermact
 *
 * Unprecondition domain-wall fermion action. The conventions used here
 * are specified in Phys.Rev.D63:094505,2001 (hep-lat/0005002).
 */

class EvenOddPrecDWFermActBaseArray : public EvenOddPrecWilsonTypeFermAct< multi1d<LatticeFermion> >
{
public:
  //! Return the quark mass
  virtual Real quark_mass() const = 0;

  //! Produce a linear operator for this action but with quark mass 1
  virtual const LinearOperator< multi1d<LatticeFermion> >* linOpPV(const ConnectState& state) const = 0;

  //! Define quark propagator routine for 4D fermions
  void qprop(LatticeFermion& psi, 
	     const ConnectState& state, 
	     const LatticeFermion& chi, 
	     enum InvType invType,
	     const Real& RsdCG, 
	     int MaxCG, int& ncg_had) const;
};

#endif
