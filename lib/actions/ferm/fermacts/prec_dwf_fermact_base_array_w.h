// -*- C++ -*-
// $Id: prec_dwf_fermact_base_array_w.h,v 1.6 2004-01-23 17:59:07 edwards Exp $
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

template<typename T>
class EvenOddPrecDWFermActBaseArray : public EvenOddPrecWilsonTypeFermAct< multi1d<T> >
{
public:
  //! Return the quark mass
  virtual Real quark_mass() const = 0;

  //! Produce a hermitian version of the linear operator
  /*! This code is generic */
  virtual const LinearOperator< multi1d<T> >* gamma5HermLinOp(Handle<const ConnectState> state) const
    {
      // Have not implemented this yet, but it is generic
      QDPIO::cerr << "EvenOddPrecDWFermActBaseArray::gamma5HermLinOp not implemented" << endl;
      QDP_abort(1);
    }

  //! Produce a linear operator for this action but with quark mass 1
  virtual const LinearOperator< multi1d<T> >* linOpPV(Handle<const ConnectState> state) const = 0;

  //! Define quark propagator routine for 4D fermions
  void qprop(T& psi, 
	     Handle<const ConnectState> state, 
	     const T& chi, 
	     enum InvType invType,
	     const Real& RsdCG, 
	     int MaxCG, int& ncg_had) const;
};

#endif
