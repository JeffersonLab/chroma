// -*- C++ -*-
// $Id: prec_dwf_fermact_base_array_w.h,v 1.19 2004-12-29 22:13:40 edwards Exp $
/*! \file
 *  \brief Base class for even-odd preconditioned domain-wall-like fermion actions
 */

#ifndef __prec_dwf_fermact_base_array_w_h__
#define __prec_dwf_fermact_base_array_w_h__

#include "fermact.h"
#include "actions/ferm/linop/prec_dwf_linop_base_array_w.h"
#include "actions/ferm/linop/unprec_dwf_linop_base_array_w.h"
#include "actions/ferm/linop/unprec_ppdwf4d_linop_w.h"
#include "actions/ferm/linop/lDeltaLs_w.h"
#include "actions/ferm/linop/llincomb.h"
#include "actions/ferm/linop/lmdagm.h"

using namespace QDP;

namespace Chroma
{
  //! Base class for unpreconditioned domain-wall-like fermion actions
  /*! \ingroup fermact
   *
   * Unprecondition domain-wall fermion action. The conventions used here
   * are specified in Phys.Rev.D63:094505,2001 (hep-lat/0005002).
   */
  template<typename T>
  class EvenOddPrecDWFermActBaseArray : public EvenOddPrecWilsonTypeFermAct5D<T, multi1d<LatticeColorMatrix> >
  {
  public:
    //! Return the quark mass
    virtual Real quark_mass() const = 0;

    //! Produce an unpreconditioned linear operator for this action with arbitrary quark mass
    virtual const UnprecDWLinOpBaseArray<T>* unprecLinOp(Handle<const ConnectState> state, 
							 const Real& m_q) const = 0;

    //! Produce an even-odd preconditioned linear operator for this action with arbitrary quark mass
    virtual const EvenOddPrecDWLinOpBaseArray<T>* precLinOp(Handle<const ConnectState> state, 
							    const Real& m_q) const = 0;

    //! Override to produce a DWF-link even-odd prec. linear operator for this action
    /*! Covariant return rule - override base class function */
    virtual const EvenOddPrecDWLinOpBaseArray<T>* linOp(Handle<const ConnectState> state) const
    {
      return precLinOp(state,quark_mass());
    }

    //! Produce a linear operator M^dag.M for this action
    virtual const LinearOperator< multi1d<LatticeFermion> >* lMdagM(Handle<const ConnectState> state) const
    {
      return new lmdagm< multi1d<LatticeFermion> >(linOp(state));
    }

    //! Produce a hermitian version of the linear operator
    /*! This code is generic */
    virtual const LinearOperator< multi1d<T> >* gamma5HermLinOp(Handle<const ConnectState> state) const
    {
      // Have not implemented this yet, but it is generic
      QDPIO::cerr << "EvenOddPrecDWFermActBaseArray::gamma5HermLinOp not implemented" << endl;
      QDP_abort(1);
      return 0;
    }

    //! Produce an unpreconditioned linear operator projecting 5D to 4D (the inverse of qprop below)
    /*! Use the fact that  linOp4D(m_q) = [P^{-1} (D^{(5)}(1))^{-1} D^{(5)}(m_q) P]_{11} */
    virtual const LinearOperator<T>* linOp4D(Handle<const ConnectState> state,
					     const Real& m_q,
					     const InvertParam_t& invParam) const
    {
      return new UnprecPPDWF4DLinOp<T>(precLinOp(state,m_q),
				       precLinOp(state,Real(1)),
				       invParam);
    }

    //! Produce a  DeltaLs = 1-epsilon^2(H) operator
    virtual const LinearOperator<LatticeFermion>* DeltaLs(Handle< const ConnectState> state,
							  const InvertParam_t& invParam) const 
    {
      Handle< const LinearOperator<LatticeFermion> >  lin(linOp4D(state,Real(0),invParam));
      return new lDeltaLs(lin);
    }

    //! Define quark propagator routine for 4D fermions
    void qprop(T& psi, 
	       Handle<const ConnectState> state, 
	       const T& chi, 
	       const InvertParam_t& invParam,
	       int& ncg_had) const;

    //! Apply the Dminus operator on a fermion.
    /*! 
     * Slightly more than a convenience function, 
     * it avoids specifying the type of the linOp. 
     * Used in the dwf_quarkProp routine.
     */
    void Dminus(T& chi,
		const T& psi,
		Handle<const ConnectState> state,
		enum PlusMinus isign,
		int s5) const
      {
	Handle< const EvenOddPrecDWLinOpBaseArray<T> > A(linOp(state));
	A->Dminus(chi,psi,isign,s5);
      }

  };

}

#endif
