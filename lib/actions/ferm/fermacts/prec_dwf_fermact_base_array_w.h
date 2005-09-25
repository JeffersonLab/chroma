// -*- C++ -*-
// $Id: prec_dwf_fermact_base_array_w.h,v 2.0 2005-09-25 21:04:26 edwards Exp $
/*! \file
 *  \brief Base class for even-odd preconditioned domain-wall-like fermion actions
 */

#ifndef __prec_dwf_fermact_base_array_w_h__
#define __prec_dwf_fermact_base_array_w_h__

#include "fermact.h"
#include "actions/ferm/linop/prec_dwflike_linop_base_array_w.h"
#include "actions/ferm/linop/unprec_dwflike_linop_base_array_w.h"
#include "actions/ferm/linop/unprec_ppdwf4d_linop_w.h"
#include "actions/ferm/linop/lDeltaLs_w.h"
#include "actions/ferm/linop/lmdagm.h"


namespace Chroma
{
  //! Base class for unpreconditioned domain-wall-like fermion actions
  /*! \ingroup fermacts
   *
   * Unprecondition domain-wall fermion action. The conventions used here
   * are specified in Phys.Rev.D63:094505,2001 (hep-lat/0005002).
   */
  template<typename T, typename P>
  class EvenOddPrecDWFermActBaseArray : public EvenOddPrecWilsonTypeFermAct5D<T,P>
  {
  public:
    //! Return the quark mass
    virtual Real getQuarkMass() const = 0;

    //! Produce an unpreconditioned linear operator for this action with arbitrary quark mass
    virtual const UnprecDWLikeLinOpBaseArray<T,P>* unprecLinOp(Handle<const ConnectState> state, 
							       const Real& m_q) const = 0;

    //! Produce an even-odd preconditioned linear operator for this action with arbitrary quark mass
    virtual const EvenOddPrecDWLikeLinOpBaseArray<T,P>* precLinOp(Handle<const ConnectState> state, 
								  const Real& m_q) const = 0;

    //! Override to produce a DWF-link even-odd prec. linear operator for this action
    /*! Covariant return rule - override base class function */
    virtual const EvenOddPrecDWLikeLinOpBaseArray<T,P>* linOp(Handle<const ConnectState> state) const
    {
      return precLinOp(state,getQuarkMass());
    }

    //! Produce a linear operator M^dag.M for this action
    virtual const LinearOperator< multi1d<T> >* lMdagM(Handle<const ConnectState> state) const
    {
      return new lmdagm< multi1d<T> >(linOp(state));
    }

    //! Override to produce a DWF-link even-odd prec. Pauli-Villars linear operator for this action
    /*! Covariant return rule - override base class function */
    virtual const EvenOddPrecDWLikeLinOpBaseArray<T,P>* linOpPV(Handle<const ConnectState> state) const
    {
      return precLinOp(state,Real(1));
    }

    //! Produce a hermitian version of the linear operator
    /*! This code is generic */
    virtual const LinearOperator< multi1d<T> >* hermitianLinOp(Handle<const ConnectState> state) const
    {
      // Have not implemented this yet, but it is generic
      //
      // BALINT:
      // Oh no it isn't or is it(?)... what is the reflection matrix
      // for Moebius with Zolo coeffs where D(1) != D(2) etc.
      // However the function is at least virtual so it can be 
      // overridden on an as needed basis
      QDPIO::cerr << "EvenOddPrecDWFermActBaseArray::hermitianLinOp not implemented" << endl;
      QDP_abort(1);
      return 0;
    }

    //! Produce an unpreconditioned linear operator projecting 5D to 4D (the inverse of qprop below)
    /*! Use the fact that  linOp4D(m_q) = [P^{-1} (D^{(5)}(1))^{-1} D^{(5)}(m_q) P]_{11} */
    virtual const LinearOperator<T>* linOp4D(Handle<const ConnectState> state,
					     const Real& m_q,
					     const InvertParam_t& invParam) const
    {
      return new UnprecPPDWF4DLinOp<T,P>(precLinOp(state,m_q),
					 precLinOp(state,Real(1)),
					 invParam);
    }

    //! Produce a  DeltaLs = 1-epsilon^2(H) operator
    virtual const LinearOperator<T>* DeltaLs(Handle< const ConnectState> state,
					     const InvertParam_t& invParam) const 
    {
      Handle< const LinearOperator<T> >  lin(linOp4D(state,Real(0),invParam));
      return new lDeltaLs(lin);
    }

    //! Define quark propagator routine for 4D fermions
    /*! Default implementation provided */
    const SystemSolver<T>* qprop(Handle<const ConnectState> state,
				 const InvertParam_t& invParam) const;

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
	Handle< const EvenOddPrecDWLikeLinOpBaseArray<T,P> > A(linOp(state));
	A->Dminus(chi,psi,isign,s5);
      }

  };

}

#endif
