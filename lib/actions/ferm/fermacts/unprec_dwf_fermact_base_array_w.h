// -*- C++ -*-
// $Id: unprec_dwf_fermact_base_array_w.h,v 3.2 2006-10-19 16:01:29 edwards Exp $
/*! \file
 *  \brief Base class for unpreconditioned domain-wall-like fermion actions
 */

#ifndef __unprec_dwf_fermact_base_array_w_h__
#define __unprec_dwf_fermact_base_array_w_h__

#include "unprec_wilstype_fermact_w.h"
#include "actions/ferm/linop/unprec_dwflike_linop_base_array_w.h"
#include "actions/ferm/linop/unprec_dwf4d_linop_w.h"
#include "actions/ferm/linop/lDeltaLs_w.h"
#include "actions/ferm/linop/llincomb.h"
#include "actions/ferm/invert/syssolver_cg_params.h"

 
namespace Chroma
{
  //! Base class for unpreconditioned domain-wall-like fermion actions
  /*! \ingroup fermacts
   *
   * Unprecondition domain-wall fermion action. The conventions used here
   * are specified in Phys.Rev.D63:094505,2001 (hep-lat/0005002).
   */
  template<typename T, typename P, typename Q>
  class UnprecDWFermActBaseArray : public UnprecWilsonTypeFermAct5D<T,P,Q>
  {
  public:
    //! Return the quark mass
    virtual Real getQuarkMass() const = 0;

    //! Produce an unpreconditioned linear operator for this action with arbitrary quark mass
    virtual UnprecDWLikeLinOpBaseArray<T,P,Q>* unprecLinOp(Handle< FermState<T,P,Q> > state, 
							   const Real& m_q) const = 0;

    //! Override to produce a DWF-link unprec. linear operator for this action
    /*! Covariant return rule - override base class function */
    virtual UnprecDWLikeLinOpBaseArray<T,P,Q>* linOp(Handle< FermState<T,P,Q> > state) const
    {
      return unprecLinOp(state,getQuarkMass());
    }

    //! Override to produce a DWF-link unprec. Pauli-Villars linear operator for this action
    /*! Covariant return rule - override base class function */
    virtual UnprecDWLikeLinOpBaseArray<T,P,Q>* linOpPV(Handle< FermState<T,P,Q> > state) const
    {
      return unprecLinOp(state,Real(1));
    }

    //! Produce a hermitian version of the linear operator
    /*! This code is generic */
    virtual LinearOperatorArray<T>* hermitianLinOp(Handle< FermState<T,P,Q> > state) const
    {
      // Have not implemented this yet, but it is generic
      QDPIO::cerr << "UnprecDWFermActBaseArray::gamma5HermLinOp not implemented" << endl;
      QDP_abort(1);
      return 0;
    }

    //! Produce an unpreconditioned linear operator projecting 5D to 4D (the inverse of qprop below)
    /*! Use the fact that  linOp4D(m_q) = [P^{-1} (D^{(5)}(1))^{-1} D^{(5)}(m_q) P]_{11} */
    virtual LinearOperator<T>* linOp4D(Handle< FermState<T,P,Q> > state,
				       const Real& m_q,
				       const GroupXML_t& invParam) const
    {
      std::istringstream  is(invParam.xml);
      XMLReader  paramtop(is);
	
      return new UnprecDWF4DLinOp<T>(unprecLinOp(state,m_q),
				     unprecLinOp(state,Real(1)),
				     SysSolverCGParams(paramtop,invParam.path));
    }

    //! Produce a  DeltaLs = 1-epsilon^2(H) operator
    virtual LinearOperator<T>* DeltaLs(Handle< FermState<T,P,Q> > state,
				       const GroupXML_t& invParam) const 
    {
      Handle< LinearOperator<T> >  lin(linOp4D(state,Real(0),invParam));
      return new lDeltaLs(lin);
    }

    //! Define quark propagator routine for 4D fermions
    /*! Default implementation provided */
    SystemSolver<T>* qprop(Handle< FermState<T,P,Q> > state,
			   const GroupXML_t& invParam) const;

    //! Apply the Dminus operator on a fermion.
    /*! 
     * Slightly more than a convenience function, 
     * it avoids specifying the type of the linOp. 
     * Used in the dwf_quarkProp routine.
     */
    void Dminus(T& chi,
		const T& psi,
		Handle< FermState<T,P,Q> > state,
		enum PlusMinus isign,
		int s5) const
      {
	Handle< UnprecDWLikeLinOpBaseArray<T,P,Q> > A(linOp(state));
	A->Dminus(chi,psi,isign,s5);
      }

  };

}


#endif
