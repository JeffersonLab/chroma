// -*- C++ -*-
// $Id: unprec_dwf_fermact_base_array_w.h,v 1.14 2004-09-19 02:39:45 edwards Exp $
/*! \file
 *  \brief Base class for unpreconditioned domain-wall-like fermion actions
 */

#ifndef __unprec_dwf_fermact_base_array_w_h__
#define __unprec_dwf_fermact_base_array_w_h__

#include "fermact.h"
#include "actions/ferm/linop/unprec_dwf_linop_base_array_w.h"

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
  class UnprecDWFermActBaseArray : public UnprecWilsonTypeFermAct< multi1d<T> >
  {
  public:
    //! Return the quark mass
    virtual Real quark_mass() const = 0;

    //! Override to produce a DWF-link unprec. linear operator for this action
    /*! Covariant return rule - override base class function */
    virtual const UnprecDWLinOpBaseArray<T>* linOp(Handle<const ConnectState> state) const = 0;

    //! Produce a hermitian version of the linear operator
    /*! This code is generic */
    virtual const LinearOperator< multi1d<T> >* gamma5HermLinOp(Handle<const ConnectState> state) const
      {
	// Have not implemented this yet, but it is generic
	QDPIO::cerr << "UnprecDWFermActBaseArray::gamma5HermLinOp not implemented" << endl;
	QDP_abort(1);
	return 0;
      }

    //! Produce a linear operator for this action but with quark mass 1
    virtual const UnprecDWLinOpBaseArray<T>* linOpPV(Handle<const ConnectState> state) const = 0;

    //! Define quark propagator routine for 4D fermions
    void qprop(T& psi, 
	       Handle<const ConnectState> state, 
	       const T& chi, 
	       const InvertParam_t& invParam,
	       int& ncg_had) const;

    //! Given a complete propagator as a source, this does all the inversions needed
    /*! \ingroup qprop
     *
     * This routine is actually generic to Domain Wall fermions (Array) fermions
     *
     * \param q_sol    quark propagator ( Write )
     * \param q_src    source ( Read )
     * \param xml_out  diagnostic output ( Modify )
     * \param state    gauge connection state ( Read )
     * \param t_src    time slice of source ( Read )
     * \param j_decay  direction of decay ( Read )
     * \param invParam inverter parameters ( Read )
     * \param ncg_had  number of CG iterations ( Write )
     */
    virtual 
    void dwf_quarkProp4(LatticePropagator& q_sol, 
			XMLWriter& xml_out,
			const LatticePropagator& q_src,
			int t_src, int j_decay,
			Handle<const ConnectState> state,
			const InvertParam_t& invParam,
			int& ncg_had)
      {
	// Simply call corresponding quarkProp4 routine
	// Assumes nonRelProp = false
	quarkProp4(q_sol, xml_out, q_src, state, invParam, false, ncg_had);
      }
    
  
    //! Apply the Dminus operator on a vector in Ls.
    /*! 
     * Slightly more than a convenience function, 
     * it avoids specifying the type of the linOp. 
     * Used in the dwf_quarkProp4 routine.
     */
    void Dminus(multi1d<T>& chi,
		const multi1d<T>& psi,
		Handle<const ConnectState> state,
		enum PlusMinus isign) const
      {
	Handle< const UnprecDWLinOpBaseArray<T> > A(linOp(state));
	A->Dminus(chi,psi,isign);
      }

    //! Apply the Dminus operator on a fermion.
    /*! 
     * Slightly more than a convenience function, 
     * it avoids specifying the type of the linOp. 
     * Used in the dwf_quarkProp4 routine.
     */
    void Dminus(T& chi,
		const T& psi,
		Handle<const ConnectState> state,
		enum PlusMinus isign) const
      {
	Handle< const UnprecDWLinOpBaseArray<T> > A(linOp(state));
	A->Dminus(chi,psi,isign);
      }

  };

}

using namespace Chroma;

#endif
