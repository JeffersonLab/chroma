// -*- C++ -*-
// $Id: overlap_fermact_base_w.h,v 1.14 2004-10-15 16:37:19 bjoo Exp $
/*! \file
 *  \brief Base class for unpreconditioned overlap-like fermion actions
 */

#ifndef __overlap_fermact_base_w_h__
#define __overlap_fermact_base_w_h__

#include "fermact.h"
#include "meas/eig/ischiral_w.h"
#include "actions/ferm/linop/lgherm_w.h"
#include "io/enum_io/enum_inner_solver_type_io.h"
#include "actions/ferm/linop/lDeltaLs_w.h"

using namespace QDP;

namespace Chroma
{

  //! Base class for unpreconditioned overlap-like fermion actions
  /*! \ingroup fermact
   *
   * Unpreconditioned overlap-like fermion action. 
   * The conventions used here are specified in some Nucl.Phys.B. article
   * by Edwards,Heller, Narayanan
   *
   * NOTE: for now we assume the kernel is a fund. rep. fermion type,
   * but that is not necessary
   */
  class OverlapFermActBase : public UnprecWilsonTypeFermAct<LatticeFermion>
  {
  public:
    //! Return the quark mass
    virtual Real quark_mass() const = 0;

    //! Does this object really satisfy the Ginsparg-Wilson relation?
    virtual bool isChiral() const = 0;

    //! Robert's way: 
    //! Produce a linear operator M^dag.M for this action to be applied
    //  to a vector of known chirality. Chirality is passed in
    virtual const LinearOperator<LatticeFermion>* lMdagM(Handle<const ConnectState> state) const = 0;

    virtual const LinearOperator<LatticeFermion>* lMdagM(Handle<const ConnectState> state, const Chirality& chirality) const = 0;

    virtual const LinearOperator<LatticeFermion>* gamma5HermLinOp(Handle< const ConnectState> state) const { 
      return new lgherm<LatticeFermion>(linOp(state));
    }

    virtual const LinearOperator<LatticeFermion>* DeltaLs(Handle< const ConnectState> state) const {
      Handle< const LinearOperator<LatticeFermion> > g5eps=lgamma5epsH(state);
      return new lDeltaLs(g5eps);
    }

    //! Produce a linear operator that gives back gamma_5 eps(H)
    virtual const LinearOperator<LatticeFermion>* lgamma5epsH(Handle<const ConnectState> state) const = 0;

   
    //! Produce a linear operator that gives back gamma_5 eps(H)
    virtual const LinearOperator<LatticeFermion>* lgamma5epsHPrecondition(Handle<const ConnectState> state) const = 0;

    virtual const LinearOperator<LatticeFermion>* linOpPrecondition(Handle<const ConnectState> state) const = 0;

    //! Robert's way: 
    //! Produce a linear operator M^dag.M for this action to be applied
    //  to a vector of known chirality. Chirality is passed in
    virtual const LinearOperator<LatticeFermion>* lMdagMPrecondition(Handle<const ConnectState> state) const = 0;

    virtual const LinearOperator<LatticeFermion>* lMdagMPrecondition(Handle<const ConnectState> state, const Chirality& chirality) const = 0;

    //! Redefine quark propagator routine for 4D fermions
    /*! 
     * NOTE: the arg ConectState MUST be in the original base because C++ 
     * requires it for a virtual func!
     * The function will have to downcast to get the correct state
     */
    void qprop(LatticeFermion& psi, 
	       Handle<const ConnectState> state, 
	       const LatticeFermion& chi, 
	       const InvertParam_t& invParam,
	       int& ncg_had) const;

  
    //! Define a multi mass qprop
    /*! this should be possible for most 4D operators of the 
     *  form   (1/2)[ (1 + Mass ) + (1 - Mass) gamma_5 psi 
     *
     */
    void multiQprop(multi1d<LatticeFermion>& psi, 
		    const multi1d<Real>& masses,
		    Handle<const ConnectState> state,
		    const LatticeFermion& chi,
		    const MultiInvertParam_t& invParam,
		    const int n_soln,
		    int & ncg_had) const;


    void quarkProp4(LatticePropagator& q_sol, 
		    XMLWriter& xml_out,
		    const LatticePropagator& q_src,
		    Handle<const ConnectState> state,
		    const InvertParam_t& invParam,
		    bool nonRelProp,
		    int& ncg_had);
    /*! \ingroup qprop
     *
     * \param q_sol    quark propagator ( Write )
     * \param q_src    source ( Read )
     * \param invParam inverter parameters ( Read )
     * \param ncg_had  number of CG iterations ( Write )
     */
    void multiQuarkProp4(multi1d<LatticePropagator>& q_sol, 
			 XMLWriter& xml_out,
			 const LatticePropagator& q_src,
			 Handle<const ConnectState> state,
			 const multi1d<Real>& masses,
			 const MultiInvertParam_t& invParam,
			 const int n_soln,
			 int& ncg_had);
  };

}

#endif
