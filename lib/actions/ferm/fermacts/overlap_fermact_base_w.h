// -*- C++ -*-
// $Id: overlap_fermact_base_w.h,v 3.2 2006-10-19 16:01:27 edwards Exp $
/*! \file
 *  \brief Base class for unpreconditioned overlap-like fermion actions
 */

#ifndef __overlap_fermact_base_w_h__
#define __overlap_fermact_base_w_h__

#include "unprec_wilstype_fermact_w.h"
#include "meas/eig/ischiral_w.h"
#include "actions/ferm/linop/lgherm_w.h"
#include "io/enum_io/enum_inner_solver_type_io.h"
#include "actions/ferm/linop/lDeltaLs_w.h"


namespace Chroma
{

  //! Base class for unpreconditioned overlap-like fermion actions
  /*! \ingroup fermacts
   *
   * Unpreconditioned overlap-like fermion action. 
   * The conventions used here are specified in some Nucl.Phys.B. article
   * by Edwards,Heller, Narayanan
   *
   * NOTE: for now we assume the kernel is a fund. rep. fermion type,
   * but that is not necessary
   */
  class OverlapFermActBase : public UnprecWilsonTypeFermAct<LatticeFermion, 
	   multi1d<LatticeColorMatrix>,
	   multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Virtual copy constructor
    virtual OverlapFermActBase* clone() const = 0;

    //! Return the quark mass
    virtual Real getQuarkMass() const = 0;

    //! Does this object really satisfy the Ginsparg-Wilson relation?
    virtual bool isChiral() const = 0;

    //! Produce an unpreconditioned linear operator for this action with arbitrary quark mass
    virtual UnprecLinearOperator<T,P,Q>* unprecLinOp(Handle< FermState<T,P,Q> > state, 
						     const Real& m_q) const = 0;
    
    //! Override to produce a DWF-link unprec. linear operator for this action
    /*! Covariant return rule - override base class function */
    virtual UnprecLinearOperator<T,P,Q>* linOp(Handle< FermState<T,P,Q> > state) const
      {
	return unprecLinOp(state,getQuarkMass());
      }

    //! Robert's way: 
    //! Produce a linear operator M^dag.M for this action to be applied
    //  to a vector of known chirality. Chirality is passed in
    virtual DiffLinearOperator<T,P,Q>* lMdagM(Handle< FermState<T,P,Q> > state) const = 0;

    virtual DiffLinearOperator<T,P,Q>* lMdagM(Handle< FermState<T,P,Q> > state, 
					      const Chirality& chirality) const = 0;

    virtual LinearOperator<T>* hermitianLinOp(Handle< FermState<T,P,Q> > state) const 
      { 
	return new lgherm<T>(linOp(state));
      }

    virtual LinearOperator<T>* DeltaLs(Handle< FermState<T,P,Q> > state) const {
      Handle< LinearOperator<T> > lin(unprecLinOp(state,Real(0)));
      return new lDeltaLs(lin);
    }

    //! Produce a linear operator that gives back gamma_5 eps(H)
    virtual LinearOperator<T>* lgamma5epsH(Handle< FermState<T,P,Q> > state) const = 0;

   
    //! Produce a linear operator that gives back gamma_5 eps(H)
    virtual LinearOperator<T>* lgamma5epsHPrecondition(Handle< FermState<T,P,Q> > state) const = 0;

    virtual LinearOperator<T>* linOpPrecondition(Handle< FermState<T,P,Q> > state) const = 0;

    //! Robert's way: 
    //! Produce a linear operator M^dag.M for this action to be applied
    //  to a vector of known chirality. Chirality is passed in
    virtual LinearOperator<T>* lMdagMPrecondition(Handle< FermState<T,P,Q> > state) const = 0;

    virtual LinearOperator<T>* lMdagMPrecondition(Handle< FermState<T,P,Q> > state, 
						  const Chirality& chirality) const = 0;

    //! Redefine quark propagator routine for 4D fermions
    /*! Default implementation provided */
    SystemSolver<T>* qprop(Handle< FermState<T,P,Q> > state,
			   const GroupXML_t& invParam) const;
  
    //! Define a multi mass qprop
    /*! this should be possible for most 4D operators of the 
     *  form   (1/2)[ (1 + Mass ) + (1 - Mass) gamma_5 psi 
     *
     */
    void multiQprop(multi1d<T>& psi, 
		    const multi1d<Real>& masses,
		    Handle< FermState<T,P,Q> > state,
		    const T& chi,
		    const GroupXML_t& invParam,
		    const int n_soln,
		    int & ncg_had) const;


    //! Define a multi mass qprop
    /*!
     * \param q_sol    quark propagator ( Write )
     * \param q_src    source ( Read )
     * \param invParam inverter parameters ( Read )
     * \param ncg_had  number of CG iterations ( Write )
     */
    void multiQuarkProp(multi1d<LatticePropagator>& q_sol, 
			XMLWriter& xml_out,
			const LatticePropagator& q_src,
			Handle< FermState<T,P,Q> > state,
			const multi1d<Real>& masses,
			const GroupXML_t& invParam,
			const int n_soln,
			int& ncg_had);
  };

}

#endif
