// -*- C++ -*-
// $Id: stagtype_fermact_s.h,v 3.3 2008-09-08 16:00:19 bjoo Exp $

/*! @file
 * @brief Staggered-like fermion actions
 */

#ifndef __stagtype_fermact_s_h__
#define __stagtype_fermact_s_h__

#include "fermact.h"
#include "eo_linop.h"

namespace Chroma
{

  //-------------------------------------------------------------------------------------------
  //! Staggered-like fermion actions
  /*! @ingroup actions
   *
   * Staggered-like fermion actions
   */
  template<typename T, typename P, typename Q>
  class StaggeredTypeFermAct : public DiffFermAct4D<T,P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~StaggeredTypeFermAct() {}

    //! Return the quark mass
    virtual const Real getQuarkMass() const = 0;

    //! Return a linear operator solver for this action to solve M*psi=chi 
    /*! Default implementation provided */
    virtual LinOpSystemSolver<T>* invLinOp(Handle< FermState<T,P,Q> > state,
					   const GroupXML_t& invParam) const;

    //! Return a linear operator solver for this action to solve MdagM*psi=chi 
    /*! Default implementation provided */
    virtual MdagMSystemSolver<T>* invMdagM(Handle< FermState<T,P,Q> > state,
					   const GroupXML_t& invParam) const;

    //! Return a multi-shift linear operator solver for this action to solve (M+shift)*psi=chi 
    /*! Default implementation */
    virtual LinOpMultiSystemSolver<T>* mInvLinOp(Handle< FermState<T,P,Q> > state,
						 const GroupXML_t& invParam) const;

    //! Return a multi-shift linear operator solver for this action to solve (MdagM+shift)*psi=chi 
    /*! Default implementation */
    virtual MdagMMultiSystemSolver<T>* mInvMdagM(Handle< FermState<T,P,Q> > state,
						 const GroupXML_t& invParam) const;

    //! Return a multi-shift linear operator solver for this action to solve (MdagM+shift)*psi=chi 
    /*! Default implementation */
    virtual MdagMMultiSystemSolverAccumulate<T>* mInvMdagMAcc(Handle< FermState<T,P,Q> > state,
						 const GroupXML_t& invParam) const;

    //! Given a complete propagator as a source, this does all the inversions needed
    /*!
     * Provides a default version
     *
     * \param q_sol         quark propagator ( Write )
     * \param q_src         source ( Read )
     * \param xml_out       diagnostic output ( Modify )
     * \param state         gauge connection state ( Read )
     * \param invParam      inverter parameters ( Read )
     * \param quarkSpinType compute only a non-relativistic prop ( Read )
     * \param ncg_had       number of solver iterations ( Write )
     */
    virtual void quarkProp(typename PropTypeTraits<T>::Type_t& q_sol,
			   XMLWriter& xml_out,
			   const typename PropTypeTraits<T>::Type_t& q_src,
			   Handle< FermState<T,P,Q> > state,
			   const GroupXML_t& invParam,
			   QuarkSpinType quarkSpinType,
			   int& ncg_had) const;
  };


  //-------------------------------------------------------------------------------------------
  //! Staggered-like fermion actions with derivatives
  /*! @ingroup actions
   *
   * Staggered-like fermion actions
   */
  template<typename T, typename P, typename Q>
  class UnprecStaggeredTypeFermAct : public StaggeredTypeFermAct<T,P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~UnprecStaggeredTypeFermAct() {}

    //! Override to produce an even-odd prec. linear operator for this action
    /*! Covariant return rule - override base class function */
    virtual UnprecLinearOperator<T,P,Q>* linOp(Handle< FermState<T,P,Q> > state) const = 0;
  };


  //! Even-odd preconditioned Staggered-like fermion actions
  /*! @ingroup actions
   *
   * Even-odd preconditioned like Staggered-like fermion actions
   */
  template<typename T, typename P, typename Q>
  class EvenOddStaggeredTypeFermAct : public StaggeredTypeFermAct<T,P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~EvenOddStaggeredTypeFermAct() {}

    //! Override to produce an even-odd prec. linear operator for this action
    /*! Covariant return rule - override base class function */
    virtual EvenOddLinearOperator<T,P,Q>* linOp(Handle< FermState<T,P,Q> > state) const = 0;

    //! Return quark prop solver, solution of unpreconditioned system
    /*! Default implementation provided */
    virtual SystemSolver<T>* qprop(Handle< FermState<T,P,Q> > state,
				   const GroupXML_t& invParam) const;
  };

}


#endif
