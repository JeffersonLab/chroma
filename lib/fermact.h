// -*- C++ -*-
// $Id: fermact.h,v 3.6 2008-09-08 16:00:19 bjoo Exp $

/*! @file
 * @brief Class structure for fermion actions
 */

#ifndef __fermact_h__
#define __fermact_h__

#include "chromabase.h"
#include "fermbc.h"
#include "state.h"
#include "create_state.h"
#include "linearop.h"
#include "lmdagm.h"
#include "syssolver.h"
#include "io/xml_group_reader.h"
#include "io/enum_io/enum_quarkspintype_io.h"
#include "actions/ferm/invert/syssolver_linop.h"
#include "actions/ferm/invert/syssolver_mdagm.h"
#include "actions/ferm/invert/multi_syssolver_linop.h"
#include "actions/ferm/invert/multi_syssolver_mdagm.h"
#include "actions/ferm/invert/multi_syssolver_mdagm_accumulate.h"
namespace Chroma
{
  //! Base class for quadratic matter actions (e.g., fermions)
  /*! @ingroup actions
   *
   * Supports creation and application for quadratic actions.
   * This is basically a foundry class with additional operations.
   *
   * The class holds info on the particulars of a bi-local action,
   * but it DOES NOT hold gauge fields. A specific dirac-operator
   * is a functional of the gauge fields, hence when a dirac-operator
   * is needed, it is created.
   *
   * The FermState holds gauge fields and whatever auxilliary info
   * is needed to create a specific dirac operator (linear operator)
   * on some background gauge field.
   *
   * The FermBC is the type of boundary conditions used for this action
   *
   * The linop and lmdagm functions create a linear operator on a 
   * fixed FermState
   *
   * The qprop solves for the full linear system undoing all preconditioning.
   * No direct functions are provided (or needed) to call a linear system
   * solver - that is a stand-alone function in the generic programming sense.
   *
   */
  template<typename T, typename P, typename Q>
  class FermionAction
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~FermionAction() {}

    //! Given links (coordinates Q) create the state needed for the linear operators
    virtual FermState<T,P,Q>* createState(const Q& q) const
    {
      return getCreateState()(q);
    }

    //! Given links (coordinates Q) create a state with additional info held by the XMLReader
    /*! 
     * THIS FUNCTION SHOULD DISAPPEAR.
     * This function is only really used by the overlap operators to hand
     * around the StateInfo stuff. We are moving to have this stuff live within
     * the FermState
     */
    virtual FermState<T,P,Q>* createState(const Q& q,
					  XMLReader& reader,
					  const string& path) const
    {
      return createState(q);
    }

    //! Return the fermion BC object for this action
    virtual const FermBC<T,P,Q>& getFermBC() const
    {
      return getCreateState().getBC();
    }

    //! Return the factory object that produces a state
    /*! 
     * The user will supply the FermState in a derived class 
     *
     * This method is public for the simple reason that
     * a fermact might be nested, in which case the container would
     * forward this function call to the children.
     */
    virtual const CreateFermState<T,P,Q>& getCreateState() const = 0;

    //! Return quark prop solver, solution of unpreconditioned system
    virtual SystemSolver<T>* qprop(Handle< FermState<T,P,Q> > state,
				   const GroupXML_t& invParam) const = 0;

    //! Given a complete propagator as a source, this does all the inversions needed
    /*!
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
			   int& ncg_had) const = 0;

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
     * \param t_src         time slice of source ( Read )
     * \param j_decay       direction of decay ( Read )
     * \param obsvP         compute currents and residual mass ( Read )
     * \param ncg_had       number of solver iterations ( Write )
     */
    virtual void quarkProp(typename PropTypeTraits<T>::Type_t& q_sol,
			   XMLWriter& xml_out,
			   const typename PropTypeTraits<T>::Type_t& q_src,
			   int t_src, int j_decay,
			   Handle< FermState<T,P,Q> > state,
			   const GroupXML_t& invParam,
			   QuarkSpinType quarkSpinType,
			   bool obsvP,
			   int& ncg_had) const
      {
	quarkProp(q_sol, xml_out, q_src, state, invParam, quarkSpinType, ncg_had);
      }

  };


  //-------------------------------------------------------------------------------------------
  //! Base class for quadratic matter actions (e.g., fermions)
  /*! @ingroup actions
   *
   * Supports creation and application for quadratic actions.
   * This is basically a foundry class with additional operations.
   *
   * The class holds info on the particulars of a bi-local action,
   * but it DOES NOT hold gauge fields. A specific dirac-operator
   * is a functional of the gauge fields, hence when a dirac-operator
   * is needed, it is created.
   *
   * The FermState<T,P,Q> holds gauge fields and whatever auxilliary info
   * is needed to create a specific dirac operator (linear operator)
   * on some background gauge field.
   *
   * The FermBC is the type of boundary conditions used for this action
   *
   * The linop and lmdagm functions create a linear operator on a 
   * fixed FermState
   *
   * The qprop solves for the full linear system undoing all preconditioning.
   * No direct functions are provided (or needed) to call a linear system
   * solver - that is a stand-alone function in the generic programming sense.
   *
   */
  template<typename T, typename P, typename Q>
  class FermAct4D : public FermionAction<T,P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~FermAct4D() {}

    //! Produce a linear operator for this action
    virtual LinearOperator<T>* linOp(Handle< FermState<T,P,Q> > state) const = 0;

    //! Produce a linear operator M^dag.M for this action
    virtual LinearOperator<T>* lMdagM(Handle< FermState<T,P,Q> > state) const = 0;

    //! Return a linear operator solver for this action to solve M*psi=chi 
    virtual LinOpSystemSolver<T>* invLinOp(Handle< FermState<T,P,Q> > state,
					   const GroupXML_t& invParam) const = 0;

    //! Return a linear operator solver for this action to solve MdagM*psi=chi 
    virtual MdagMSystemSolver<T>* invMdagM(Handle< FermState<T,P,Q> > state,
					   const GroupXML_t& invParam) const = 0;

    //! Return a multi-shift linear operator solver for this action to solve (M+shift)*psi=chi 
    virtual LinOpMultiSystemSolver<T>* mInvLinOp(Handle< FermState<T,P,Q> > state,
						 const GroupXML_t& invParam) const = 0;

    //! Return a multi-shift linear operator solver for this action to solve (MdagM+shift)*psi=chi 
    virtual MdagMMultiSystemSolver<T>* mInvMdagM(Handle< FermState<T,P,Q> > state,
						 const GroupXML_t& invParam) const = 0;

    //! Return a multi-shift accumulate solver.
    //  for solving  A \sum_j p_j (MdagM + shift_j)^{-1} psi_j
    virtual MdagMMultiSystemSolverAccumulate<T>* mInvMdagMAcc(Handle< FermState<T,P,Q> > state,
						 const GroupXML_t& invParam) const = 0;

    //! Return quark prop solver, solution of unpreconditioned system
    /*! Default implementation provided */
    virtual SystemSolver<T>* qprop(Handle< FermState<T,P,Q> > state,
				   const GroupXML_t& invParam) const;
  };


  //-------------------------------------------------------------------------------------------
  //! Base class for quadratic matter actions (e.g., fermions)
  /*! @ingroup actions
   *
   * Supports creation and application for quadratic actions with derivatives
   * This is basically a foundry class with additional operations.
   */
  template<typename T, typename P, typename Q>
  class DiffFermAct4D : public FermAct4D<T,P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~DiffFermAct4D() {}

    //! Produce a linear operator for this action
    virtual DiffLinearOperator<T,Q,P>* linOp(Handle< FermState<T,P,Q> > state) const = 0;

    //! Produce a linear operator M^dag.M for this action
    virtual DiffLinearOperator<T,Q,P>* lMdagM(Handle< FermState<T,P,Q> > state) const = 0;
  };


  //-------------------------------------------------------------------------------------------
  //! Base class of quadratic matter actions (e.g., fermions) living in an extra dimension
  /*! @ingroup actions
   *
   * Supports creation and application for quadratic actions, specialized
   * to support arrays of matter fields
   *
   * The class holds info on the particulars of a bi-local action,
   * but it DOES NOT hold gauge fields. A specific dirac-operator
   * is a functional of the gauge fields, hence when a dirac-operator
   * is needed, it is created.
   *
   * The FermState holds gauge fields and whatever auxilliary info
   * is needed to create a specific dirac operator (linear operator)
   * on some background gauge field.
   *
   * The FermBC is the type of boundary conditions used for this action
   *
   * The linop and lmdagm functions create a linear operator on a 
   * fixed FermState
   *
   * The qprop solves for the full linear system undoing all preconditioning.
   * No direct functions are provided (or needed) to call a linear system
   * solver - that is a stand-alone function in the generic programming sense.
   *
   */
  template<typename T, typename P, typename Q>
  class FermAct5D : public FermionAction<T,P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~FermAct5D() {}

    //! Return the quark mass
    virtual Real getQuarkMass() const = 0;

    //! Expected length of array index
    virtual int size() const = 0;

    //! Produce a linear operator for this action
    virtual LinearOperatorArray<T>* linOp(Handle< FermState<T,P,Q> > state) const = 0;

    //! Produce a linear operator M^dag.M for this action
    virtual LinearOperatorArray<T >* lMdagM(Handle< FermState<T,P,Q> > state) const = 0;

    //! Produce a Pauli-Villars linear operator for this action
    virtual LinearOperatorArray<T >* linOpPV(Handle< FermState<T,P,Q> > state) const = 0;

    //! Return a linear operator solver for this action to solve M*psi=chi 
    virtual LinOpSystemSolverArray<T>* invLinOp(Handle< FermState<T,P,Q> > state,
						const GroupXML_t& invParam) const = 0;

    //! Return a linear operator solver for this action to solve MdagM*psi=chi 
    virtual MdagMSystemSolverArray<T>* invMdagM(Handle< FermState<T,P,Q> > state,
						const GroupXML_t& invParam) const = 0;

    //! Return a linear operator solver for this action to solve PV*psi=chi 
    /*! Do we need this critter? */
    virtual LinOpSystemSolverArray<T>* invLinOpPV(Handle< FermState<T,P,Q> > state,
						  const GroupXML_t& invParam) const = 0;

    //! Return a linear operator solver for this action to solve PV^dag*PV*psi=chi 
    /*! This is used in two-flavor molecdyn monomials */
    virtual MdagMSystemSolverArray<T>* invMdagMPV(Handle< FermState<T,P,Q> > state,
						  const GroupXML_t& invParam) const = 0;
 
    //! Return a multi-shift linear operator solver for this action to solve (MdagM+shift)*psi=chi 
    virtual MdagMMultiSystemSolverArray<T>* mInvMdagM(Handle< FermState<T,P,Q> > state,
						      const GroupXML_t& invParam) const = 0;



    
    //! Return a multi-shift accumulate solver.
    //  for solving  A \sum_j p_j (MdagM + shift_j)^{-1} psi_j
    virtual MdagMMultiSystemSolverAccumulateArray<T>* mInvMdagMAcc(Handle< FermState<T,P,Q> > state,
								   const GroupXML_t& invParam) const = 0;


    //! Return a multi-shift linear operator solver for this action to solve (PV^dag*PV+shift)*psi=chi 
    virtual MdagMMultiSystemSolverArray<T>* mInvMdagMPV(Handle< FermState<T,P,Q> > state,
							const GroupXML_t& invParam) const = 0;


    //! Return a multi-shift accumulate solver.
    //  for solving  A \sum_j p_j (MdagM + shift_j)^{-1} psi_j
    //  but for the PV.
    virtual MdagMMultiSystemSolverAccumulateArray<T>* mInvMdagMPVAcc(Handle< FermState<T,P,Q> > state,
							const GroupXML_t& invParam) const = 0;


    //! Return quark prop solver, solution of unpreconditioned system
    /*! Default implementation provided */
    virtual SystemSolverArray<T>* qpropT(Handle< FermState<T,P,Q> > state,
					 const GroupXML_t& invParam) const;
  };


  //-------------------------------------------------------------------------------------------
  //! Base class of quadratic matter actions (e.g., fermions) living in an extra dimension
  /*! @ingroup actions
   *
   * Supports creation and application for quadratic actions, specialized
   * to support arrays of matter fields including derivatives
   */
  template<typename T, typename P, typename Q>
  class DiffFermAct5D : public FermAct5D<T,P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~DiffFermAct5D() {}

    //! Produce a linear operator for this action
    virtual DiffLinearOperatorArray<T,P,Q>* linOp(Handle< FermState<T,P,Q> > state) const = 0;

    //! Produce a linear operator M^dag.M for this action
    virtual DiffLinearOperatorArray<T,P,Q>* lMdagM(Handle< FermState<T,P,Q> > state) const = 0;

    //! Produce a Pauli-Villars linear operator for this action
    virtual DiffLinearOperatorArray<T,P,Q>* linOpPV(Handle< FermState<T,P,Q> > state) const = 0;
  };

}


#endif
