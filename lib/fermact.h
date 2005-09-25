// -*- C++ -*-
// $Id: fermact.h,v 2.0 2005-09-25 21:04:24 edwards Exp $

/*! @file
 * @brief Class structure for fermion actions
 */

#ifndef __fermact_h__
#define __fermact_h__

#include "chromabase.h"
#include "invtype.h"
#include "state.h"
#include "linearop.h"
#include "syssolver.h"
#include "fermbc.h"


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
   * The ConnectState holds gauge fields and whatever auxilliary info
   * is needed to create a specific dirac operator (linear operator)
   * on some background gauge field.
   *
   * The FermBC is the type of boundary conditions used for this action
   *
   * The linop and lmdagm functions create a linear operator on a 
   * fixed ConnectState
   *
   * The qprop solves for the full linear system undoing all preconditioning.
   * No direct functions are provided (or needed) to call a linear system
   * solver - that is a stand-alone function in the generic programming sense.
   *
   */
  template<typename T>
  class FermionAction
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~FermionAction() {}

    //! Given links, create the state needed for the linear operators
    virtual const ConnectState* createState(const multi1d<LatticeColorMatrix>& u) const = 0;

    //! Given the links create a state with additional info held by the XMLReader
    virtual const ConnectState* createState(const multi1d<LatticeColorMatrix>& u,
					    XMLReader& reader,
					    const string& path) const = 0;

    //! Return quark prop solver, solution of unpreconditioned system
    virtual const SystemSolver<T>* qprop(Handle<const ConnectState> state,
					 const InvertParam_t& invParam) const = 0;

    //! Given a complete propagator as a source, this does all the inversions needed
    /*!
     * \param q_sol    quark propagator ( Write )
     * \param q_src    source ( Read )
     * \param xml_out  diagnostic output ( Modify )
     * \param state    gauge connection state ( Read )
     * \param invParam inverter parameters ( Read )
     * \param nonRelProp compute only a non-relativistic prop ( Read )
     * \param ncg_had  number of solver iterations ( Write )
     */
    virtual void quarkProp(typename PropTypeTraits<T>::Type_t& q_sol,
			   XMLWriter& xml_out,
			   const typename PropTypeTraits<T>::Type_t& q_src,
			   Handle<const ConnectState> state,
			   const InvertParam_t& invParam,
			   bool nonRelProp,
			   int& ncg_had) = 0;

    //! Given a complete propagator as a source, this does all the inversions needed
    /*!
     * Provides a default version
     *
     * \param q_sol    quark propagator ( Write )
     * \param q_src    source ( Read )
     * \param xml_out  diagnostic output ( Modify )
     * \param t_src    time slice of source ( Read )
     * \param j_decay  direction of decay ( Read )
     * \param state    gauge connection state ( Read )
     * \param invParam inverter parameters ( Read )
     * \param nonRelProp compute only a non-relativistic prop ( Read )
     * \param obsvP    compute currents and residual mass ( Read )
     * \param ncg_had  number of solver iterations ( Write )
     */
    virtual void quarkProp(typename PropTypeTraits<T>::Type_t& q_sol,
			   XMLWriter& xml_out,
			   const typename PropTypeTraits<T>::Type_t& q_src,
			   int t_src, int j_decay,
			   Handle<const ConnectState> state,
			   const InvertParam_t& invParam,
			   bool nonRelProp,
			   bool obsvP,
			   int& ncg_had)
    {
      quarkProp(q_sol, xml_out, q_src, state, invParam, nonRelProp, ncg_had);
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
   * The ConnectState holds gauge fields and whatever auxilliary info
   * is needed to create a specific dirac operator (linear operator)
   * on some background gauge field.
   *
   * The FermBC is the type of boundary conditions used for this action
   *
   * The linop and lmdagm functions create a linear operator on a 
   * fixed ConnectState
   *
   * The qprop solves for the full linear system undoing all preconditioning.
   * No direct functions are provided (or needed) to call a linear system
   * solver - that is a stand-alone function in the generic programming sense.
   *
   */
  template<typename T>
  class FermAct4D : public FermionAction<T>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~FermAct4D() {}

    //! Given links, create the state needed for the linear operators
    /*! Default version uses a SimpleConnectState */
    virtual const ConnectState* createState(const multi1d<LatticeColorMatrix>& u) const
    {
      multi1d<LatticeColorMatrix> u_tmp = u;
      getFermBC().modifyU(u_tmp);
      return new SimpleConnectState(u_tmp);
    }

    //! Given the links create a state with additional info held by the XMLReader
    /*! This is a default version, which ignores the reader and just calls the default
      createState() function */
    virtual const ConnectState* createState(const multi1d<LatticeColorMatrix>& u,
					    XMLReader& reader,
					    const string& path) const
    {
      return createState(u);
    }

    //! Return the fermion BC object for this action
    /*! The user will supply the FermBC in a derived class */
    virtual const FermBC<T>& getFermBC() const = 0;

    //! Produce a linear operator for this action
    virtual const LinearOperator<T>* linOp(Handle<const ConnectState> state) const = 0;

    //! Produce a linear operator M^dag.M for this action
    virtual const LinearOperator<T>* lMdagM(Handle<const ConnectState> state) const = 0;

    //! Return quark prop solver, solution of unpreconditioned system
    /*! Default implementation provided */
    virtual const SystemSolver<T>* qprop(Handle<const ConnectState> state,
					 const InvertParam_t& invParam) const;
  };


  //-------------------------------------------------------------------------------------------
  //! Base class for quadratic matter actions (e.g., fermions)
  /*! @ingroup actions
   *
   * Supports creation and application for quadratic actions with derivatives
   * This is basically a foundry class with additional operations.
   */
  template<typename T, typename P>
  class DiffFermAct4D : public FermAct4D<T>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~DiffFermAct4D() {}

    //! Produce a linear operator for this action
    virtual const DiffLinearOperator<T,P>* linOp(Handle<const ConnectState> state) const = 0;
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
   * The ConnectState holds gauge fields and whatever auxilliary info
   * is needed to create a specific dirac operator (linear operator)
   * on some background gauge field.
   *
   * The FermBC is the type of boundary conditions used for this action
   *
   * The linop and lmdagm functions create a linear operator on a 
   * fixed ConnectState
   *
   * The qprop solves for the full linear system undoing all preconditioning.
   * No direct functions are provided (or needed) to call a linear system
   * solver - that is a stand-alone function in the generic programming sense.
   *
   */
  template<typename T>
  class FermAct5D : public FermionAction<T>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~FermAct5D() {}

    //! Return the quark mass
    virtual Real getQuarkMass() const = 0;

    //! Expected length of array index
    virtual int size() const = 0;

    //! Given links, create the state needed for the linear operators
    /*! Default version uses a SimpleConnectState */
    virtual const ConnectState* createState(const multi1d<LatticeColorMatrix>& u) const
    {
      multi1d<LatticeColorMatrix> u_tmp = u;
      getFermBC().modifyU(u_tmp);
      return new SimpleConnectState(u_tmp);
    }

    //! Given the links create a state with additional info held by the XMLReader
    /*! This is a default version, which ignores the reader and just calls the default
      createState() function */
    virtual const ConnectState* createState(const multi1d<LatticeColorMatrix>& u,
					    XMLReader& reader,
					    const string& path) const
    {
      return createState(u);
    }

    //! Return the fermion BC object for this action
    /*! The user will supply the FermBC in a derived class */
    virtual const FermBC< multi1d<T> >& getFermBC() const = 0;

    //! Produce a linear operator for this action
    virtual const LinearOperator< multi1d<T> >* linOp(Handle<const ConnectState> state) const = 0;

    //! Produce a linear operator M^dag.M for this action
    virtual const LinearOperator< multi1d<T> >* lMdagM(Handle<const ConnectState> state) const = 0;

    //! Produce a Pauli-Villars linear operator for this action
    virtual const LinearOperator< multi1d<T> >* linOpPV(Handle<const ConnectState> state) const = 0;

    //! Return quark prop solver, solution of unpreconditioned system
    /*! Default implementation provided */
    virtual const SystemSolver< multi1d<T> >* qpropT(Handle<const ConnectState> state,
						     const InvertParam_t& invParam) const;
  };


  //-------------------------------------------------------------------------------------------
  //! Base class of quadratic matter actions (e.g., fermions) living in an extra dimension
  /*! @ingroup actions
   *
   * Supports creation and application for quadratic actions, specialized
   * to support arrays of matter fields including derivatives
   */
  template<typename T, typename P>
  class DiffFermAct5D : public FermAct5D<T>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~DiffFermAct5D() {}

    //! Produce a linear operator for this action
    virtual const DiffLinearOperator<multi1d<T>, P>* linOp(Handle<const ConnectState> state) const = 0;

    //! Produce a Pauli-Villars linear operator for this action
    virtual const DiffLinearOperator<multi1d<T>, P>* linOpPV(Handle<const ConnectState> state) const = 0;
  };


  //-------------------------------------------------------------------------------------------
  //! Wilson-like fermion actions
  /*! @ingroup actions
   *
   * Wilson-like fermion actions
   */
  template<typename T, typename P>
  class WilsonTypeFermAct : public DiffFermAct4D<T,P>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~WilsonTypeFermAct() {}

    //! Produce a hermitian version of the linear operator
    virtual const LinearOperator<T>* hermitianLinOp(Handle<const ConnectState> state) const = 0;

    //! Given a complete propagator as a source, this does all the inversions needed
    /*!
     * Provides a default version
     *
     * \param q_sol    quark propagator ( Write )
     * \param q_src    source ( Read )
     * \param xml_out  diagnostic output ( Modify )
     * \param state    gauge connection state ( Read )
     * \param invParam inverter parameters ( Read )
     * \param nonRelProp compute only a non-relativistic prop ( Read )
     * \param ncg_had  number of solver iterations ( Write )
     */
    virtual void quarkProp(typename PropTypeTraits<T>::Type_t& q_sol,
			   XMLWriter& xml_out,
			   const typename PropTypeTraits<T>::Type_t& q_src,
			   Handle<const ConnectState> state,
			   const InvertParam_t& invParam,
			   bool nonRelProp,
			   int& ncg_had);
  };


  //-------------------------------------------------------------------------------------------
  //! Wilson-like fermion actions
  /*! @ingroup actions
   *
   * Wilson-like fermion actions
   */
  template<typename T, typename P>
  class WilsonTypeFermAct5D : public DiffFermAct5D<T,P>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~WilsonTypeFermAct5D() {}

    //! Produce a hermitian version of the linear operator
    virtual const LinearOperator< multi1d<T> >* hermitianLinOp(Handle<const ConnectState> state) const = 0;

    //! Produce an unpreconditioned linear operator projecting 5D to 4D (the inverse of qprop below)
    virtual const LinearOperator<T>* linOp4D(Handle<const ConnectState> state,
					     const Real& m_q,
					     const InvertParam_t& invParam) const = 0;

    //! Produce a  DeltaLs = 1-epsilon^2(H) operator
    virtual const LinearOperator<T>* DeltaLs(Handle< const ConnectState> state,
					     const InvertParam_t& invParam) const = 0;

    //! Given a complete propagator as a source, this does all the inversions needed
    /*!
     * Provides a default version
     *
     * \param q_sol    quark propagator ( Write )
     * \param q_src    source ( Read )
     * \param xml_out  diagnostic output ( Modify )
     * \param state    gauge connection state ( Read )
     * \param invParam inverter parameters ( Read )
     * \param nonRelProp compute only a non-relativistic prop ( Read )
     * \param ncg_had  number of solver iterations ( Write )
     */
    virtual void quarkProp(typename PropTypeTraits<T>::Type_t& q_sol,
			   XMLWriter& xml_out,
			   const typename PropTypeTraits<T>::Type_t& q_src,
			   Handle<const ConnectState> state,
			   const InvertParam_t& invParam,
			   bool nonRelProp,
			   int& ncg_had);
  };


  //-------------------------------------------------------------------------------------------
  //! Unpreconditioned Wilson-like fermion actions with derivatives
  /*! @ingroup actions
   *
   * Unpreconditioned like Wilson-like fermion actions
   */
  template<typename T, typename P>
  class UnprecWilsonTypeFermAct : public WilsonTypeFermAct<T,P>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~UnprecWilsonTypeFermAct() {}

    //! Produce a linear operator for this action
    virtual const UnprecLinearOperator<T,P>* linOp(Handle<const ConnectState> state) const = 0;
  };



  //! Even-odd preconditioned Wilson-like fermion actions including derivatives
  /*! @ingroup actions
   *
   * Even-odd preconditioned like Wilson-like fermion actions
   */
  template<typename T, typename P>
  class EvenOddPrecWilsonTypeFermAct : public WilsonTypeFermAct<T,P>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~EvenOddPrecWilsonTypeFermAct() {}

    //! Override to produce an even-odd prec. linear operator for this action
    /*! Covariant return rule - override base class function */
    virtual const EvenOddPrecLinearOperator<T,P>* linOp(Handle<const ConnectState> state) const = 0;

    //! Return quark prop solver, solution of unpreconditioned system
    /*! Default implementation provided */
    virtual const SystemSolver<T>* qprop(Handle<const ConnectState> state,
					 const InvertParam_t& invParam) const;
  };



  //-------------------------------------------------------------------------------------------
  //! Unpreconditioned Wilson-like fermion actions in extra dims with derivatives
  /*! @ingroup actions
   *
   * Unpreconditioned like Wilson-like fermion actions
   * Here, use arrays of matter fields.
   */
  template<typename T, typename P>
  class UnprecWilsonTypeFermAct5D : public WilsonTypeFermAct5D<T,P>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~UnprecWilsonTypeFermAct5D() {}

    //! Produce a linear operator for this action
    virtual const UnprecLinearOperator<multi1d<T>, P>* linOp(Handle<const ConnectState> state) const = 0;

    //! Produce a Pauli-Villars linear operator for this action
    virtual const UnprecLinearOperator<multi1d<T>, P>* linOpPV(Handle<const ConnectState> state) const = 0;
  };



  //! Even-odd preconditioned Wilson-like fermion actions including derivatives
  /*! @ingroup actions
   *
   * Even-odd preconditioned like Wilson-like fermion actions
   * Here, use arrays of matter fields.
   */
  template<typename T, typename P>
  class EvenOddPrecWilsonTypeFermAct5D : public WilsonTypeFermAct5D<T,P>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~EvenOddPrecWilsonTypeFermAct5D() {}

    //! Override to produce an even-odd prec. linear operator for this action
    /*! Covariant return rule - override base class function */
    virtual const EvenOddPrecLinearOperator<multi1d<T>, P>* linOp(Handle<const ConnectState> state) const = 0;

    //! Override to produce an even-odd prec. Pauli-Villars linear operator for this action
    /*! Covariant return rule - override base class function */
    virtual const EvenOddPrecLinearOperator<multi1d<T>, P>* linOpPV(Handle<const ConnectState> state) const = 0;

    //! Return quark prop solver, solution of unpreconditioned system
    /*! Default implementation provided */
    virtual const SystemSolver< multi1d<T> >* qpropT(Handle<const ConnectState> state,
						     const InvertParam_t& invParam) const;
  };



  //-------------------------------------------------------------------------------------------
  //! Staggered-like fermion actions
  /*! @ingroup actions
   *
   * Staggered-like fermion actions
   */
  template<typename T, typename P>
  class StaggeredTypeFermAct : public DiffFermAct4D<T,P>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~StaggeredTypeFermAct() {}

    //! Return the quark mass
    virtual const Real getQuarkMass() const = 0;

    //! Given a complete propagator as a source, this does all the inversions needed
    /*!
     * Provides a default version
     *
     * \param q_sol    quark propagator ( Write )
     * \param q_src    source ( Read )
     * \param xml_out  diagnostic output ( Modify )
     * \param state    gauge connection state ( Read )
     * \param invParam inverter parameters ( Read )
     * \param nonRelProp compute only a non-relativistic prop ( Read )
     * \param ncg_had  number of solver iterations ( Write )
     */
    virtual void quarkProp(typename PropTypeTraits<T>::Type_t& q_sol,
			   XMLWriter& xml_out,
			   const typename PropTypeTraits<T>::Type_t& q_src,
			   Handle<const ConnectState> state,
			   const InvertParam_t& invParam,
			   bool nonRelProp,
			   int& ncg_had);
  };


  //-------------------------------------------------------------------------------------------
  //! Staggered-like fermion actions with derivatives
  /*! @ingroup actions
   *
   * Staggered-like fermion actions
   */
  template<typename T, typename P>
  class UnprecStaggeredTypeFermAct : public StaggeredTypeFermAct<T,P>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~UnprecStaggeredTypeFermAct() {}

    //! Override to produce an even-odd prec. linear operator for this action
    /*! Covariant return rule - override base class function */
    virtual const UnprecLinearOperator<T,P>* linOp(Handle<const ConnectState> state) const = 0;
  };


  //! Even-odd preconditioned Staggered-like fermion actions
  /*! @ingroup actions
   *
   * Even-odd preconditioned like Staggered-like fermion actions
   */
  template<typename T, typename P>
  class EvenOddStaggeredTypeFermAct : public StaggeredTypeFermAct<T,P>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~EvenOddStaggeredTypeFermAct() {}

    //! Override to produce an even-odd prec. linear operator for this action
    /*! Covariant return rule - override base class function */
    virtual const EvenOddLinearOperator<T,P>* linOp(Handle<const ConnectState> state) const = 0;

    //! Return quark prop solver, solution of unpreconditioned system
    /*! Default implementation provided */
    virtual const SystemSolver<T>* qprop(Handle<const ConnectState> state,
					 const InvertParam_t& invParam) const;
  };

}


#endif
