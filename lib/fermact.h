// -*- C++ -*-
// $Id: fermact.h,v 1.15 2004-12-29 22:13:40 edwards Exp $

/*! @file
 * @brief Class structure for fermion actions
 */

#ifndef __fermact_h__
#define __fermact_h__

using namespace QDP;

#include "invtype.h"
#include "state.h"
#include "linearop.h"
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

    //! Compute quark propagator over simpler type
    /*! 
     * Solves  M.psi = chi
     *
     * \param psi      quark propagator ( Modify )
     * \param state    gauge connection state ( Read )
     * \param chi      source ( Read )
     * \param invParam inverter parameters ( Read (
     * \param ncg_had  number of CG iterations ( Write )
     *
     */
    virtual void qprop(T& psi, 
		       Handle<const ConnectState> state, 
		       const T& chi, 
		       const InvertParam_t& invParam,
		       int& ncg_had) const = 0;

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

    //! Compute quark propagator over simpler type
    /*! 
     * Solves  M.psi = chi
     *
     * \param psi      quark propagator ( Modify )
     * \param state    gauge connection state ( Read )
     * \param chi      source ( Read )
     * \param invParam inverter parameters ( Read (
     * \param ncg_had  number of CG iterations ( Write )
     *
     */
    virtual void qprop(T& psi, 
		       Handle<const ConnectState> state, 
		       const T& chi, 
		       const InvertParam_t& invParam,
		       int& ncg_had) const;
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
    virtual Real quark_mass() const = 0;

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

    //! Compute quark propagator over base type
    /*! 
     * Solves  M.psi = chi
     *
     * \param psi      quark propagator ( Modify )
     * \param state    gauge connection state ( Read )
     * \param chi      source ( Modify )
     * \param invParam inverter parameters ( Read (
     * \param ncg_had  number of CG iterations ( Write )
     *
     */
    virtual void qpropT(multi1d<T>& psi, 
			Handle<const ConnectState> state, 
			const multi1d<T>& chi, 
			const InvertParam_t& invParam,
			int& ncg_had) const;
  };


  //-------------------------------------------------------------------------------------------
  //! Wilson-like fermion actions
  /*! @ingroup actions
   *
   * Wilson-like fermion actions
   */
  template<typename T>
  class WilsonTypeFermAct : public FermAct4D<T>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~WilsonTypeFermAct() {}

    //! Produce a hermitian version of the linear operator
    virtual const LinearOperator<T>* gamma5HermLinOp(Handle<const ConnectState> state) const = 0;

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
  template<typename T>
  class WilsonTypeFermAct5D : public FermAct5D<T>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~WilsonTypeFermAct5D() {}

    //! Produce a hermitian version of the linear operator
    virtual const LinearOperator< multi1d<T> >* gamma5HermLinOp(Handle<const ConnectState> state) const = 0;

    //! Produce an unpreconditioned linear operator projecting 5D to 4D (the inverse of qprop below)
    virtual const LinearOperator<T>* linOp4D(Handle<const ConnectState> state,
					     const Real& m_q,
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
  //! Unpreconditioned Wilson-like fermion actions
  /*! @ingroup actions
   *
   * Unpreconditioned like Wilson-like fermion actions
   */
  template<typename T>
  class UnprecWilsonTypeFermActBase : public WilsonTypeFermAct<T>
  {
  public:
    //! Produce a linear operator for this action
    virtual const UnprecLinearOperatorBase<T>* linOp(Handle<const ConnectState> state) const = 0;

    //! Virtual destructor to help with cleanup;
    virtual ~UnprecWilsonTypeFermActBase() {}
  };



  //-------------------------------------------------------------------------------------------
  //! Even-odd preconditioned Wilson-like fermion actions
  /*! @ingroup actions
   *
   * Even-odd preconditioned like Wilson-like fermion actions
   */
  template<typename T>
  class EvenOddPrecWilsonTypeFermActBase : public WilsonTypeFermAct<T>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~EvenOddPrecWilsonTypeFermActBase() {}

    //! Override to produce an even-odd prec. linear operator for this action
    /*! Covariant return rule - override base class function */
    virtual const EvenOddPrecLinearOperatorBase<T>* linOp(Handle<const ConnectState> state) const = 0;

    //! Compute quark propagator over base type
    /*! 
     * Solves  M.psi = chi
     *
     * Provides a default version
     *
     * \param psi      quark propagator ( Modify )
     * \param state    gauge connection state ( Read )
     * \param chi      source ( Modify )
     * \param invParam inverter parameters ( Read (
     * \param ncg_had  number of CG iterations ( Write )
     */
    virtual void qprop(T& psi, 
		       Handle<const ConnectState> state, 
		       const T& chi, 
		       const InvertParam_t& invParam,
		       int& ncg_had) const;
  };


  //-------------------------------------------------------------------------------------------
  //! Unpreconditioned Wilson-like fermion actions over extra dim arrays
  /*! @ingroup actions
   *
   * Unpreconditioned like Wilson-like fermion actions
   * Here, use arrays of matter fields.
   */
  template<typename T>
  class UnprecWilsonTypeFermActBase5D : public WilsonTypeFermAct5D<T>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~UnprecWilsonTypeFermActBase5D() {}

    //! Produce a linear operator for this action
    virtual const UnprecLinearOperatorBase< multi1d<T> >* linOp(Handle<const ConnectState> state) const = 0;
  };



  //-------------------------------------------------------------------------------------------
  //! Even-odd preconditioned Wilson-like fermion actions over extra dim arrays
  /*! @ingroup actions
   *
   * Even-odd preconditioned like Wilson-like fermion actions.
   * Here, use arrays of matter fields.
   */
  template<typename T>
  class EvenOddPrecWilsonTypeFermActBase5D : public WilsonTypeFermAct5D<T>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~EvenOddPrecWilsonTypeFermActBase5D() {}

    //! Override to produce an even-odd prec. linear operator for this action
    /*! Covariant return rule - override base class function */
    virtual const EvenOddPrecLinearOperatorBase< multi1d<T> >* linOp(Handle<const ConnectState> state) const = 0;

    //! Compute quark propagator over base type
    /*! 
     * Solves  M.psi = chi
     *
     * Provides a default version
     *
     * \param psi      quark propagator ( Modify )
     * \param state    gauge connection state ( Read )
     * \param chi      source ( Modify )
     * \param invParam inverter parameters ( Read (
     * \param ncg_had  number of CG iterations ( Write )
     */
    virtual void qpropT(multi1d<T>& psi, 
			Handle<const ConnectState> state, 
			const multi1d<T>& chi, 
			const InvertParam_t& invParam,
			int& ncg_had) const;
  };


  //-------------------------------------------------------------------------------------------
  //! Unpreconditioned Wilson-like fermion actions with derivatives
  /*! @ingroup actions
   *
   * Unpreconditioned like Wilson-like fermion actions
   */
  template<typename T, typename P>
  class UnprecWilsonTypeFermAct : public UnprecWilsonTypeFermActBase<T>
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
  class EvenOddPrecWilsonTypeFermAct : public EvenOddPrecWilsonTypeFermActBase<T>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~EvenOddPrecWilsonTypeFermAct() {}

    //! Override to produce an even-odd prec. linear operator for this action
    /*! Covariant return rule - override base class function */
    virtual const EvenOddPrecLinearOperator<T,P>* linOp(Handle<const ConnectState> state) const = 0;
  };



  //-------------------------------------------------------------------------------------------
  //! Unpreconditioned Wilson-like fermion actions in extra dims with derivatives
  /*! @ingroup actions
   *
   * Unpreconditioned like Wilson-like fermion actions
   * Here, use arrays of matter fields.
   */
  template<typename T, typename P>
  class UnprecWilsonTypeFermAct5D : public UnprecWilsonTypeFermActBase5D<T>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~UnprecWilsonTypeFermAct5D() {}

    //! Produce a linear operator for this action
    virtual const UnprecLinearOperator< multi1d<T>, P>* linOp(Handle<const ConnectState> state) const = 0;
  };



  //! Even-odd preconditioned Wilson-like fermion actions including derivatives
  /*! @ingroup actions
   *
   * Even-odd preconditioned like Wilson-like fermion actions
   * Here, use arrays of matter fields.
   */
  template<typename T, typename P>
  class EvenOddPrecWilsonTypeFermAct5D : public EvenOddPrecWilsonTypeFermActBase5D<T>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~EvenOddPrecWilsonTypeFermAct5D() {}

    //! Override to produce an even-odd prec. linear operator for this action
    /*! Covariant return rule - override base class function */
    virtual const EvenOddPrecLinearOperator< multi1d<T>, P>* linOp(Handle<const ConnectState> state) const = 0;
  };



  //-------------------------------------------------------------------------------------------
  //! Staggered-like fermion actions
  /*! @ingroup actions
   *
   * Staggered-like fermion actions
   */
  template<typename T>
  class StaggeredTypeFermAct : public FermAct4D<T>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~StaggeredTypeFermAct() {}

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


  //! Staggered-like fermion actions
  /*! @ingroup actions
   *
   * Staggered-like fermion actions
   */
  template<typename T>
  class UnprecStaggeredTypeFermActBase : public StaggeredTypeFermAct<T>
  {
  public:
    //! Override to produce an even-odd prec. linear operator for this action
    /*! Covariant return rule - override base class function */
    virtual const UnprecLinearOperatorBase<T>* linOp(Handle<const ConnectState> state) const = 0;

    //! Virtual destructor to help with cleanup;
    virtual ~UnprecStaggeredTypeFermActBase() {}
  };


  //! Even-odd preconditioned Staggered-like fermion actions
  /*! @ingroup actions
   *
   * Even-odd preconditioned like Staggered-like fermion actions
   */
  template<typename T>
  class EvenOddStaggeredTypeFermActBase : public StaggeredTypeFermAct<T>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~EvenOddStaggeredTypeFermActBase() {}

    //! Override to produce an even-odd prec. linear operator for this action
    /*! Covariant return rule - override base class function */
    virtual const EvenOddLinearOperatorBase<T>* linOp(Handle<const ConnectState> state) const = 0;

    //! Provide a default version of qprop
    void qprop(T& psi,
	       Handle<const ConnectState> state,
	       const T& chi,
	       const InvertParam_t& invParam,
	       int& ncg_had);

    //! Return the quark mass
    virtual const Real  getQuarkMass() const = 0;
  };



  //-------------------------------------------------------------------------------------------
  //! Staggered-like fermion actions with derivatives
  /*! @ingroup actions
   *
   * Staggered-like fermion actions
   */
  template<typename T, typename P>
  class UnprecStaggeredTypeFermAct : public UnprecStaggeredTypeFermActBase<T>
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
  class EvenOddStaggeredTypeFermAct : public EvenOddStaggeredTypeFermActBase<T>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~EvenOddStaggeredTypeFermAct() {}

    //! Override to produce an even-odd prec. linear operator for this action
    /*! Covariant return rule - override base class function */
    virtual const EvenOddLinearOperator<T,P>* linOp(Handle<const ConnectState> state) const = 0;
  };

}

using namespace Chroma;

#endif
