// -*- C++ -*-
// $Id: fermact.h,v 1.12 2003-11-23 05:51:05 edwards Exp $

/*! @file
 * @brief Class structure for fermion actions
 */

#ifndef __fermact_h__
#define __fermact_h__

using namespace QDP;

#include "invtype.h"
#include "linearop.h"

//! Base class for quadratic matter actions (e.g., fermions)
/*! @ingroup actions
 *
 * Supports creation and application for quadratic actions
 */
template<typename T>
class FermionAction
{
public:
#if 0
  //! Produce a foundry class for linear operators
  /*! NOTE: maybe this should be abstracted to a foundry class */
//  virtual LinOpFoundry* linop(const multi1d<LatticeColorMatrix>& _u) const = 0;
#endif

  //! Produce a linear operator for this action
  /*! NOTE: maybe this should be abstracted to a foundry class object */
  virtual const LinearOperator<T>* linOp(const multi1d<LatticeColorMatrix>& u) const = 0;

  //! Produce a linear operator M^dag.M for this action
  /*! NOTE: maybe this should be abstracted to a foundry class object */
  virtual const LinearOperator<T>* lMdagM(const multi1d<LatticeColorMatrix>& u) const = 0;

  //! Compute quark propagator over base type
  /*! 
   * Solves  M.psi = chi
   *
   * \param psi      quark propagator ( Modify )
   * \param u        gauge field ( Read )
   * \param chi      source ( Modify )
   * \param invType  inverter type ( Read (
   * \param RsdCG    CG (or MR) residual used here ( Read )
   * \param MaxCG    maximum number of CG iterations ( Read )
   * \param ncg_had  number of CG iterations ( Write )
   *
   * NOTE: maybe this should produce a quark prop foundry class object 
   */
  virtual void qpropT(T& psi, 
		      const multi1d<LatticeColorMatrix>& u, 
		      const T& chi, 
		      enum InvType invType,
		      const Real& RsdCG, 
		      int MaxCG, int& ncg_had) const;

  //! Compute quark propagator over simpler type
  /*! 
   * Solves  M.psi = chi
   *
   * \param psi      quark propagator ( Modify )
   * \param u        gauge field ( Read )
   * \param chi      source ( Modify )
   * \param invType  inverter type ( Read (
   * \param RsdCG    CG (or MR) residual used here ( Read )
   * \param MaxCG    maximum number of CG iterations ( Read )
   * \param ncg_had  number of CG iterations ( Write )
   *
   * NOTE: maybe this should produce a quark prop foundry class object 
   */
  virtual void qprop(typename BaseType<T>::Type_t& psi, 
		     const multi1d<LatticeColorMatrix>& u, 
		     const typename BaseType<T>::Type_t& chi, 
		     enum InvType invType,
		     const Real& RsdCG, 
		     int MaxCG, int& ncg_had) const;
  
  //! Compute dS_f/dU
  /*! NOTE: maybe this should produce a derivative foundry class object */
  virtual void dsdu(multi1d<LatticeColorMatrix>& result,
		    const multi1d<LatticeColorMatrix>& u,
		    const T& psi) const
    {
      QDPIO::cerr << "FermionAction::dsdu not implemented" << endl;
      QDP_abort(1);
    }

  //! Virtual destructor to help with cleanup;
  virtual ~FermionAction() {}
};


//! Partial specialization for base class of quadratic matter actions (e.g., fermions)
/*! @ingroup actions
 *
 * Supports creation and application for quadratic actions, specialized
 * to support arrays of matter fields
 */
template<typename T>
class FermionAction< multi1d<T> >
{
public:
  //! Expected length of array index
  virtual int size() const = 0;

  //! Produce a linear operator for this action
  /*! NOTE: maybe this should be abstracted to a foundry class object */
  virtual const LinearOperator< multi1d<T> >* linOp(const multi1d<LatticeColorMatrix>& u) const = 0;

  //! Produce a linear operator M^dag.M for this action
  /*! NOTE: maybe this should be abstracted to a foundry class object */
  virtual const LinearOperator< multi1d<T> >* lMdagM(const multi1d<LatticeColorMatrix>& u) const = 0;

  //! Compute quark propagator over base type
  /*! 
   * Solves  M.psi = chi
   *
   * \param psi      quark propagator ( Modify )
   * \param u        gauge field ( Read )
   * \param chi      source ( Modify )
   * \param invType  inverter type ( Read (
   * \param RsdCG    CG (or MR) residual used here ( Read )
   * \param MaxCG    maximum number of CG iterations ( Read )
   * \param ncg_had  number of CG iterations ( Write )
   *
   * NOTE: maybe this should produce a quark prop foundry class object 
   */
  virtual void qpropT(multi1d<T>& psi, 
		      const multi1d<LatticeColorMatrix>& u, 
		      const multi1d<T>& chi, 
		      enum InvType invType,
		      const Real& RsdCG, 
		      int MaxCG, int& ncg_had) const;

  //! Compute quark propagator over base type
  /*! 
   * Solves  M.psi = chi
   *
   * \param psi      quark propagator ( Modify )
   * \param u        gauge field ( Read )
   * \param chi      source ( Modify )
   * \param invType  inverter type ( Read (
   * \param RsdCG    CG (or MR) residual used here ( Read )
   * \param MaxCG    maximum number of CG iterations ( Read )
   * \param ncg_had  number of CG iterations ( Write )
   *
   * NOTE: maybe this should produce a quark prop foundry class object 
   */
  virtual void qprop(T& psi, 
		     const multi1d<LatticeColorMatrix>& u, 
		     const T& chi, 
		     enum InvType invType,
		     const Real& RsdCG, 
		     int MaxCG, int& ncg_had) const = 0;

  //! Compute dS_f/dU
  /*! NOTE: maybe this should produce a derivative foundry class object */
  virtual void dsdu(multi1d<LatticeColorMatrix>& result,
		    const multi1d<LatticeColorMatrix>& u,
		    const multi1d<T>& psi) const
    {
      QDPIO::cerr << "FermionAction::dsdu not implemented" << endl;
      QDP_abort(1);
    }

  //! Virtual destructor to help with cleanup;
  virtual ~FermionAction() {}
};


//! Wilson-like fermion actions
/*! @ingroup actions
 *
 * Wilson-like fermion actions
 */
template<typename T>
class WilsonTypeFermAct : public FermionAction<T>
{
public:
};


//! Unpreconditioned Wilson-like fermion actions
/*! @ingroup actions
 *
 * Unpreconditioned like Wilson-like fermion actions
 */
template<typename T>
class UnprecWilsonTypeFermAct : public WilsonTypeFermAct<T>
{
public:
#if 0
  //! Compute quark propagator
  /*! NOTE: maybe this should produce a quark prop foundry */
  void qprop(T& psi, 
	     const multi1d<LatticeColorMatrix>& u, 
	     const T& chi, 
	     enum InvType invType,
	     const Real& RsdCG, 
	     int MaxCG, int& ncg_had) const;
#endif
};


//! Even-odd preconditioned Wilson-like fermion actions
/*! @ingroup actions
 *
 * Even-odd preconditioned like Wilson-like fermion actions
 */
template<typename T>
class EvenOddPrecWilsonTypeFermAct : public WilsonTypeFermAct<T>
{
public:
  //! Override to produce an even-odd prec. linear operator for this action
  /*! Covariant return rule - override base class function */
  virtual const EvenOddPrecLinearOperator<T>* linOp(const multi1d<LatticeColorMatrix>& u) const = 0;

  //! Compute quark propagator over base type
  /*! 
   * Solves  M.psi = chi
   *
   * Provides a default version
   *
   * \param psi      quark propagator ( Modify )
   * \param u        gauge field ( Read )
   * \param chi      source ( Modify )
   * \param invType  inverter type ( Read (
   * \param RsdCG    CG (or MR) residual used here ( Read )
   * \param MaxCG    maximum number of CG iterations ( Read )
   * \param ncg_had  number of CG iterations ( Write )
   */
  virtual void qpropT(T& psi, 
		      const multi1d<LatticeColorMatrix>& u, 
		      const T& chi, 
		      enum InvType invType,
		      const Real& RsdCG, 
		      int MaxCG, int& ncg_had) const;

  //! Compute quark propagator over simpler type
  /*! 
   * Solves  M.psi = chi
   *
   * Provides a default version
   *
   * \param psi      quark propagator ( Modify )
   * \param u        gauge field ( Read )
   * \param chi      source ( Modify )
   * \param invType  inverter type ( Read (
   * \param RsdCG    CG (or MR) residual used here ( Read )
   * \param MaxCG    maximum number of CG iterations ( Read )
   * \param ncg_had  number of CG iterations ( Write )
   */
  virtual void qprop(typename BaseType<T>::Type_t& psi, 
		     const multi1d<LatticeColorMatrix>& u, 
		     const typename BaseType<T>::Type_t& chi, 
		     enum InvType invType,
		     const Real& RsdCG, 
		     int MaxCG, int& ncg_had) const;
  
};


//! Partial specialization of even-odd preconditioned Wilson-like fermion actions over arrays
/*! @ingroup actions
 *
 * Even-odd preconditioned like Wilson-like fermion actions.
 * Here, use arrays of matter fields.
 */
template<typename T>
class EvenOddPrecWilsonTypeFermAct< multi1d<T> > : public WilsonTypeFermAct< multi1d<T> >
{
public:
  //! Override to produce an even-odd prec. linear operator for this action
  /*! Covariant return rule - override base class function */
  virtual const EvenOddPrecLinearOperator< multi1d<T> >* linOp(const multi1d<LatticeColorMatrix>& u) const = 0;

  //! Compute quark propagator over base type
  /*! 
   * Solves  M.psi = chi
   *
   * Provides a default version
   *
   * \param psi      quark propagator ( Modify )
   * \param u        gauge field ( Read )
   * \param chi      source ( Modify )
   * \param invType  inverter type ( Read (
   * \param RsdCG    CG (or MR) residual used here ( Read )
   * \param MaxCG    maximum number of CG iterations ( Read )
   * \param ncg_had  number of CG iterations ( Write )
   */
  virtual void qpropT(multi1d<T>& psi, 
		      const multi1d<LatticeColorMatrix>& u, 
		      const multi1d<T>& chi, 
		      enum InvType invType,
		      const Real& RsdCG, 
		      int MaxCG, int& ncg_had) const;

  //! Compute quark propagator over simpler type
  /*! 
   * Solves  M.psi = chi
   *
   * Provides a default version
   *
   * \param psi      quark propagator ( Modify )
   * \param u        gauge field ( Read )
   * \param chi      source ( Modify )
   * \param invType  inverter type ( Read (
   * \param RsdCG    CG (or MR) residual used here ( Read )
   * \param MaxCG    maximum number of CG iterations ( Read )
   * \param ncg_had  number of CG iterations ( Write )
   */
  virtual void qprop(typename BaseType<T>::Type_t& psi, 
		     const multi1d<LatticeColorMatrix>& u, 
		     const typename BaseType<T>::Type_t& chi, 
		     enum InvType invType,
		     const Real& RsdCG, 
		     int MaxCG, int& ncg_had) const;
  
};


//! Even-odd/unit-diag preconditioned Wilson-like fermion actions
/*! @ingroup actions
 *
 * Even-odd with a unit diagonal preconditioned like Wilson-like fermion actions
 */
template<typename T>
class DiagEvenOddPrecWilsonTypeFermAct : public EvenOddPrecWilsonTypeFermAct<T>
{
public:
#if 0
  //! Compute quark propagator
  /*! NOTE: maybe this should produce a quark prop foundry */
  void qprop(T& psi, 
	     const multi1d<LatticeColorMatrix>& u, 
	     const T& chi, 
	     enum InvType invType,
	     const Real& RsdCG, 
	     int MaxCG, int& ncg_had) const;
#endif
};


//! Staggered-like fermion actions
/*! @ingroup actions
 *
 * Staggered-like fermion actions
 */
template<typename T>
class StaggeredTypeFermAct : public FermionAction<T>
{
public:
};


//! Unpreconditioned Staggered-like fermion actions
/*! @ingroup actions
 *
 * Unpreconditioned like Staggered-like fermion actions
 */
template<typename T>
class UnprecStaggeredTypeFermAct : public StaggeredTypeFermAct<T>
{
public:
#if 0
  //! Compute quark propagator
  /*! NOTE: maybe this should produce a quark prop foundry */
  void qprop(T& psi, 
	     const multi1d<LatticeColorMatrix>& u, 
	     const T& chi, 
	     enum InvType invType,
	     const Real& RsdCG, 
	     int MaxCG, int& ncg_had) const;
#endif
};


//! Even-odd preconditioned Staggered-like fermion actions
/*! @ingroup actions
 *
 * Even-odd preconditioned like Staggered-like fermion actions
 */
template<typename T>
class EvenOddPrecStaggeredTypeFermAct : public StaggeredTypeFermAct<T>
{
public:
};


//! Even-odd/unit-diag preconditioned Staggered-like fermion actions
/*! @ingroup actions
 *
 * Even-odd with a unit diagonal preconditioned like Staggered-like fermion actions
 */
template<typename T>
class DiagEvenOddPrecStaggeredTypeFermAct : public EvenOddPrecStaggeredTypeFermAct<T>
{
public:
#if 0
  //! Compute quark propagator
  /*! NOTE: maybe this should produce a quark prop foundry */
  void qprop(T& psi, 
	     const multi1d<LatticeColorMatrix>& u, 
	     const T& chi, 
	     enum InvType invType,
	     const Real& RsdCG, 
	     int MaxCG, int& ncg_had) const;
#endif
};


#endif
