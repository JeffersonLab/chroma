// -*- C++ -*-
// $Id: fermactarray.h,v 1.1 2003-11-12 22:09:55 edwards Exp $

/*! @file
 * @brief Class structures for bi-local matter actions - here in array form
 */

#ifndef __fermactarrayact_h__
#define __fermactarrayact_h__

using namespace QDP;

enum InvType {
  CG_INVERTER = 21, 
  MR_INVERTER = 22,
  BICG_INVERTER = 23,
  CR_INVERTER = 24};

#include "linearop.h"

//! Base class for quadratic matter actions (e.g., fermions)
/*! @ingroup actions
 *
 * Supports creation and application for quadratic actions
 */

template<typename T>
class FermionActionArray
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
  virtual void qprop(multi1d<T>& psi, 
		     const multi1d<LatticeColorMatrix>& u, 
		     const multi1d<T>& chi, 
		     enum InvType invType,
		     const Real& RsdCG, 
		     int MaxCG, int& ncg_had) const;

  //! Compute dS_f/dU
  /*! NOTE: maybe this should produce a derivative foundry class object */
  virtual multi1d<LatticeColorMatrix> dsdu(const multi1d<LatticeColorMatrix>& u,
					   const multi1d<T>& psi) const = 0;

  //! Virtual destructor to help with cleanup;
  virtual ~FermionActionArray() {}
};


//! Wilson-like fermion actions
/*! @ingroup actions
 *
 * Wilson-like fermion actions
 */
template<typename T>
class WilsonTypeFermActArray : public FermionActionArray<T>
{
public:
};


//! Unpreconditioned Wilson-like fermion actions
/*! @ingroup actions
 *
 * Unpreconditioned like Wilson-like fermion actions
 */
template<typename T>
class UnprecWilsonTypeFermActArray : public WilsonTypeFermActArray<T>
{
public:
};


//! Even-odd preconditioned Wilson-like fermion actions
/*! @ingroup actions
 *
 * Even-odd preconditioned like Wilson-like fermion actions
 */
template<typename T>
class EvenOddPrecWilsonTypeFermActArray : public WilsonTypeFermActArray<T>
{
public:
};



//! Staggered-like fermion actions
/*! @ingroup actions
 *
 * Staggered-like fermion actions
 */
template<typename T>
class StaggeredTypeFermActArray : public FermionActionArray<T>
{
public:
};


//! Unpreconditioned Staggered-like fermion actions
/*! @ingroup actions
 *
 * Unpreconditioned like Staggered-like fermion actions
 */
template<typename T>
class UnprecStaggeredTypeFermActArray : public StaggeredTypeFermActArray<T>
{
public:
};


//! Even-odd preconditioned Staggered-like fermion actions
/*! @ingroup actions
 *
 * Even-odd preconditioned like Staggered-like fermion actions
 */
template<typename T>
class EvenOddPrecStaggeredTypeFermActArray : public StaggeredTypeFermActArray<T>
{
public:
};


#endif
