// -*- C++ -*-
// $Id: fermact.h,v 1.5 2003-10-10 03:46:46 edwards Exp $

/*! @file
 * @brief Class structure for fermion actions
 */

#ifndef __fermact_h__
#define __fermact_h__

using namespace QDP;

#include "linearop.h"

//! Base class for fermion actions
/*! @ingroup actions
 *
 * Supports creation and application for fermion actions
 */

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
  virtual const LinearOperator* linOp(const multi1d<LatticeColorMatrix>& _u) const = 0;

  //! Produce a linear operator M^dag.M for this action
  /*! NOTE: maybe this should be abstracted to a foundry class object */
  virtual const LinearOperator* lMdagM(const multi1d<LatticeColorMatrix>& _u) const = 0;

  //! Compute quark propagator
  /*! 
   * Solves  M.psi = chi
   *
   * \param psi      quark propagator ( Modify )
   * \param u        gauge field ( Read )
   * \param chi      source ( Modify )
   * \param RsdCG    CG (or MR) residual used here ( Read )
   * \param MaxCG    maximum number of CG iterations ( Read )
   * \param ncg_had  number of CG iterations ( Write )
   *
   * NOTE: maybe this should produce a quark prop foundry class object 
   */
  virtual void qprop(LatticeFermion& psi, 
		     const multi1d<LatticeColorMatrix>& u, 
		     const LatticeFermion& chi, 
		     const Real& RsdCG, 
		     int MaxCG, int& ncg_had) const = 0;

  //! Compute dS_f/dU
  /*! NOTE: maybe this should produce a derivative foundry class object */
  virtual multi1d<LatticeColorMatrix> dsdu(const multi1d<LatticeColorMatrix>& u,
					   const LatticeFermion& psi) const = 0;

  //! Virtual destructor to help with cleanup;
  virtual ~FermionAction() {}
};


//! Wilson-like fermion actions
/*! @ingroup actions
 *
 * Wilson-like fermion actions
 */
class WilsonTypeFermAct : public FermionAction
{
public:
};


//! Unpreconditioned Wilson-like fermion actions
/*! @ingroup actions
 *
 * Unpreconditioned like Wilson-like fermion actions
 */
class UnprecWilsonTypeFermAct : public WilsonTypeFermAct
{
public:
  //! Compute quark propagator
  /*! NOTE: maybe this should produce a quark prop foundry */
  void qprop(LatticeFermion& psi, 
	     const multi1d<LatticeColorMatrix>& u, 
	     const LatticeFermion& chi, 
	     const Real& RsdCG, 
	     int MaxCG, int& ncg_had) const;
};


//! Even-odd preconditioned Wilson-like fermion actions
/*! @ingroup actions
 *
 * Even-odd preconditioned like Wilson-like fermion actions
 */
class EvenOddPrecWilsonTypeFermAct : public WilsonTypeFermAct
{
public:
};


//! Even-odd/unit-diag preconditioned Wilson-like fermion actions
/*! @ingroup actions
 *
 * Even-odd with a unit diagonal preconditioned like Wilson-like fermion actions
 */
class DiagEvenOddPrecWilsonTypeFermAct : public EvenOddPrecWilsonTypeFermAct
{
public:
  //! Compute quark propagator
  /*! NOTE: maybe this should produce a quark prop foundry */
  void qprop(LatticeFermion& psi, 
	     const multi1d<LatticeColorMatrix>& u, 
	     const LatticeFermion& chi, 
	     const Real& RsdCG, 
	     int MaxCG, int& ncg_had) const;
};


//! Staggered-like fermion actions
/*! @ingroup actions
 *
 * Staggered-like fermion actions
 */
class StaggeredTypeFermAct : public FermionAction
{
public:
};


//! Unpreconditioned Staggered-like fermion actions
/*! @ingroup actions
 *
 * Unpreconditioned like Staggered-like fermion actions
 */
class UnprecStaggeredTypeFermAct : public StaggeredTypeFermAct
{
public:
  //! Compute quark propagator
  /*! NOTE: maybe this should produce a quark prop foundry */
  void qprop(LatticeFermion& psi, 
	     const multi1d<LatticeColorMatrix>& u, 
	     const LatticeFermion& chi, 
	     const Real& RsdCG, 
	     int MaxCG, int& ncg_had) const;
};


//! Even-odd preconditioned Staggered-like fermion actions
/*! @ingroup actions
 *
 * Even-odd preconditioned like Staggered-like fermion actions
 */
class EvenOddPrecStaggeredTypeFermAct : public StaggeredTypeFermAct
{
public:
};


//! Even-odd/unit-diag preconditioned Staggered-like fermion actions
/*! @ingroup actions
 *
 * Even-odd with a unit diagonal preconditioned like Staggered-like fermion actions
 */
class DiagEvenOddPrecStaggeredTypeFermAct : public EvenOddPrecStaggeredTypeFermAct
{
public:
  //! Compute quark propagator
  /*! NOTE: maybe this should produce a quark prop foundry */
  void qprop(LatticeFermion& psi, 
	     const multi1d<LatticeColorMatrix>& u, 
	     const LatticeFermion& chi, 
	     const Real& RsdCG, 
	     int MaxCG, int& ncg_had) const;
};


#endif
