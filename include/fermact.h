// -*- C++ -*-
// $Id: fermact.h,v 1.1 2003-04-09 01:21:01 edwards Exp $

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
//  virtual LinOpFoundry* linop() const = 0;
#endif

  //! Produce a linear operator for this action
  /*! NOTE: maybe this should be abstracted to a foundry class object */
  virtual LinearOperator* linop() const = 0;

  //! Produce a linear operator M^dag.M for this action
  /*! NOTE: maybe this should be abstracted to a foundry class object */
  virtual LinearOperator* lmdagm() const = 0;

  //! Compute quark propagator
  /*! NOTE: maybe this should produce a quark prop foundry class object */
  virtual LatticeFermion Qprop(const multi1d<LatticeComplex>& u) const = 0;

  //! Compute dS_f/dU
  /*! NOTE: maybe this should produce a derivative foundry class object */
  virtual multi1d<LatticeColorMatrix> dsdu(const multi1d<LatticeComplex>& u) const = 0;

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
  virtual LatticeFermion Qprop(const multi1d<LatticeComplex>& u) const;
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
  virtual LatticeFermion Qprop(const multi1d<LatticeComplex>& u) const;
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
  virtual LatticeFermion Qprop(const multi1d<LatticeComplex>& u) const;
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
  virtual LatticeFermion Qprop(const multi1d<LatticeComplex>& u) const;
};


#endif
