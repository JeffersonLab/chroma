// -*- C++ -*-
// $Id: zolotarev4d_fermact_w.h,v 1.2 2003-10-10 03:46:46 edwards Exp $

/*! \file
 *  \brief 4D Zolotarev variant of Overlap-Dirac operator
 */

#ifndef __zolotarev4d_fermact_w_h__
#define __zolotarev4d_fermact_w_h__

#include "fermact.h"
#include "linearop.h"

using namespace QDP;

//! 4D Zolotarev variant of Overlap-Dirac operator
/*!
 * \ingroup fermact
 *
 * This routine is specific to Wilson fermions!
 *
 */

class Zolotarev4DFermAct : UnprecWilsonTypeFermAct 
{
public:
  //! Full constructor
  Zolotarev4DFermAct(const Real& _m_q, const UnprecWilsonTypeFermAct& _M);

  //! Full constructor including eigenvectors
  /*! 
   * NOTE: here must include gauge field, otherwise eigenv. make no sense 
   *
   * NOTE: could/should generalize this to some object holding underlying kernel
   * and eigenmodes
   */
  Zolotarev4dLinOp(const multi1d<LatticeColorMatrix>& _u, 
		   const Real& _m_q, const UnprecWilsonTypeFermAct& _M,
		   const multi1d<LatticeFermion>& _EigVec, const multi1d<Real>& _EigVal);

  //! Produce a linear operator for this action
  const LinearOperator* linOp(const multi1d<LatticeColorMatrix>& u) const;

  //! Produce a linear operator M^dag.M for this action
  const LinearOperator* lMdagM(const multi1d<LatticeColorMatrix>& u) const;

  //! Compute dS_f/dU
  multi1d<LatticeColorMatrix> dsdu(const multi1d<LatticeColorMatrix>& u) const;

  //! Destructor is automatic
  ~Zolotarev4DFermAct() {}

private:
  //! Partial constructor is hidden
  Zolotarev4DFermAct() {}

private:
  Real m_q;
  UnprecWilsonTypeFermAct& M;

#endif
