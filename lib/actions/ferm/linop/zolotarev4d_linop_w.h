// -*- C++ -*-
// $Id: zolotarev4d_linop_w.h,v 1.1 2003-04-09 05:57:16 edwards Exp $

/*! \file
 *  \brief 4D Zolotarev variant of Overlap-Dirac operator
 */

#ifndef __zolotarev4d_linop_w_h__
#define __zolotarev4d_linop_w_h__

#include "actions/ferm/fermacts/unprec_wilson_fermact_w.h"

using namespace QDP;

//! 4D Zolotarev variant of Overlap-Dirac operator
/*!
 * \ingroup linop
 *
 * This routine is specific to Wilson fermions!
 *
 */

class Zolotarev4DLinOp : public LinearOperator
{
public:
  //! Full constructor
  Zolotarev4dLinOp(const multi1d<LatticeColorMatrix>& _u, const Real& _m_q)
    {create(_u,_m_q);}

  //! Full constructor including eigenvectors
  Zolotarev4dLinOp(const multi1d<LatticeColorMatrix>& _u, const Real& _m_q,
		   const multi1d<LatticeFermion>& _EigVec, const multi1d<Real>& _EigVal)
    {create(_u,_m_q,_EigVec,_EigVal);}

  //! Destructor
  ~Zolotarev4dLinOp();

  //! Only defined on the odd subset
  const Subset& subset() const {return all;}

  //! Creation routine
  void create(const multi1d<LatticeColorMatrix>& _u, const Real& _m_q);

  //! Full creation routine including eigenvectors
  void create(const multi1d<LatticeColorMatrix>& _u, const Real& _m_q,
	      const multi1d<LatticeFermion>& _EigVec, const multi1d<Real>& _EigVal);

  //! Apply the operator onto a source vector
  LatticeFermion operator() (const LatticeFermion& psi, enum LinOpSign isign) const;

protected:
  //! Internal creation routine
  void make();

private:
  Real m_q;
  multi1d<LatticeColorMatrix> u;

  Real RsdCGinner;   // inner CG accuracy

  // Keep a private copy (actually a link) to eigenvectors/values
  int NEigVal;
  multi1d<Real> EigVal;
  multi1d<LatticeFermion> EigVec;

  /* For now, use standard Wilson. 
   * This can be generalized to pointer to a generic linear operator
   */
  UnpreconditionedWilson M;
  
  /* Underlying lovlapms linear operator */
  LinearOperator* A;
};

#endif
