// -*- C++ -*-
// $Id: zolotarev4d_linop_w.h,v 1.6 2003-11-20 05:43:41 edwards Exp $

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

class Zolotarev4DLinOp : public LinearOperator<LatticeFermion>
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
  const OrderedSubset& subset() const {return all;}

  //! Creation routine
  void create(const multi1d<LatticeColorMatrix>& _u, const Real& _m_q);

  //! Full creation routine including eigenvectors
  void create(const multi1d<LatticeColorMatrix>& _u, const Real& _m_q,
	      const multi1d<LatticeFermion>& _EigVec, const multi1d<Real>& _EigVal);

  //! Apply the operator onto a source vector
  void operator() (LatticeFermion& chi, const LatticeFermion& psi, enum PlusMinus isign) const;

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
  LinearOperator<LatticeFermion>* A;
};

#endif
