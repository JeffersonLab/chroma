// -*- C++ -*-
// $Id: zolotarev4d_linop_w.h,v 1.7 2003-12-02 22:35:46 edwards Exp $

/*! \file
 *  \brief 4D Zolotarev variant of Overlap-Dirac operator
 */

#ifndef __zolotarev4d_linop_w_h__
#define __zolotarev4d_linop_w_h__

#include "fermact.h"
#include "actions/ferm/fermacts/ev_state.h"

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
  Zolotarev4DLinOp(const EVConnectState<LatticeFermion>& state_, 
		   const UnprecWilsonTypeFermAct<LatticeFermion>& M_, 
		   const Real& m_q_) : state(state_), M(M_), m_q(m_q_)
    {init();}

  //! Destructor
  ~Zolotarev4DLinOp();

  //! Only defined on the odd subset
  const OrderedSubset& subset() const {return all;}

  //! Apply the operator onto a source vector
  void operator() (LatticeFermion& chi, const LatticeFermion& psi, enum PlusMinus isign) const;

protected:
  //! Internal creation routine
  void init();

private:
  Real m_q;
  Real RsdCGinner;   // inner CG accuracy

  // Keep a private copy (actually a link) to eigenvectors/values
  const EVConnectState<LatticeFermion>  state;
  int NEigVal;

  // Underlying fermion operator
  const UnprecWilsonTypeFermAct<LatticeFermion>& M;
};

#endif
