// -*- C++ -*-
// $Id: zolotarev4d_fermact_w.h,v 1.7 2003-12-02 22:35:26 edwards Exp $

/*! \file
 *  \brief 4D Zolotarev variant of Overlap-Dirac operator
 */

#ifndef __zolotarev4d_fermact_w_h__
#define __zolotarev4d_fermact_w_h__

#include "fermact.h"
#include "actions/ferm/fermacts/ev_state.h"

using namespace QDP;

//! 4D Zolotarev variant of Overlap-Dirac operator
/*!
 * \ingroup fermact
 *
 * This routine is specific to Wilson-like fermions!
 *
 * NOTE: for now we assume the kernel is a fund. rep. fermion type,
 * but that is not necessary
 */

class Zolotarev4DFermAct : public UnprecWilsonTypeFermAct<LatticeFermion>
{
public:
  //! Full constructor
  Zolotarev4DFermAct(const UnprecWilsonTypeFermAct<LatticeFermion>& M_, const Real& m_q_) :
    M(M_), m_q(m_q_)
    {init();}

  //! Default version, given links create the state needed for the linear operators
  const EVConnectState<LatticeFermion>* createState(const multi1d<LatticeColorMatrix>& u) const;

  //! Full version, given links create the state needed for the linear operators
  const EVConnectState<LatticeFermion>* createState(const multi1d<LatticeColorMatrix>& u, 
						    const multi1d<LatticeFermion>& EigVec, 
						    const multi1d<Real>& EigVal) const;

  //! Produce a linear operator for this action
  const LinearOperator<LatticeFermion>* linOp(const EVConnectState<LatticeFermion>& state) const;

  //! Produce a linear operator M^dag.M for this action
  const LinearOperator<LatticeFermion>* lMdagM(const EVConnectState<LatticeFermion>& state) const;

  //! Destructor is automatic
  ~Zolotarev4DFermAct() {}

protected:
  //! Helper in construction
  void init();

private:
  //! Partial constructor not allowed
  Zolotarev4DFermAct();

private:
  Real m_q;
  const UnprecWilsonTypeFermAct<LatticeFermion>& M;
};

#endif
