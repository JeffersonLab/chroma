// -*- C++ -*-
// $Id: zolotarev4d_fermact_w.h,v 1.8 2003-12-03 03:04:44 edwards Exp $

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
  /*! Covariant Return Rule - overload base class virtual function with new return type */
  const EVConnectStateBase<LatticeFermion>* createState(const multi1d<LatticeColorMatrix>& u) const;

  //! Full version, given links create the state needed for the linear operators
  /*! New function */
  const EVConnectStateBase<LatticeFermion>* createState(const multi1d<LatticeColorMatrix>& u, 
							const multi1d<LatticeFermion>& eigVec, 
							const multi1d<Real>& eigVal,
							const Real& eigValMax) const;

  //! Produce a linear operator for this action
  /*! 
   * NOTE: the arg MUST be the original base because C++ requires it for a virtual func!
   * The function will have to downcast to get the correct state
   */
  const LinearOperator<LatticeFermion>* linOp(const ConnectState& state) const;

  //! Produce a linear operator M^dag.M for this action
  /*! 
   * NOTE: the arg MUST be the original base because C++ requires it for a virtual func!
   * The function will have to downcast to get the correct state
   */
  const LinearOperator<LatticeFermion>* lMdagM(const ConnectState& state) const;

  //! Destructor is automatic
  ~Zolotarev4DFermAct() {}

protected:
  //! Helper in construction
  void init();

private:
  //! Partial constructor not allowed
  Zolotarev4DFermAct();

private:
  // Auxilliary action used for kernel of operator
  const UnprecWilsonTypeFermAct<LatticeFermion>& M;
  Real m_q;
};

#endif
