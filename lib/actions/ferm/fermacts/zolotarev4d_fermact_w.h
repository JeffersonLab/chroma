// -*- C++ -*-
// $Id: zolotarev4d_fermact_w.h,v 1.9 2003-12-03 06:12:02 edwards Exp $

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
  Zolotarev4DFermAct(const UnprecWilsonTypeFermAct<LatticeFermion>& Mact_, const Real& m_q_,
		     int RatPolyDeg_) :
    Mact(Mact_), m_q(m_q_), RatPolyDeg(RatPolyDeg_)
    {}

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
  void init(int& numroot, Real& coeffP, multi1d<Real>& resP, multi1d<Real>& rootQ, 
	    int& NEig, 
	    multi1d<Real>& EigValFunc,
	    const EVConnectState<LatticeFermion>& state,
	    int& MaxCG, Real& RsdCGinner) const;

private:
  //! Partial constructor not allowed
  Zolotarev4DFermAct();

private:
  // Auxilliary action used for kernel of operator
  const UnprecWilsonTypeFermAct<LatticeFermion>& Mact;
  Real m_q;
  int RatPolyDeg;   // degree of polynomial-ratio
};

#endif
