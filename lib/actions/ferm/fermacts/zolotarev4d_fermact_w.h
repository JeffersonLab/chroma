// -*- C++ -*-
// $Id: zolotarev4d_fermact_w.h,v 1.12 2003-12-30 17:27:15 bjoo Exp $

/*! \file
 *  \brief 4D Zolotarev variant of Overlap-Dirac operator
 */

#ifndef __zolotarev4d_fermact_w_h__
#define __zolotarev4d_fermact_w_h__

#include "fermact.h"
#include "actions/ferm/fermacts/overlap_fermact_base_w.h"
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

class Zolotarev4DFermAct : public OverlapFermActBase
{
public:
  //! Full constructor
  Zolotarev4DFermAct(const UnprecWilsonTypeFermAct<LatticeFermion>& Mact_, 
		     const Real& m_q_,
		     int RatPolyDeg_) :
    Mact(Mact_), m_q(m_q_), RatPolyDeg(RatPolyDeg_)
    {}

  //! Return the quark mass
  Real quark_mass() const {return m_q;}

  //! Does this object really satisfy the Ginsparg-Wilson relation?
  /*! HACK - NEED TO FIX THIS */
  bool isChiral() const {return false;}

  //! Default version, given links create the state needed for the linear 
  //  operators
  /*! Covariant Return Rule - overload base class virtual function with new
    return type */
  const EVConnectStateBase<LatticeFermion>* 
    createState(const multi1d<LatticeColorMatrix>& u) const;

  //! Full version, given links create the state needed for the linear operators
  /*! New function */

  const EVConnectStateBase<LatticeFermion>* 
    createState(const multi1d<LatticeColorMatrix>& u, 
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
  void init(int& numroot, 
	    Real& coeffP, 
	    multi1d<Real>& resP, 
	    multi1d<Real>& rootQ, 
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
