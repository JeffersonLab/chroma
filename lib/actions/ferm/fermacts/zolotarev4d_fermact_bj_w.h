// -*- C++ -*-
// $Id: zolotarev4d_fermact_bj_w.h,v 1.1 2003-12-09 17:44:47 bjoo Exp $

/*! \file
 *  \brief 4D Zolotarev variant of Overlap-Dirac operator
 */

#ifndef __zolotarev4d_fermact_w_bj_h__
#define __zolotarev4d_fermact_w_bj_h__

#include "fermact.h"
#include "actions/ferm/fermacts/overlap_fermact_base_w.h"
#include "actions/ferm/fermacts/zolotarev_state.h"

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

class Zolotarev4DFermActBj : public OverlapFermActBase
{
public:
  //! Full constructor
  Zolotarev4DFermActBj(const UnprecWilsonTypeFermAct<LatticeFermion>& Mact_, 
		       const Real& m_q_,
		       int RatPolyDeg_,
		       Real RsdCGinner_,
		       int MaxCGinner_,
		       XMLBufferWriter& writer_) :
    Mact(Mact_), m_q(m_q_), RatPolyDeg(RatPolyDeg_), RsdCGinner(RsdCGinner_), MaxCGinner(MaxCGinner_), writer(writer_)
    {}

  //! Return the quark mass
  Real quark_mass() const {return m_q;}

  //! Does this object really satisfy the Ginsparg-Wilson relation?
  /*! HACK - NEED TO FIX THIS */
  bool isChiral() const {return false;}

  XMLBufferWriter& getWriter() const { 
    return writer;
  }

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
  ~Zolotarev4DFermActBj() {}

protected:
  //! Helper in construction
  void init(int& numroot, 
	    Real& coeffP, 
 	    multi1d<Real>& resP, 
	    multi1d<Real>& rootQ, 
	    int& NEig, 
	    multi1d<Real>& EigValFunc,
	    const ZolotarevConnectState<LatticeFermion>& state) const;

private:
  //! Partial constructor not allowed
  Zolotarev4DFermActBj();

private:
  // Auxilliary action used for kernel of operator
  const UnprecWilsonTypeFermAct<LatticeFermion>& Mact;
  Real m_q;
  int RatPolyDeg;   // degree of polynomial-ratio
  Real RsdCGinner;  // Residuum for inner iterations -- property of action
                    // not state. (eg 5D op has none of this)
  int MaxCGinner;
  XMLBufferWriter& writer;
};

#endif
