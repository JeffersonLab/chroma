// -*- C++ -*-
// $Id: zolotarev4d_fermact_bj_w.h,v 1.4 2004-01-02 03:19:41 edwards Exp $

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
  Zolotarev4DFermActBj(Handle<FermBC<LatticeFermion> > fbc_,
		       Handle<UnprecWilsonTypeFermAct<LatticeFermion> > Mact_, 
		       const Real& m_q_,
		       int RatPolyDeg_,
		       const Real& RsdCGinner_,
		       int MaxCGinner_,
		       XMLBufferWriter& writer_) :
    fbc(fbc_), Mact(Mact_), m_q(m_q_), RatPolyDeg(RatPolyDeg_), 
    RsdCGinner(RsdCGinner_), MaxCGinner(MaxCGinner_), writer(writer_)
    {}

  //! Return the fermion BC object for this action
  const FermBC<LatticeFermion>& getFermBC() const {return *fbc;}

  //! Return the quark mass
  Real quark_mass() const {return m_q;}

  //! Does this object really satisfy the Ginsparg-Wilson relation?
  /*! HACK - NEED TO FIX THIS */
  bool isChiral() const {return false;}

  XMLBufferWriter& getWriter() const { 
    return writer;
  }

  // Create state functions

  // Just gauge field and epsilon -- Approx Max is 2*Nd 
  const ZolotarevConnectState<LatticeFermion>*
  createState(const multi1d<LatticeColorMatrix>& u_, 
	      const Real& approxMin_) const ;
 
  // Gauge field, epsilon, approx min, approx max
  const ZolotarevConnectState<LatticeFermion>*
  createState(const multi1d<LatticeColorMatrix>& u_, 
	      const Real& approxMin_, 
	      const Real& approxMax_) const;


  // Gauge field, e-values, e-vectors
  const ZolotarevConnectState<LatticeFermion>*
  createState(const multi1d<LatticeColorMatrix>& u_,
	      const multi1d<Real>& lambda_lo_,
	      const multi1d<LatticeFermion>& evecs_lo_, 
	      const Real& lambda_hi_) const;

  // Gauge field, e-values, e-vectors and explicit approx min and approx max
  /*
  const ZolotarevConnectState<LatticeFermion>*
  createState(const multi1d<LatticeColorMatrix>& u_, 
	      const multi1d<Real>& lambda_lo_,
	      const multi1d<LatticeFermion>& evecs_lo_,
	      const Real& lambda_hi_, 
	      const Real& approxMin_,
	      const Real& approxMax_) const;
  */

  //! Produce a linear operator for this action
  /*! 
   * NOTE: the arg MUST be the original base because C++ requires it for a virtual func!
   * The function will have to downcast to get the correct state
   */
  const LinearOperator<LatticeFermion>* linOp(Handle<const ConnectState> state) const;

  //! Produce a linear operator M^dag.M for this action
  /*! 
   * NOTE: the arg MUST be the original base because C++ requires it for a virtual func!
   * The function will have to downcast to get the correct state
   */
  //const LinearOperator<LatticeFermion>* lMdagM(Handle<const ConnectState> state) const;
  const LinearOperator<LatticeFermion>* lMdagM(Handle<const ConnectState> state) const {};

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
  Handle<FermBC<LatticeFermion> >  fbc;   // fermion BC
  // Auxilliary action used for kernel of operator
  Handle<UnprecWilsonTypeFermAct<LatticeFermion> > Mact;
  Real m_q;
  int RatPolyDeg;   // degree of polynomial-ratio
  Real RsdCGinner;  // Residuum for inner iterations -- property of action
                    // not state. (eg 5D op has none of this)
  int MaxCGinner;
  XMLBufferWriter& writer;
};

#endif
