// -*- C++ -*-
// $Id: zolotarev4d_fermact_w.h,v 1.17 2004-01-13 10:00:57 bjoo Exp $

/*! \file
 *  \brief 4D Zolotarev variant of Overlap-Dirac operator
 */

#ifndef __zolotarev4d_fermact_w_h__
#define __zolotarev4d_fermact_w_h__

#include "fermact.h"
#include "actions/ferm/fermacts/overlap_fermact_base_w.h"
#include "actions/ferm/fermacts/overlap_state.h"
#include "meas/eig/eig_w.h"

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
  Zolotarev4DFermAct(Handle<FermBC<LatticeFermion> > fbc_,
		       Handle<UnprecWilsonTypeFermAct<LatticeFermion> > Mact_, 
		       const Real& m_q_,
		       const int RatPolyDeg_,
		       const Real& RsdCGinner_,
		       int MaxCGinner_,
		       XMLBufferWriter& writer_) :
    fbc(fbc_), Mact(Mact_), m_q(m_q_), RatPolyDeg(RatPolyDeg_), 
    RsdCGinner(RsdCGinner_), MaxCGinner(MaxCGinner_), writer(writer_)
    {}

  //! Copy Constructor
  Zolotarev4DFermAct(const Zolotarev4DFermAct& a) :
    fbc(a.fbc), Mact(a.Mact), m_q(a.m_q), RatPolyDeg(a.RatPolyDeg),
    RsdCGinner(a.RsdCGinner), MaxCGinner(a.MaxCGinner), writer(a.writer)
  {};
  

  // Assignment
  /* The writer makes this duff not supporting copying or assignment
  Zolotarev4DFermAct& operator=(const Zolotarev4DFermAct& a) 
  {
    fbc = a.fbc;
    Mact = a.Mact;
    m_q = a.m_q;
    RatPolyDeg=a.RatPolyDeg;
    RsdCGinner=a.RsdCGinner;
    MaxCGinner=a.MaxCGinner;
    writer = a.writer;
    return *this;
  }
  */

  //! Return the fermion BC object for this action
  const FermBC<LatticeFermion>& getFermBC() const {return *fbc;}

  //! Return the quark mass
  Real quark_mass() const {return m_q;}


  //! Is the operator Chiral 
  /*! The operator is chiral if it satisfies the GW 
   *  relation (or a massive version of it). It is certainly 
   *  the intention that this operator be chiral in this sense.
   *  However, setting it up wrongly may make it non chiral.
   *  that would need a run-time check. So this is a hack below,
   *  signifying intent */
  bool isChiral() const { return true; }

  XMLBufferWriter& getWriter() const { 
    return writer;
  }

  // Create state functions

  // Just gauge field and epsilon -- Approx Max is 2*Nd 
  const OverlapConnectState<LatticeFermion>*
  createState(const multi1d<LatticeColorMatrix>& u_, 
	      const Real& approxMin_) const ;
 
  // Gauge field, epsilon, approx min, approx max
  const OverlapConnectState<LatticeFermion>*
  createState(const multi1d<LatticeColorMatrix>& u_, 
	      const Real& approxMin_, 
	      const Real& approxMax_) const;


  // Gauge field, e-values, e-vectors
  const OverlapConnectState<LatticeFermion>*
  createState(const multi1d<LatticeColorMatrix>& u_,
	      const multi1d<Real>& lambda_lo_,
	      const multi1d<LatticeFermion>& evecs_lo_, 
	      const Real& lambda_hi_) const;


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
  const LinearOperator<LatticeFermion>* lMdagM(Handle<const ConnectState> state) const;

  //! Produce a linear operator M^dag.M for this action to be applied
  //  to a vector of known chirality. Chirality is passed in
  const LinearOperator<LatticeFermion>* lMdagM(Handle<const ConnectState> state, const Chirality& chirality) const;


  // Special qprop for now
  /*
  void qprop(LatticeFermion& psi, 
	     Handle<const ConnectState> state, 
	     const LatticeFermion& chi, 
	     enum InvType invType,
	     const Real& RsdCG, 
	     int MaxCG, int& ncg_had) const;
  */
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
	    const OverlapConnectState<LatticeFermion>& state) const;

private:
  //!  Partial constructor not allowed
  Zolotarev4DFermAct();

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
