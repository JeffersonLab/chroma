// -*- C++ -*-
// $Id: zolotarev5d_fermact_array_w.h,v 1.3 2004-01-13 13:14:49 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned extended-Overlap (5D) (Naryanan&Neuberger) action
 */

#ifndef __zolotarev5d_fermact_array_w_h__
#define __zolotarev5d_fermact_array_w_h__

#include "fermact.h"
#include "actions/ferm/fermacts/overlap_state.h"
using namespace QDP;

//! 5D continued fraction overlap action (Borici,Wenger, Edwards)
/*!
 * \ingroup fermact
 *
 * This operator applies the extended version of the hermitian overlap operator
 *   Chi  =   ((1+m_q)/(1-m_q)*gamma_5 + B) . Psi
 *  where  B  is the continued fraction of the zolotarev approx. to eps(H(m))
 */

class Zolotarev5DFermActArray : public UnprecWilsonTypeFermAct< multi1d<LatticeFermion> >
{
public:

  Zolotarev5DFermActArray(Handle< FermBC< multi1d< LatticeFermion> > > fbc_, 
			  Handle<UnprecWilsonTypeFermAct<LatticeFermion> > S_aux_,
			  Real& m_q_,
	 		  int RatPolyDeg_,
			  XMLBufferWriter& writer_) :
    fbc(fbc_), S_aux(S_aux_), m_q(m_q_), RatPolyDeg(RatPolyDeg_), writer(writer_)  {
    
    if ( RatPolyDeg_ % 2 == 0 ) { 
      QDP_error_exit("For Now (and possibly forever), 5D Operators can only be constructed with ODD approximation order. You gave an even one: =%d\n", RatPolyDeg_);
    }
    N5 = RatPolyDeg_;
  }



  //! Copy constructor
  Zolotarev5DFermActArray(const Zolotarev5DFermActArray& a) : 
    fbc(a.fbc), S_aux(a.S_aux), m_q(a.m_q), RatPolyDeg(a.RatPolyDeg), writer(a.writer), N5(a.N5) {};

  //! Assignment
  /* Writer screws this up -- I might get rid of that 
  Zolotarev5DFermActArray& operator=(const Zolotarev5DFermActArray& a) {
    fbc=a.fbc; 
    S_aux=a.S_aux;
    m_q=a.m_q;
    RatPolyDeg = a.RatPolyDeg;
    writer=a.writer;
    return *this;
  }
  */

  //! Return the fermion BC object for this action
  const FermBC< multi1d< LatticeFermion> >& getFermBC() const {return *fbc;}

  int size(void) const { return N5; }

  //! Return the quark mass
  Real quark_mass() const {return m_q;}

  //! Produce a linear operator for this action
  const LinearOperator< multi1d<LatticeFermion> >* linOp(Handle<const ConnectState> state) const;

  //! Produce a linear operator M^dag.M for this action
  const LinearOperator< multi1d<LatticeFermion> >* lMdagM(Handle<const ConnectState> state) const;

  //! Compute quark propagator over base type
  /*! 
   * Solves  M.psi = chi
   *
   * \param psi      quark propagator ( Modify )
   * \param u        gauge field ( Read )
   * \param chi      source ( Modify )
   * \param invType  inverter type ( Read (
   * \param RsdCG    CG (or MR) residual used here ( Read )
   * \param MaxCG    maximum number of CG iterations ( Read )
   * \param ncg_had  number of CG iterations ( Write )
   *
   * NOTE: maybe this should produce a quark prop foundry class object 
   */

  void qprop(LatticeFermion& psi, 
	     Handle<const ConnectState> state, 
	     const LatticeFermion& chi, 
	     enum InvType invType,
	     const Real& RsdCG, 
    int MaxCG, int& ncg_had) const {}

  //! Destructor is automatic
  ~Zolotarev5DFermActArray() {}


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


protected:
  //! Helper in construction
  void init(Real& scale_fac,
	    multi1d<Real>& alpha,
	    multi1d<Real>& beta,
	    int& NEig,
	    multi1d<Real>& EigValFunc,
	    const OverlapConnectState<LatticeFermion>& state) const;
private:
  // Hide partial constructor
  Zolotarev5DFermActArray();

private:
  Handle< FermBC< multi1d<LatticeFermion> > >  fbc;
  Handle< UnprecWilsonTypeFermAct<LatticeFermion> > S_aux;

  Real m_q;
  int RatPolyDeg;
  int  N5;
  
  XMLBufferWriter& writer;
};

#endif
