// -*- C++ -*-
// $Id: zolotarev5d_fermact_array_w.h,v 1.11 2004-09-09 15:51:31 edwards Exp $
/*! \file
 *  \brief Unpreconditioned extended-Overlap (5D) (Naryanan&Neuberger) action
 */

#ifndef __zolotarev5d_fermact_array_w_h__
#define __zolotarev5d_fermact_array_w_h__

#include "fermact.h"
#include "actions/ferm/fermacts/overlap_state.h"
#include "actions/ferm/linop/lgherm_w.h"
#include "io/overlap_state_info.h"

using namespace QDP;

namespace Chroma
{
  //! Name and registration
  namespace Zolotarev5DFermActEnv
  {
    extern const std::string name;
    extern const bool registered;
  }


  //! Params for 5D overlap ferm acts
  struct Zolotarev5DFermActParams
  {
    Zolotarev5DFermActParams(XMLReader& in, const std::string& path);
  
    Real Mass;
    int RatPolyDeg;
    ZolotarevStateInfo StateInfo;
  
    std::string AuxFermAct;
  };


  // Reader/writers
  void read(XMLReader& xml, const string& path, Zolotarev5DFermActParams& param);
  void write(XMLReader& xml, const string& path, const Zolotarev5DFermActParams& param);


  //! 5D continued fraction overlap action (Borici,Wenger, Edwards)

  /*!
   * \ingroup fermact
   *
   * This operator applies the extended version of the hermitian overlap operator
   *   Chi  =   ((1+Mass)/(1-Mass)*gamma_5 + B) . Psi
   *  where  B  is the continued fraction of the zolotarev approx. to eps(H(m))
   */
  class Zolotarev5DFermActArray : public UnprecWilsonTypeFermAct< multi1d<LatticeFermion> >
  {
  public:

    Zolotarev5DFermActArray(Handle< FermBC< multi1d< LatticeFermion> > > fbc_, 
			    Handle<UnprecWilsonTypeFermAct<LatticeFermion> > S_aux_,
			    Real& Mass_,
			    int RatPolyDeg_,
			    XMLWriter& writer_) :
      fbc(fbc_), S_aux(S_aux_), Mass(Mass_), RatPolyDeg(RatPolyDeg_), writer(writer_)  {
    
      if ( RatPolyDeg_ % 2 == 0 ) { 
	QDP_error_exit("For Now (and possibly forever), 5D Operators can only be constructed with ODD approximation order. You gave an even one: =%d\n", RatPolyDeg_);
      }
      N5 = RatPolyDeg_;
    }

    // Construct the action out of a parameter structure
    Zolotarev5DFermActArray(Handle< FermBC< multi1d< LatticeFermion> > > fbc_a_,
			    Handle< FermBC< LatticeFermion > > fbc_,
			    const Zolotarev5DFermActParams& param,
			    XMLWriter& writer_);
  

    //! Copy constructor
    Zolotarev5DFermActArray(const Zolotarev5DFermActArray& a) : 
      fbc(a.fbc), S_aux(a.S_aux), Mass(a.Mass), RatPolyDeg(a.RatPolyDeg), writer(a.writer), N5(a.N5) {};

    //! Assignment
    /* Writer screws this up -- I might get rid of that 
       Zolotarev5DFermActArray& operator=(const Zolotarev5DFermActArray& a) {
       fbc=a.fbc; 
       S_aux=a.S_aux;
       Mass=a.Mass;
       RatPolyDeg = a.RatPolyDeg;
       writer=a.writer;
       return *this;
       }
    */

    //! Return the fermion BC object for this action
    const FermBC< multi1d< LatticeFermion> >& getFermBC() const {return *fbc;}

    int size(void) const { return N5; }

    //! Return the quark mass
    Real quark_mass() const {return Mass;}

    //! Produce a linear operator for this action
    const LinearOperator< multi1d<LatticeFermion> >* linOp(Handle<const ConnectState> state) const;

    //! Produce a linear operator for this action
    const LinearOperator< multi1d<LatticeFermion> >* lnonHermLinOp(Handle<const ConnectState> state) const;

    //! Produce a linear operator M^dag.M for this action
    const LinearOperator< multi1d<LatticeFermion> >* lMdagM(Handle<const ConnectState> state) const;

    //! Produce a linear operator M^dag.M for this action
    const LinearOperator< multi1d<LatticeFermion> >* lnonHermMdagM(Handle<const ConnectState> state) const;

    //! produce gamma_5 times M 
    const LinearOperator< multi1d<LatticeFermion> >* gamma5HermLinOp(Handle<const ConnectState> state) const {
      return new lgherm<multi1d<LatticeFermion> >(linOp(state));
    }
    //! Compute quark propagator over base type
    /*! 
     * Solves  M.psi = chi
     *
     * \param psi      quark propagator ( Modify )
     * \param u        gauge field ( Read )
     * \param chi      source ( Modify )
     * \param invParam inverter parameters ( Read (
     * \param ncg_had  number of CG iterations ( Write )
     *
     * NOTE: maybe this should produce a quark prop foundry class object 
     */

    void qprop(LatticeFermion& psi, 
	       Handle<const ConnectState> state, 
	       const LatticeFermion& chi, 
	       const InvertParam_t& invParam,
	       int& ncg_had) const;

    //! Destructor is automatic
    ~Zolotarev5DFermActArray() {}


    // Create state functions

    //! Given links, create the state needed for the linear operators
    /*! Override the parent */
    const OverlapConnectState<LatticeFermion>*
    createState(const multi1d<LatticeColorMatrix>& u_) const ;


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


    // Gauge field and whatever (min/max, e-vectors)
    const OverlapConnectState<LatticeFermion>*
    createState(const multi1d<LatticeColorMatrix>& u_,
		const OverlapStateInfo& state_info,
		XMLWriter& xml_out,
		Real wilsonMass=16) const;


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

    Real Mass;
    int RatPolyDeg;
    int  N5;
  
    XMLWriter& writer;
  };

}

#endif
