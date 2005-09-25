// -*- C++ -*-
// $Id: ovlap_partfrac4d_fermact_w.h,v 2.0 2005-09-25 21:04:25 edwards Exp $

/*! \file
 *  \brief 4D Zolotarev variant of Overlap-Dirac operator
 */

#ifndef __ovlap_partfrac4d_fermact_w_h__
#define __ovlap_partfrac4d_fermact_w_h__

#include "fermact.h"
#include "actions/ferm/fermacts/overlap_fermact_base_w.h"
#include "actions/ferm/fermacts/overlap_state.h"
#include "meas/eig/eig_w.h"
#include "io/overlap_state_info.h"
#include "io/enum_io/enum_io.h"


namespace Chroma
{
  //! Name and registration
  namespace OvlapPartFrac4DFermActEnv
  {
    extern const std::string name;
    extern const bool registered;
  }


  //! Params for overlap ferm acts
  struct OvlapPartFrac4DFermActParams
  {
    OvlapPartFrac4DFermActParams() : ReorthFreqInner(10), inner_solver_type(OVERLAP_INNER_CG_SINGLE_PASS) {};
    OvlapPartFrac4DFermActParams(XMLReader& in, const std::string& path);
    
    Real Mass;
    Real AuxMass;  // Filled when AuxFermAct is processed
    int RatPolyDeg;
    int RatPolyDegPrecond;
    int ReorthFreqInner;

    //! Type of approximation ZOLOTAREV or TANH
    CoeffType approximation_type;
    Real approxMin;
    Real approxMax;

    InvertParam_t invParamInner;
    OverlapInnerSolverType inner_solver_type;

    bool isChiralP;    
    std::string AuxFermAct;
    std::string AuxFermActGrp;
  };


  // Reader/writers
  void read(XMLReader& xml, const string& path, OvlapPartFrac4DFermActParams& param);
  void write(XMLWriter& xml, const string& path, const OvlapPartFrac4DFermActParams& param);


  //! 4D Zolotarev variant of Overlap-Dirac operator
  /*!
   * \ingroup fermacts
   *
   * This routine is specific to Wilson-like fermions!
   *
   * NOTE: for now we assume the kernel is a fund. rep. fermion type,
   * but that is not necessary
   */
  class OvlapPartFrac4DFermAct : public OverlapFermActBase
  {
  public:
    //! Full constructor
    /*
    OvlapPartFrac4DFermAct(Handle<FermBC<LatticeFermion> > fbc_,
		       Handle<UnprecWilsonTypeFermAct<LatticeFermion> > Mact_, 		       const Real& Mass_,
		       const int RatPolyDeg_,
		       const Real& RsdCGinner_,
		       int MaxCGinner_,
		       XMLWriter& writer_,
		       const int ReorthFreqInner_=10,
		       const OverlapInnerSolverType inner_solver_type_=OVERLAP_INNER_CG_SINGLE_PASS
      ) : writer(writer_)
      {
	fbc = fbc_;
	Mact = Mact_;
	params.Mass = Mass_;
	params.RatPolyDeg = RatPolyDeg_; 
	params.invParamInner.RsdCG = RsdCGinner_; 
	params.invParamInner.MaxCG = MaxCGinner_;
	params.ReorthFreqInner = ReorthFreqInner_; 
	params.inner_solver_type = inner_solver_type_;

	// Default Preconditioner degree is RatPolyDeg
	params.RatPolyDegPrecond = RatPolyDeg_;
      }
    */

    /*
    OvlapPartFrac4DFermAct(Handle<FermBC<LatticeFermion> > fbc_,
		       Handle<UnprecWilsonTypeFermAct<LatticeFermion> > Mact_, 
		       const Real& Mass_,
		       const int RatPolyDeg_,
		       const int RatPolyDegPrecond_,
		       const Real& RsdCGinner_,
		       int MaxCGinner_,
		       XMLWriter& writer_,
		       const int ReorthFreqInner_=10,
		       const OverlapInnerSolverType inner_solver_type_=OVERLAP_INNER_CG_SINGLE_PASS
      ) : writer(writer_)
      {
	fbc = fbc_; Mact = Mact_; params.Mass = Mass_; params.RatPolyDeg = RatPolyDeg_; 
	params.RatPolyDegPrecond = RatPolyDegPrecond_;
	params.invParamInner.RsdCG = RsdCGinner_; 
	params.invParamInner.MaxCG = MaxCGinner_;
	params.ReorthFreqInner = ReorthFreqInner_; params.inner_solver_type = inner_solver_type_;
      }
    */

    //! Construct from param struct
    OvlapPartFrac4DFermAct(Handle<FermBC<LatticeFermion> > fbc_,
		       const OvlapPartFrac4DFermActParams& params);


    //! Copy Constructor
    OvlapPartFrac4DFermAct(const OvlapPartFrac4DFermAct& a) : 
      fbc(a.fbc), params(a.params)
      {
	Mact = a.Mact;
      }
  

    //! Virtual copy constructor
    OvlapPartFrac4DFermAct* clone() const {return new OvlapPartFrac4DFermAct(*this);}

    // Assignment
    OvlapPartFrac4DFermAct& operator=(const OvlapPartFrac4DFermAct& a) 
    {
      fbc = a.fbc;
      Mact = a.Mact;
      params = a.params;
      return *this;
    }
   

    //! Return the fermion BC object for this action
    const FermBC<LatticeFermion>& getFermBC() const {return *fbc;}

    //! Return the quark mass
    Real getQuarkMass() const {return params.Mass;}


    //! Is the operator Chiral 
    /*! The operator is chiral if it satisfies the GW 
     *  relation (or a massive version of it). It is certainly 
     *  the intention that this operator be chiral in this sense.
     *  However, setting it up wrongly may make it non chiral.
     *  that would need a run-time check. So this is a hack below,
     *  signifying intent */
    bool isChiral() const { return params.isChiralP; }
    
    // Create state functions
    
    //! Create OverlapConnectState from XML
    const OverlapConnectState* 
    createState(const multi1d<LatticeColorMatrix>& u, 
		XMLReader& state_info_xml,
		const string& state_info_path) const;
    
    //! Given links, create the state needed for the linear operators
    /*! Override the parent */
    
    //! Create a ConnectState with just the gauge fields
    const OverlapConnectState*
    createState(const multi1d<LatticeColorMatrix>& u_) const ;

    //! Create a ConnectState with just the gauge fields, and a lower
    //!  approximation bound
    const OverlapConnectState*
    createState(const multi1d<LatticeColorMatrix>& u_, 
		const Real& approxMin_) const ;
 
    
    //! Create a connect State with just approximation range bounds
    const OverlapConnectState*
    createState(const multi1d<LatticeColorMatrix>& u_, 
		const Real& approxMin_, 
		const Real& approxMax_) const;


    //! Create OverlapConnectState with eigenvalues/vectors 
    const OverlapConnectState*
    createState(const multi1d<LatticeColorMatrix>& u_,
		const multi1d<Real>& lambda_lo_,
		const multi1d<LatticeFermion>& evecs_lo_, 
		const Real& lambda_hi_) const;


    //! Create from OverlapStateInfo Structure
    const OverlapConnectState*
    createState(const multi1d<LatticeColorMatrix>& u_,
		const OverlapStateInfo& state_info) const;

    //! Produce a linear operator for this action
    /*! 
     * NOTE: the arg MUST be the original base because C++ requires it for a virtual func!
     * The function will have to downcast to get the correct state
     */
    const UnprecLinearOperator< LatticeFermion, multi1d<LatticeColorMatrix> >* 
    unprecLinOp(Handle<const ConnectState> state, const Real& m_q) const;

    const LinearOperator<LatticeFermion>* 
    linOpPrecondition(Handle<const ConnectState > state) const;

    //! Produce a linear operator M^dag.M for this action
    /*! 
     * NOTE: the arg MUST be the original base because C++ requires it for a
     * virtual func!
     * The function will have to downcast to get the correct state
     */
    const LinearOperator<LatticeFermion>* 
    lMdagM(Handle<const ConnectState> state) const;

    //! Produce a linear operator M^dag.M for this action to be applied
    //  to a vector of known chirality. Chirality is passed in
    const LinearOperator<LatticeFermion>* 
    lMdagM(Handle<const ConnectState> state, const Chirality& chirality) const;

    //! Produce a linear operator M^dag.M for this action
    /*! 
     * NOTE: the arg MUST be the original base because C++ requires 
     * it for a virtual func!
     * The function will have to downcast to get the correct state
     */
    const LinearOperator<LatticeFermion>* 
    lMdagMPrecondition(Handle<const ConnectState> state) const;

    //! Produce a linear operator M^dag.M for this action to be applied
    //  to a vector of known chirality. Chirality is passed in
    const LinearOperator<LatticeFermion>* 
    lMdagMPrecondition(Handle<const ConnectState> state, 
		       const Chirality& chirality) const;


    //! Produce a linear operator that gives back gamma_5 eps(H)
    const LinearOperator<LatticeFermion>* 
    lgamma5epsH(Handle<const ConnectState> state) const;

    //! Produce a linear operator that gives back gamma_5 eps(H)
    const LinearOperator<LatticeFermion>* 
    lgamma5epsHPrecondition(Handle<const ConnectState> state) const;
    
    //! Destructor is automatic
    ~OvlapPartFrac4DFermAct() {}

  protected:
    //! Helper in construction
    void init(int& numroot, 
	      Real& coeffP, 
	      multi1d<Real>& resP, 
	      multi1d<Real>& rootQ, 
	      int& NEig, 
	      multi1d<Real>& EigValFunc,
	      const OverlapConnectState& state) const;

    //! Construct stuff but use RatPolyDegPrec in the polynomial
    void initPrec(int& numroot, 
		  Real& coeffP, 
		  multi1d<Real>& resP, 
		  multi1d<Real>& rootQ, 
		  int& NEig, 
		  multi1d<Real>& EigValFunc,
		  const OverlapConnectState& state) const;


  private:
    //!  Partial constructor not allowed
    OvlapPartFrac4DFermAct();

  private:
    Handle<FermBC<LatticeFermion> >  fbc;   // fermion BC
    // Auxilliary action used for kernel of operator
    Handle<UnprecWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> > > Mact;   
    OvlapPartFrac4DFermActParams params;
  };

}
#endif
