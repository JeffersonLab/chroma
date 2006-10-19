// -*- C++ -*-
// $Id: ovlap_partfrac4d_fermact_w.h,v 3.4 2006-10-19 16:01:27 edwards Exp $

/*! \file
 *  \brief 4D Zolotarev variant of Overlap-Dirac operator
 */

#ifndef __ovlap_partfrac4d_fermact_w_h__
#define __ovlap_partfrac4d_fermact_w_h__

#include "actions/ferm/fermacts/overlap_fermact_base_w.h"
#include "actions/ferm/fermstates/eigen_state.h"
#include "meas/eig/eig_w.h"
// #include "io/overlap_state_info.h"
#include "io/enum_io/enum_io.h"


namespace Chroma
{
  //! Name and registration
  /*! \ingroup fermacts */
  namespace OvlapPartFrac4DFermActEnv
  {
    extern const std::string name;
    bool registerAll();
  }


  //! Params for overlap ferm acts
  /*! \ingroup fermacts */
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

    struct InvParamInner
    {
      Real RsdCG;
      int  MaxCG;
    } invParamInner;
    OverlapInnerSolverType inner_solver_type;

    bool isChiralP;    
    std::string AuxFermAct;
    std::string AuxFermActGrp;
  };


  // Reader/writers
  /*! \ingroup fermacts */
  void read(XMLReader& xml, const string& path, OvlapPartFrac4DFermActParams& param);
  /*! \ingroup fermacts */
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
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Construct from param struct
    OvlapPartFrac4DFermAct(Handle< FermBC<T,P,Q> > fbc_,
			   const OvlapPartFrac4DFermActParams& params);


    //! Virtual copy constructor
    OvlapPartFrac4DFermAct* clone() const {return new OvlapPartFrac4DFermAct(*this);}

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
    
    //! Create OverlapFermState<T,P,Q> from XML
    EigenConnectState* 
    createState(const multi1d<LatticeColorMatrix>& u, 
		XMLReader& state_info_xml,
		const string& state_info_path) const;
    
    //! Given links, create the state needed for the linear operators
    /*! Override the parent */
    
    //! Create a ConnectState with just the gauge fields
    EigenConnectState* createState(const multi1d<LatticeColorMatrix>& u_) const;


    //! Produce a linear operator for this action
    /*! 
     * NOTE: the arg MUST be the original base because C++ requires it for a virtual func!
     * The function will have to downcast to get the correct state
     */
    UnprecLinearOperator<T,P,Q>* unprecLinOp(Handle< FermState<T,P,Q> > state, 
					     const Real& m_q) const;

    LinearOperator<T>* linOpPrecondition(Handle< FermState<T,P,Q> > state) const;

    //! Produce a linear operator M^dag.M for this action
    /*! 
     * NOTE: the arg MUST be the original base because C++ requires it for a
     * virtual func!
     * The function will have to downcast to get the correct state
     */
    DiffLinearOperator<T,P,Q>* lMdagM(Handle< FermState<T,P,Q> > state) const;

    //! Produce a linear operator M^dag.M for this action to be applied
    //  to a vector of known chirality. Chirality is passed in
    DiffLinearOperator<T,P,Q>* lMdagM(Handle< FermState<T,P,Q> > state, 
				      const Chirality& chirality) const;

    //! Produce a linear operator M^dag.M for this action
    /*! 
     * NOTE: the arg MUST be the original base because C++ requires 
     * it for a virtual func!
     * The function will have to downcast to get the correct state
     */
    LinearOperator<T>* lMdagMPrecondition(Handle< FermState<T,P,Q> > state) const;

    //! Produce a linear operator M^dag.M for this action to be applied
    //  to a vector of known chirality. Chirality is passed in
    LinearOperator<T>* lMdagMPrecondition(Handle< FermState<T,P,Q> > state, 
					  const Chirality& chirality) const;


    //! Produce a linear operator that gives back gamma_5 eps(H)
    LinearOperator<T>* lgamma5epsH(Handle< FermState<T,P,Q> > state) const;

    //! Produce a linear operator that gives back gamma_5 eps(H)
    LinearOperator<T>* lgamma5epsHPrecondition(Handle< FermState<T,P,Q> > state) const;
    
    //! Destructor is automatic
    ~OvlapPartFrac4DFermAct() {}

  protected:
    //! Return the factory object that produces a state
    const CreateFermState<T,P,Q>& getCreateState() const {return *cfs;}

    //! Helper in construction
    void init(int& numroot, 
	      Real& coeffP, 
	      multi1d<Real>& resP, 
	      multi1d<Real>& rootQ, 
	      int& NEig, 
	      multi1d<Real>& EigValFunc,
	      const EigenConnectState& state) const;

    //! Construct stuff but use RatPolyDegPrec in the polynomial
    void initPrec(int& numroot, 
		  Real& coeffP, 
		  multi1d<Real>& resP, 
		  multi1d<Real>& rootQ, 
		  int& NEig, 
		  multi1d<Real>& EigValFunc,
		  const EigenConnectState& state) const;


  private:
    //!  Partial constructor not allowed
    OvlapPartFrac4DFermAct();
    // Assignment
    void operator=(const OvlapPartFrac4DFermAct& a) {}

  private:
    Handle< FermBC<T,P,Q> >           fbc;   // fermion bc
    Handle< CreateFermState<T,P,Q> >  cfs;   // fermion state creator
    // Auxilliary action used for kernel of operator
    Handle< UnprecWilsonTypeFermAct<T,P,Q> > Mact;   
    OvlapPartFrac4DFermActParams params;
  };

}
#endif
