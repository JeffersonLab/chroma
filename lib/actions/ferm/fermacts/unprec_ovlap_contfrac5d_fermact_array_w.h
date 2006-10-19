// -*- C++ -*-
// $Id: unprec_ovlap_contfrac5d_fermact_array_w.h,v 3.4 2006-10-19 16:01:29 edwards Exp $
/*! \file
 *  \brief Unpreconditioned extended-Overlap (5D) (Naryanan&Neuberger) action
 */

#ifndef __unprec_ovlap_contfrac5d_fermact_array_w_h__
#define __unprec_ovlap_contfrac5d_fermact_array_w_h__

#include "unprec_wilstype_fermact_w.h"
#include "actions/ferm/fermstates/overlap_state.h"
#include "actions/ferm/linop/lgherm_w.h"
#include "actions/ferm/linop/lDeltaLs_w.h"
#include "io/enum_io/enum_coeffs_io.h"

namespace Chroma
{
  //! Name and registration
  namespace UnprecOvlapContFrac5DFermActArrayEnv
  {
    extern const std::string name;
    bool registerAll();
  }


  //! Params for 5D overlap ferm acts
  struct UnprecOvlapContFrac5DFermActParams
  {
    //! Default empty construction
    UnprecOvlapContFrac5DFermActParams() {};

    //! Read params from XML
    UnprecOvlapContFrac5DFermActParams(XMLReader& in, const std::string& path);
  
    
    Real Mass;     //!< Fermion Mass
    int RatPolyDeg; //!<  Degree of the Rational Poly
    CoeffType approximation_type;  //!< ZOLOTAREV | TANH | Other approximation coeffs
    std::string AuxFermAct;        //!<  The auxiliary ferm act
    std::string AuxFermActGrp;     //!<  The group name for the auxiliary fermion action
    Real ApproxMin;                //!< Approximate min eigenvalue of H_T
    Real ApproxMax;                //!< Approximate max eigenvalue of H_T
  };


  // Reader/writers

  //! Read the Continued Fraction parameters 
  void read(XMLReader& xml, const string& path, UnprecOvlapContFrac5DFermActParams& param);

  //! Write the Continued Fraction parameters
  void write(XMLWriter& xml, const string& path, const UnprecOvlapContFrac5DFermActParams& param);


  //! 5D continued fraction overlap action (Borici,Wenger, Edwards)

  /*!
   * \ingroup fermacts
   *
   * This operator applies the extended version of the hermitian overlap operator
   *   Chi  =   ((1+Mass)/(1-Mass)*gamma_5 + B) . Psi
   *  where  B  is the continued fraction of the zolotarev approx. to eps(H(m))
   */
  class UnprecOvlapContFrac5DFermActArray : public UnprecWilsonTypeFermAct5D<LatticeFermion, 
					    multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    // Construct the action out of a parameter structure
    UnprecOvlapContFrac5DFermActArray(Handle< FermBC<T,P,Q> > fbc_,
				      const UnprecOvlapContFrac5DFermActParams& param);
  
    //! 5D size
    int size(void) const { return N5; }

    //! Return the quark mass
    Real getQuarkMass() const {return params.Mass;}

    //! Produce a linear operator for this action
    UnprecLinearOperatorArray<T,P,Q>* linOp(Handle< FermState<T,P,Q> > state) const;

    //! Produce a Pauli-Villars linear operator for this action
    UnprecLinearOperatorArray<T,P,Q>* linOpPV(Handle< FermState<T,P,Q> > state) const;

    //! Produce a linear operator for this action
    LinearOperatorArray<T>* lnonHermLinOp(Handle< FermState<T,P,Q> > state) const;

    //! Produce a linear operator M^dag.M for this action
    LinearOperatorArray<T>* lnonHermMdagM(Handle< FermState<T,P,Q> > state) const;

    //! Matrix is itself hermitian so just return linOp here.
    LinearOperatorArray<T>* hermitianLinOp(Handle< FermState<T,P,Q> > state) const 
      {
	return linOp(state);
      }

    //! Produce an unpreconditioned linear operator projecting 5D to 4D (the inverse of qprop below)
    LinearOperator<T>* linOp4D(Handle< FermState<T,P,Q> > state,
			       const Real& m_q,
			       const GroupXML_t& invParam) const
    {
      QDPIO::cerr << "linOp4D not implemented" << endl;
      QDP_abort(1);
      return 0;
    }
    
    //! Produce a  DeltaLs = 1-epsilon^2(H) operator
    LinearOperator<T>* DeltaLs(Handle< FermState<T,P,Q> > state,
			       const GroupXML_t& invParam) const 
      {
	Handle< LinearOperator<T> >  lin(linOp4D(state,Real(0),invParam));
	return new lDeltaLs(lin);
      }

    //! Compute quark propagator over base type
    SystemSolver<LatticeFermion>* qprop(Handle< FermState<T,P,Q> > state,
					const GroupXML_t& invParam) const;

    //! Destructor is automatic
    ~UnprecOvlapContFrac5DFermActArray() {}


    // Create state functions

    //! Create OverlapConnectState from XML
    OverlapConnectState* 
    createState(const multi1d<LatticeColorMatrix>& u, 
		XMLReader& state_info_xml,
		const string& state_info_path) const;
    
    //! Given links, create the state needed for the linear operators
    /*! Override the parent */
    
    //! Create a ConnectState with just the gauge fields
    OverlapConnectState*
    createState(const multi1d<LatticeColorMatrix>& u_) const ;

    //! Create a ConnectState with just the gauge fields, and a lower
    //!  approximation bound
    OverlapConnectState*
    createState(const multi1d<LatticeColorMatrix>& u_, 
		const Real& approxMin_) const ;
 
    
    //! Create a connect State with just approximation range bounds
    OverlapConnectState*
    createState(const multi1d<LatticeColorMatrix>& u_, 
		const Real& approxMin_, 
		const Real& approxMax_) const;


    //! Create OverlapConnectState with eigenvalues/vectors 
    OverlapConnectState*
    createState(const multi1d<LatticeColorMatrix>& u_,
		const multi1d<Real>& lambda_lo_,
		const multi1d<LatticeFermion>& evecs_lo_, 
		const Real& lambda_hi_) const;


    //! Create from OverlapStateInfo Structure
    OverlapConnectState*
    createState(const multi1d<LatticeColorMatrix>& u_,
		const OverlapStateInfo& state_info) const;


  protected:
    //! Return the fermion create state for this action
    const CreateFermState<T,P,Q>& getCreateState() const {return *cfs;}

    //! Helper in construction
    void init(Real& scale_fac,
	      multi1d<Real>& alpha,
	      multi1d<Real>& beta,
	      int& NEig,
	      multi1d<Real>& EigValFunc,
	      const OverlapConnectState& state) const;
  private:
    // Hide partial constructor
    UnprecOvlapContFrac5DFermActArray();
    //! Hide =
    void operator=(const UnprecOvlapContFrac5DFermActArray& a) {}

  private:
    Handle< FermBC<T,P,Q> >           fbc;   // fermion bc
    Handle< CreateFermState<T,P,Q> >  cfs;   // fermion state creator
    Handle< UnprecWilsonTypeFermAct<T,P,Q> > S_aux;
    UnprecOvlapContFrac5DFermActParams params;
    int  N5;
    bool isLastZeroP;
  };

}

#endif
