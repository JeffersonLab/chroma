// -*- C++ -*-
// $Id: unprec_ovlap_contfrac5d_fermact_array_w.h,v 1.2 2004-09-29 21:48:34 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned extended-Overlap (5D) (Naryanan&Neuberger) action
 */

#ifndef __unprec_ovlap_contfrac5d_fermact_array_w_h__
#define __unprec_ovlap_contfrac5d_fermact_array_w_h__

#include "fermact.h"
#include "actions/ferm/fermacts/overlap_state.h"
#include "actions/ferm/linop/lgherm_w.h"
#include "io/overlap_state_info.h"

using namespace QDP;

namespace Chroma
{
  //! Name and registration
  namespace UnprecOvlapContFrac5DFermActArrayEnv
  {
    extern const std::string name;
    extern const bool registered;
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
  };


  // Reader/writers

  //! Read the Continued Fraction parameters 
  void read(XMLReader& xml, const string& path, UnprecOvlapContFrac5DFermActParams& param);

  //! Write the Continued Fraction parameters
  void write(XMLReader& xml, const string& path, const UnprecOvlapContFrac5DFermActParams& param);


  //! 5D continued fraction overlap action (Borici,Wenger, Edwards)

  /*!
   * \ingroup fermact
   *
   * This operator applies the extended version of the hermitian overlap operator
   *   Chi  =   ((1+Mass)/(1-Mass)*gamma_5 + B) . Psi
   *  where  B  is the continued fraction of the zolotarev approx. to eps(H(m))
   */
  class UnprecOvlapContFrac5DFermActArray : public UnprecWilsonTypeFermAct< multi1d<LatticeFermion> >
  {
  public:

    // Construct the action out of a parameter structure
    UnprecOvlapContFrac5DFermActArray(
			     Handle< FermBC< multi1d< LatticeFermion> > > fbc_,
			     const UnprecOvlapContFrac5DFermActParams& param);
  

    //! Copy constructor
    UnprecOvlapContFrac5DFermActArray(const UnprecOvlapContFrac5DFermActArray& a) : 
      fbc(a.fbc), S_aux(a.S_aux), params(a.params), N5(a.N5), isLastZeroP(a.isLastZeroP)  {};

    //! Assignment
    /* Writer screws this up -- I might get rid of that 
       UnprecOvlapContFrac5DFermActArray& operator=(const UnprecOvlapContFrac5DFermActArray& a) {
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
    Real quark_mass() const {return params.Mass;}

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
    ~UnprecOvlapContFrac5DFermActArray() {}


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


  protected:
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

  private:
    Handle< FermBC< multi1d<LatticeFermion> > >  fbc;
    Handle< UnprecWilsonTypeFermAct<LatticeFermion> > S_aux;
    UnprecOvlapContFrac5DFermActParams params;
    int  N5;
    bool isLastZeroP;
  };

}

#endif
