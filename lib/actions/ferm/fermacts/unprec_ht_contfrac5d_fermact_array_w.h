// -*- C++ -*-
// $Id: unprec_ht_contfrac5d_fermact_array_w.h,v 3.3 2006-10-19 16:01:29 edwards Exp $
/*! \file
 *  \brief Unpreconditioned H_T kernel continued fraction (5D) action
 */

#ifndef __unprec_ht_contfrac5d_fermact_array_w_h__
#define __unprec_ht_contfrac5d_fermact_array_w_h__

#include "unprec_wilstype_fermact_w.h"
#include "io/enum_io/enum_coeffs_io.h"


namespace Chroma
{
  //! Name and registration
  namespace UnprecHTContFrac5DFermActArrayEnv
  {
    extern const std::string name;
    bool registerAll();
  }


  //! Params for 5D overlap ferm acts
  struct UnprecHTContFrac5DFermActParams
  {
    //! Default empty construction
    UnprecHTContFrac5DFermActParams() {};

    //! Read params from XML
    UnprecHTContFrac5DFermActParams(XMLReader& in, const std::string& path);
  
    Real OverMass;  //!<  Wilson-operator mass
    Real Mass;      //!<  Fermion Mass
    Real b5;        //!<  Mobius b5
    Real c5;        //!<  Mobius c5
    Real ApproxMin; //!<  Approximate min eigenvalue of H_T
    Real ApproxMax; //!<  Approximate max eigenvalue of H_T
    int RatPolyDeg; //!<  Degree of the Rational Poly
    CoeffType approximation_type;  //!< ZOLOTAREV | TANH | Other approximation coeffs
  };


  // Reader/writers

  //! Read the Continued Fraction parameters 
  void read(XMLReader& xml, const string& path, UnprecHTContFrac5DFermActParams& param);

  //! Write the Continued Fraction parameters
  void write(XMLWriter& xml, const string& path, const UnprecHTContFrac5DFermActParams& param);


  //! 5D continued fraction overlap action using H_T kernel
  /*!
   * \ingroup fermacts
   *
   * This operator applies the extended version of the hermitian overlap operator
   *   Chi  =   ((1+Mass)/(1-Mass)*gamma_5 + B) . Psi
   *  where  B  is the continued fraction of the zolotarev approx. to eps(H_T(m))
   */
  class UnprecHTContFrac5DFermActArray : public UnprecWilsonTypeFermAct5D<LatticeFermion, 
					 multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    // Construct the action out of a parameter structure
    UnprecHTContFrac5DFermActArray(Handle< CreateFermState<T,P,Q> > cfs_,
				   const UnprecHTContFrac5DFermActParams& param);
  
    //! Copy constructor
    UnprecHTContFrac5DFermActArray(const UnprecHTContFrac5DFermActArray& a) : 
      cfs(a.cfs), params(a.params), N5(a.N5), isLastZeroP(a.isLastZeroP)  {};

    int size(void) const { return N5; }

    //! Return the quark mass
    Real getQuarkMass() const {return params.Mass;}

    //! Produce a linear operator for this action
    UnprecLinearOperatorArray<T,P,Q>* linOp(Handle< FermState<T,P,Q> > state) const;

    //! Produce a Pauli-Villars linear operator for this action
    UnprecLinearOperatorArray<T,P,Q>* linOpPV(Handle< FermState<T,P,Q> > state) const
    {
      QDPIO::cerr << "UnprecHTCFZ::linOpPV not implemented" << endl;
      QDP_abort(1);
      return 0;
    }

    //! produce gamma_5 times M 
    LinearOperatorArray<T>* hermitianLinOp(Handle< FermState<T,P,Q> > state) const
    {
      // The matrix itself is apparently hermitian
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
	QDPIO::cerr << "DeltaLs not implemented" << endl;
	QDP_abort(1);
	return 0;
      }

    //! Compute quark propagator over base type
    SystemSolver<LatticeFermion>* qprop(Handle< FermState<T,P,Q> > state,
					const GroupXML_t& invParam) const;

    //! Destructor is automatic
    ~UnprecHTContFrac5DFermActArray() {}

  protected:
    //! Return the fermion BC object for this action
    const CreateFermState<T,P,Q>& getCreateState() const {return *cfs;}

    //! Helper in construction
    void init(multi1d<Real>& alpha,
	      multi1d<Real>& beta) const;

  private:
    // Hide partial constructor
    UnprecHTContFrac5DFermActArray();
    void operator=(const UnprecHTContFrac5DFermActArray& a) {}

  private:
    Handle< CreateFermState<T,P,Q> >  cfs;
    UnprecHTContFrac5DFermActParams params;
    int  N5;
    bool isLastZeroP;
  };

}

#endif
