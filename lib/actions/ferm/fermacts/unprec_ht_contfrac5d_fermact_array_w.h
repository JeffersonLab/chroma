// -*- C++ -*-
// $Id: unprec_ht_contfrac5d_fermact_array_w.h,v 1.4 2005-01-14 20:13:04 edwards Exp $
/*! \file
 *  \brief Unpreconditioned H_T kernel continued fraction (5D) action
 */

#ifndef __unprec_ht_contfrac5d_fermact_array_w_h__
#define __unprec_ht_contfrac5d_fermact_array_w_h__

#include "fermact.h"


namespace Chroma
{
  //! Name and registration
  namespace UnprecHTContFrac5DFermActArrayEnv
  {
    extern const std::string name;
    extern const bool registered;
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
   * \ingroup fermact
   *
   * This operator applies the extended version of the hermitian overlap operator
   *   Chi  =   ((1+Mass)/(1-Mass)*gamma_5 + B) . Psi
   *  where  B  is the continued fraction of the zolotarev approx. to eps(H_T(m))
   */
  class UnprecHTContFrac5DFermActArray : public UnprecWilsonTypeFermAct5D<LatticeFermion, multi1d<LatticeColorMatrix> >
  {
  public:

    // Construct the action out of a parameter structure
    UnprecHTContFrac5DFermActArray(
			     Handle< FermBC< multi1d< LatticeFermion> > > fbc_,
			     const UnprecHTContFrac5DFermActParams& param);
  

    //! Copy constructor
    UnprecHTContFrac5DFermActArray(const UnprecHTContFrac5DFermActArray& a) : 
      fbc(a.fbc), params(a.params), N5(a.N5), isLastZeroP(a.isLastZeroP)  {};

    //! Return the fermion BC object for this action
    const FermBC< multi1d< LatticeFermion> >& getFermBC() const {return *fbc;}

    int size(void) const { return N5; }

    //! Return the quark mass
    Real getQuarkMass() const {return params.Mass;}

    //! Produce a linear operator for this action
    const UnprecLinearOperator< multi1d<LatticeFermion>, multi1d<LatticeColorMatrix> >* linOp(Handle<const ConnectState> state) const;

    //! Produce a Pauli-Villars linear operator for this action
    const UnprecLinearOperator< multi1d<LatticeFermion>, multi1d<LatticeColorMatrix> >* linOpPV(Handle<const ConnectState> state) const
    {
      QDPIO::cerr << "UnprecHTCFZ::linOpPV not implemented" << endl;
      QDP_abort(1);
      return 0;
    }

    //! Produce a linear operator M^dag.M for this action
    const LinearOperator< multi1d<LatticeFermion> >* lMdagM(Handle<const ConnectState> state) const;

    //! produce gamma_5 times M 
    const LinearOperator< multi1d<LatticeFermion> >* gamma5HermLinOp(Handle<const ConnectState> state) const
    {
      QDPIO::cerr << "UnprecHTCFZ::linOpPV not implemented" << endl;
      QDP_abort(1);
      return 0;
    }

    //! Produce an unpreconditioned linear operator projecting 5D to 4D (the inverse of qprop below)
    const LinearOperator<LatticeFermion>* linOp4D(Handle<const ConnectState> state,
						  const Real& m_q,
						  const InvertParam_t& invParam) const
    {
      QDPIO::cerr << "linOp4D not implemented" << endl;
      QDP_abort(1);
      return 0;
    }
    
    //! Produce a  DeltaLs = 1-epsilon^2(H) operator
    const LinearOperator<LatticeFermion>* DeltaLs(Handle< const ConnectState> state,
						  const InvertParam_t& invParam) const 
    {
      QDPIO::cerr << "DeltaLs not implemented" << endl;
      QDP_abort(1);
      return 0;
    }

    //! Compute quark propagator over base type
    const SystemSolver<LatticeFermion>* qprop(Handle<const ConnectState> state,
					      const InvertParam_t& invParam) const;

    //! Destructor is automatic
    ~UnprecHTContFrac5DFermActArray() {}

  protected:
    //! Helper in construction
    void init(multi1d<Real>& alpha,
	      multi1d<Real>& beta) const;

  private:
    // Hide partial constructor
    UnprecHTContFrac5DFermActArray();
    void operator=(const UnprecHTContFrac5DFermActArray& a) {}

  private:
    Handle< FermBC< multi1d<LatticeFermion> > >  fbc;
    UnprecHTContFrac5DFermActParams params;
    int  N5;
    bool isLastZeroP;
  };

}

#endif
