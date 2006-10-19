// -*- C++ -*-
// $Id: eoprec_ht_contfrac5d_fermact_array_w.h,v 3.1 2006-10-19 16:01:27 edwards Exp $
/*! \file
 *  \brief Even-odd preconditioned Continued Fraction 5D
 */

#ifndef __prec_ht_contfrac5d_fermact_array_w_h__
#define __prec_ht_contfrac5d_fermact_array_w_h__

#include "eoprec_constdet_wilstype_fermact_w.h"
#include "actions/ferm/linop/lgherm_w.h"
#include "actions/ferm/linop/lDeltaLs_w.h"
#include "io/enum_io/enum_coeffs_io.h"


namespace Chroma
{
  //! Name and registration
  namespace EvenOddPrecHtContFrac5DFermActArrayEnv
  {
    extern const std::string name;
    bool registerAll();
  }


  //! Params for 5D overlap ferm acts
  struct EvenOddPrecHtContFrac5DFermActParams
  {
    //! Default empty construction
    EvenOddPrecHtContFrac5DFermActParams() {};

    //! Read params from XML
    EvenOddPrecHtContFrac5DFermActParams(XMLReader& in, const std::string& path);
  
    Real Mass;     //!< Fermion Mass

    int RatPolyDeg; //!<  Degree of the Rational Poly
    CoeffType approximation_type;  //!< ZOLOTAREV | TANH | Other approximation coeffs
    // AuxFermAct is Gone. Use only Unprec Wilson (and the linOp creates its own Dslash
    Real OverMass;  //!< Mass of auxiliary Wilson action
    Real ApproxMin; //!< Approximate min eigenvalue of H_T
    Real ApproxMax; //!< Approximate max eigenvalue of H_T
    Real b5;        //!< b5 Moebius parameter
    Real c5;        //!< c5 Moebius parameter
  };


  // Reader/writers

  //! Read the Continued Fraction parameters 
  void read(XMLReader& xml, const string& path, EvenOddPrecHtContFrac5DFermActParams& param);

  //! Write the Continued Fraction parameters
  void write(XMLWriter& xml, const string& path, const EvenOddPrecHtContFrac5DFermActParams& param);


  //! 5D continued fraction overlap action (Borici,Wenger, Edwards)

  /*!
   * \ingroup fermacts
   *
   * This operator applies the extended version of the hermitian overlap operator
   *   Chi  =   ((1+Mass)/(1-Mass)*gamma_5 + B) . Psi
   *  where  B  is the continued fraction of the zolotarev approx. to eps(H(m))
   */
  class EvenOddPrecHtContFrac5DFermActArray : public EvenOddPrecConstDetWilsonTypeFermAct5D<LatticeFermion, 
					      multi1d<LatticeColorMatrix>,
					      multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    // Construct the action out of a parameter structure
    EvenOddPrecHtContFrac5DFermActArray(Handle< CreateFermState<T,P,Q> > fs_,
					const EvenOddPrecHtContFrac5DFermActParams& param);
  

    //! Copy constructor
    EvenOddPrecHtContFrac5DFermActArray(const EvenOddPrecHtContFrac5DFermActArray& a) : 
      fs(a.fs), params(a.params), N5(a.N5), isLastZeroP(a.isLastZeroP)  {};

    int size(void) const { return N5; }

    //! Return the quark mass
    Real getQuarkMass() const {return params.Mass;}

    //! Produce a linear operator for this action
    EvenOddPrecConstDetLinearOperatorArray<T,P,Q>* linOp(Handle< FermState<T,P,Q> > state) const;

    //! Produce a Pauli-Villars linear operator for this action
    EvenOddPrecConstDetLinearOperatorArray<T,P,Q>* linOpPV(Handle< FermState<T,P,Q> > state) const;

    //! produce hermitian version of linOp 
    LinearOperatorArray<T>* hermitianLinOp(Handle< FermState<T,P,Q> > state) const 
      {
	QDPIO::cerr << "Hermitian version of this operator is not yet implemented" << endl << flush;
	QDP_abort(1);
	
	return 0;
      }

    //! Produce an unpreconditioned linear operator projecting 5D to 4D (the inverse of qprop below)
    LinearOperator<LatticeFermion>* linOp4D(Handle< FermState<T,P,Q> > state,
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
    ~EvenOddPrecHtContFrac5DFermActArray() {}

  protected:
    //! Return the fermion BC object for this action
    const CreateFermState<T,P,Q>& getCreateState() const {return *fs;}

    //! Helper in construction
    void init(Real& scale_fac,
	      multi1d<Real>& alpha,
	      multi1d<Real>& beta) const;

  private:
    // Hide partial constructor
    EvenOddPrecHtContFrac5DFermActArray();

  private:
    Handle< CreateFermState<T,P,Q> >  fs;
    EvenOddPrecHtContFrac5DFermActParams params;
    int  N5;
    bool isLastZeroP;
  };

}

#endif
