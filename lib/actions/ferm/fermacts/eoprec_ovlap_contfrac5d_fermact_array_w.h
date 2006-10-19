// -*- C++ -*-
// $Id: eoprec_ovlap_contfrac5d_fermact_array_w.h,v 3.1 2006-10-19 16:01:27 edwards Exp $
/*! \file
 *  \brief Even-odd preconditioned Continued Fraction 5D
 */

#ifndef __prec_ovlap_contfrac5d_fermact_array_w_h__
#define __prec_ovlap_contfrac5d_fermact_array_w_h__

#include "eoprec_constdet_wilstype_fermact_w.h"
#include "actions/ferm/linop/lgherm_w.h"
#include "actions/ferm/linop/lDeltaLs_w.h"
#include "io/enum_io/enum_coeffs_io.h"

namespace Chroma
{
  //! Name and registration
  namespace EvenOddPrecOvlapContFrac5DFermActArrayEnv
  {
    extern const std::string name;
    bool registerAll();
  }


  //! Params for 5D overlap ferm acts
  struct EvenOddPrecOvlapContFrac5DFermActParams
  {
    //! Default empty construction
    EvenOddPrecOvlapContFrac5DFermActParams() {};

    //! Read params from XML
    EvenOddPrecOvlapContFrac5DFermActParams(XMLReader& in, const std::string& path);
  
    
    Real Mass;     //!< Fermion Mass

    int RatPolyDeg; //!<  Degree of the Rational Poly
    CoeffType approximation_type;  //!< ZOLOTAREV | TANH | Other approximation coeffs
    // AuxFermAct is Gone. Use only Unprec Wilson (and the linOp creates its own Dslash
    Real OverMass;  //!< Mass of auxiliary Wilson action
    Real ApproxMin; //!< Approximate min eigenvalue of H_T
    Real ApproxMax; //!< Approximate max eigenvalue of H_T
  };


  // Reader/writers

  //! Read the Continued Fraction parameters 
  void read(XMLReader& xml, const string& path, EvenOddPrecOvlapContFrac5DFermActParams& param);

  //! Write the Continued Fraction parameters
  void write(XMLWriter& xml, const string& path, const EvenOddPrecOvlapContFrac5DFermActParams& param);


  //! 5D continued fraction overlap action (Borici,Wenger, Edwards)

  /*!
   * \ingroup fermacts
   *
   * This operator applies the extended version of the hermitian overlap operator
   *   Chi  =   ((1+Mass)/(1-Mass)*gamma_5 + B) . Psi
   *  where  B  is the continued fraction of the zolotarev approx. to eps(H(m))
   */
  class EvenOddPrecOvlapContFrac5DFermActArray : public EvenOddPrecConstDetWilsonTypeFermAct5D<LatticeFermion, 
						 multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Construct the action out of a parameter structure
    EvenOddPrecOvlapContFrac5DFermActArray(
      Handle< CreateFermState<T,P,Q> > cfs_,
      const EvenOddPrecOvlapContFrac5DFermActParams& param);
  
    //! Copy constructor
    EvenOddPrecOvlapContFrac5DFermActArray(const EvenOddPrecOvlapContFrac5DFermActArray& a) : 
      cfs(a.cfs), params(a.params), N5(a.N5), isLastZeroP(a.isLastZeroP)  {};

    int size(void) const { return N5; }

    //! Return the quark mass
    Real getQuarkMass() const {return params.Mass;}

    //! Produce a linear operator for this action
    EvenOddPrecConstDetLinearOperatorArray<T,P,Q>* linOp(Handle< FermState<T,P,Q> > state) const;

    //! Produce a Pauli-Villars linear operator for this action
    EvenOddPrecConstDetLinearOperatorArray<T,P,Q>* linOpPV(Handle< FermState<T,P,Q> > state) const;

    //! produce a hermitian version of this operator
    //  but it is already hermitian
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
    ~EvenOddPrecOvlapContFrac5DFermActArray() {}


  protected:
    //! Return the fermion BC object for this action
    const CreateFermState<T,P,Q>& getCreateState() const {return *cfs;}

    //! Helper in construction
    void init(Real& scale_fac,
	      multi1d<Real>& alpha,
	      multi1d<Real>& beta) const;

  private:
    // Hide partial constructor
    EvenOddPrecOvlapContFrac5DFermActArray();

  private:
    Handle< CreateFermState<T,P,Q> >  cfs;
    EvenOddPrecOvlapContFrac5DFermActParams params;
    int  N5;
    bool isLastZeroP;
  };

}

#endif
