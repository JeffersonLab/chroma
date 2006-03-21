// -*- C++ -*-
// $Id: prec_ovlap_contfrac5d_fermact_array_w.h,v 2.3 2006-03-21 04:42:49 edwards Exp $
/*! \file
 *  \brief Even-odd preconditioned Continued Fraction 5D
 */

#ifndef __prec_ovlap_contfrac5d_fermact_array_w_h__
#define __prec_ovlap_contfrac5d_fermact_array_w_h__

#include "fermact.h"
#include "actions/ferm/linop/lgherm_w.h"
#include "actions/ferm/linop/lDeltaLs_w.h"
#include "io/enum_io/enum_coeffs_io.h"

namespace Chroma
{
  //! Name and registration
  namespace EvenOddPrecOvlapContFrac5DFermActArrayEnv
  {
    extern const std::string name;
    extern const bool registered;
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
  class EvenOddPrecOvlapContFrac5DFermActArray : public EvenOddPrecConstDetWilsonTypeFermAct5D<LatticeFermion, multi1d<LatticeColorMatrix> >
  {
  public:
    //! Construct the action out of a parameter structure
    EvenOddPrecOvlapContFrac5DFermActArray(
      Handle< FermBC< multi1d< LatticeFermion> > > fbc_,
      const EvenOddPrecOvlapContFrac5DFermActParams& param);
  
    //! Copy constructor
    EvenOddPrecOvlapContFrac5DFermActArray(const EvenOddPrecOvlapContFrac5DFermActArray& a) : 
      fbc(a.fbc), params(a.params), N5(a.N5), isLastZeroP(a.isLastZeroP)  {};

    //! Return the fermion BC object for this action
    const FermBC< multi1d< LatticeFermion> >& getFermBC() const {return *fbc;}

    int size(void) const { return N5; }

    //! Return the quark mass
    Real getQuarkMass() const {return params.Mass;}

    //! Produce a linear operator for this action
    const EvenOddPrecConstDetLinearOperator< multi1d<LatticeFermion>, multi1d<LatticeColorMatrix> >* linOp(Handle<const ConnectState> state) const;

    //! Produce a Pauli-Villars linear operator for this action
    const EvenOddPrecConstDetLinearOperator< multi1d<LatticeFermion>, multi1d<LatticeColorMatrix> >* linOpPV(Handle<const ConnectState> state) const;

    //! produce a hermitian version of this operator
    //  but it is already hermitian
    const LinearOperator< multi1d<LatticeFermion> >* hermitianLinOp(Handle<const ConnectState> state) const {
      return linOp(state);
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
      Handle< const LinearOperator<LatticeFermion> >  lin(linOp4D(state,Real(0),invParam));
      return new lDeltaLs(lin);
    }

    //! Compute quark propagator over base type
    const SystemSolver<LatticeFermion>* qprop(Handle<const ConnectState> state,
					      const InvertParam_t& invParam) const;

    //! Destructor is automatic
    ~EvenOddPrecOvlapContFrac5DFermActArray() {}


//  protected:
    //! Helper in construction
    void init(Real& scale_fac,
	      multi1d<Real>& alpha,
	      multi1d<Real>& beta) const;

  private:
    // Hide partial constructor
    EvenOddPrecOvlapContFrac5DFermActArray();

  private:
    Handle< FermBC< multi1d<LatticeFermion> > >  fbc;
    EvenOddPrecOvlapContFrac5DFermActParams params;
    int  N5;
    bool isLastZeroP;
  };

}

#endif
