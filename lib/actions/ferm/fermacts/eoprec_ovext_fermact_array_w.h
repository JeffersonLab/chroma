// -*- C++ -*-
// $Id: eoprec_ovext_fermact_array_w.h,v 3.1 2006-10-19 16:01:27 edwards Exp $
/*! \file
 *  \brief Unpreconditioned extended-Overlap (5D) (Naryanan&Neuberger) action
 */

#ifndef __prec_ovext_fermact_array_w_h__
#define __prec_ovext_fermact_array_w_h__

#include "eoprec_constdet_wilstype_fermact_w.h"
#include "handle.h"
#include "eoprec_constdet_linop.h"
#include "io/enum_io/enum_coeffs_io.h"

#include "actions/ferm/fermacts/ovext_tuning_strategy.h"

namespace Chroma
{
  //! Name and registration
  namespace EvenOddPrecOvExtFermActArrayEnv
  {
    extern const std::string name;
    bool registerAll();
  }
  

  //! Params for NEFF
  struct EvenOddPrecOvExtFermActArrayParams
  {
    EvenOddPrecOvExtFermActArrayParams() {}
    EvenOddPrecOvExtFermActArrayParams(XMLReader& in, const std::string& path);
    
    Real OverMass;
    Real b5;
    Real c5;
    Real Mass;
    int  RatPolyDeg;
    Real ApproxMin;
    Real ApproxMax;
    CoeffType approximation_type;
    std::string tuning_strategy_xml;
  };


  // Reader/writers
  void read(XMLReader& xml, const string& path, EvenOddPrecOvExtFermActArrayParams& param);
  void write(XMLWriter& xml, const string& path, const EvenOddPrecOvExtFermActArrayParams& param);


  //! EvenOddPreconditioned Extended-Overlap (N&N) linear operator
  /*!
   * \ingroup fermacts
   *
   * This operator applies the extended version of the hermitian overlap operator
   *   Chi  =   ((1+Mass)/(1-Mass)*gamma_5 + B) . Psi
   *  where  B  is the continued fraction of the pole approx. to eps(H(m))
   */
  class EvenOddPrecOvExtFermActArray : public EvenOddPrecConstDetWilsonTypeFermAct5D<LatticeFermion, 
				       multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! General FermBC
    EvenOddPrecOvExtFermActArray(Handle< CreateFermState<T,P,Q> > cfs_, 
				 const EvenOddPrecOvExtFermActArrayParams& param_);

    //! Copy constructor
    EvenOddPrecOvExtFermActArray(const EvenOddPrecOvExtFermActArray& a) : 
      cfs(a.cfs), param(a.param), theTuningStrategy(a.theTuningStrategy) {}

    //! Length of DW flavor index/space
    int size() const {return getN5FromRatPolyDeg(param.RatPolyDeg);}

    //! Return the quark mass
    Real getQuarkMass() const {return param.Mass;}

    //! Produce a linear operator for this action
    EvenOddPrecConstDetLinearOperatorArray<T,P,Q>* linOp(Handle< FermState<T,P,Q> > state) const;

    //! Produce a Pauli-Villars linear operator for this action
    EvenOddPrecConstDetLinearOperatorArray<T,P,Q>* linOpPV(Handle< FermState<T,P,Q> > state) const
      {
	QDPIO::cerr << "Ovext::linOpPV not implemented" << endl;
	QDP_abort(1);
	return 0;
      }

    //! Produce a hermitian version of the linear operator
    LinearOperatorArray<T>* hermitianLinOp(Handle< FermState<T,P,Q> > state) const
      {
	QDPIO::cerr << "EvenOddPrecOvExtFermActArray::gamma5HermLinOp not implemented" << endl;
	QDP_abort(1);
	return 0;
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
	QDPIO::cout << "NOt yet implemented" << endl;
	QDP_abort(1);
	return 0x0;
      }


    //! Compute quark propagator over base type
    SystemSolver<LatticeFermion>* qprop(Handle< FermState<T,P,Q> > state,
					const GroupXML_t& invParam) const;

    //! Destructor is automatic
    ~EvenOddPrecOvExtFermActArray() {}

  protected:
    //! Return the fermion BC object for this action
    const CreateFermState<T,P,Q>& getCreateState() const {return *cfs;}

    //! Default constructor
    EvenOddPrecOvExtFermActArray() {}
    //! Hide =
    void operator=(const EvenOddPrecOvExtFermActArray& a) {}

  private:
    //! Part of initializer
    int getN5FromRatPolyDeg(const int& RatPolyDeg) const;

    void init(int& Npoles, 
	      Real& coeffP, 
	      multi1d<Real>& resP,
	      multi1d<Real>& rootQ) const;

  private:
    Handle< CreateFermState<T,P,Q> >  cfs;
    Handle< AbsOvExtTuningStrategy > theTuningStrategy;
    EvenOddPrecOvExtFermActArrayParams param;
  };

}



#endif
