// -*- C++ -*-
// $Id: prec_ovext_fermact_array_w.h,v 2.3 2006-03-21 19:14:36 edwards Exp $
/*! \file
 *  \brief Unpreconditioned extended-Overlap (5D) (Naryanan&Neuberger) action
 */

#ifndef __prec_ovext_fermact_array_w_h__
#define __prec_ovext_fermact_array_w_h__

#include "chromabase.h"
#include "fermact.h"
#include "fermbc.h"
#include "handle.h"
#include "linearop.h"
#include "prec_constdet_linop.h"
#include "state.h"
#include "actions/ferm/fermacts/overlap_state.h"
#include "io/enum_io/enum_coeffs_io.h"

#include "actions/ferm/fermacts/ovext_tuning_strategy.h"
namespace Chroma
{
  //! Name and registration
  namespace EvenOddPrecOvExtFermActArrayEnv
  {
    extern const std::string name;
    extern const bool registered;
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
  class EvenOddPrecOvExtFermActArray : public EvenOddPrecConstDetWilsonTypeFermAct5D<LatticeFermion, multi1d<LatticeColorMatrix> >
  {
  public:
    /*
    //! General FermBC
    EvenOddPrecOvExtFermActArray(Handle< FermBC< multi1d<LatticeFermion> > > fbc_, 
			    const Real& OverMass_,
			    const Real& b5_,
			    const Real& c5_,
			    const Real& Mass_, 
			    const int RatPolyDeg_,
			    const CoeffType approximation_type_) :  fbc(fbc_), param.OverMass(OverMass_),  param.b5(b5_), param.c5(c5_), param.Mass(Mass_),  param.RatPolyDeg(RatPolyDeg_), param.approximation_type(approximation_type_) {}
    */
    //! General FermBC
    EvenOddPrecOvExtFermActArray(Handle< FermBC< multi1d<LatticeFermion> > > fbc_, 
			    const EvenOddPrecOvExtFermActArrayParams& param_) :
      fbc(fbc_), param(param_) {

      // Set up the strategy for tuning the betas
      std::istringstream ts_is(param.tuning_strategy_xml);
      XMLReader tuning_xml(ts_is);
      std::string strategy_name;
      try { 
	read(tuning_xml, "/TuningStrategy/Name", strategy_name);
      }
      catch(const std::string& e) { 
	QDPIO::cerr << "Caught exception processing TuningStrategy: " << e << endl;
      }
      
      
      theTuningStrategy = TheAbsOvExtTuningStrategyFactory::Instance().createObject(strategy_name, tuning_xml, "/TuningStrategy");
      
    }


    //! Copy constructor
    EvenOddPrecOvExtFermActArray(const EvenOddPrecOvExtFermActArray& a) : 
      fbc(a.fbc), param(a.param), theTuningStrategy(a.theTuningStrategy) {}

    //! Assignment
    EvenOddPrecOvExtFermActArray& operator=(const EvenOddPrecOvExtFermActArray& a) {
      fbc=a.fbc; 
      param =a.param;
      theTuningStrategy=a.theTuningStrategy;

      return *this;
    }

    //! Return the fermion BC object for this action
    const FermBC< multi1d<LatticeFermion> >& getFermBC() const {return *fbc;}

    //! Length of DW flavor index/space
    int size() const {return getN5FromRatPolyDeg(param.RatPolyDeg);}

    //! Return the quark mass
    Real getQuarkMass() const {return param.Mass;}

    //! Produce a linear operator for this action
    const EvenOddPrecConstDetLinearOperator< multi1d<LatticeFermion>, multi1d<LatticeColorMatrix> >* linOp(Handle<const ConnectState> state) const;

    //! Produce a Pauli-Villars linear operator for this action
    const EvenOddPrecConstDetLinearOperator< multi1d<LatticeFermion>, multi1d<LatticeColorMatrix> >* linOpPV(Handle<const ConnectState> state) const
    {
      QDPIO::cerr << "Ovext::linOpPV not implemented" << endl;
      QDP_abort(1);
      return 0;
    }

    //! Produce a hermitian version of the linear operator
    const LinearOperator< multi1d<LatticeFermion> >* hermitianLinOp(Handle<const ConnectState> state) const
      {
	QDPIO::cerr << "EvenOddPrecOvExtFermActArray::gamma5HermLinOp not implemented" << endl;
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
      Handle< const LinearOperator<LatticeFermion> >  lin(linOp4D(state,Real(0),invParam));
      QDPIO::cout << "NOt yet implemented" << endl;
      QDP_abort(1);
      return 0x0;
    }


    //! Compute quark propagator over base type
    const SystemSolver<LatticeFermion>* qprop(Handle<const ConnectState> state,
					      const InvertParam_t& invParam) const;

    //! Destructor is automatic
    ~EvenOddPrecOvExtFermActArray() {}

  private:
    // Hide partial constructor
    EvenOddPrecOvExtFermActArray() {}
    int getN5FromRatPolyDeg(const int& RatPolyDeg) const {

      // Type 0 and Tanh approximations: 

      // If RatPolyDeg is even: => 2*(RatPolyDeg/2) + 1 = RatPolyDeg+1
      // If RatPolyDeg is odd: =>  2*((RatPolyDeg-1)/2 + 1 = RatPolyDeg
      if( RatPolyDeg % 2 == 0 ) { 
	return RatPolyDeg+1; 
      }
      else { 
	return RatPolyDeg;
      }
    }

    void init(int& Npoles, 
	      Real& coeffP, 
	      multi1d<Real>& resP,
	      multi1d<Real>& rootQ) const;

  private:
    Handle< FermBC< multi1d<LatticeFermion> > >  fbc;
    Handle< AbsOvExtTuningStrategy > theTuningStrategy;
    EvenOddPrecOvExtFermActArrayParams param;
  };

}



#endif
