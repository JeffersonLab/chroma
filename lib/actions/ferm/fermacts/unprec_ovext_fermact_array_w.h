// -*- C++ -*-
// $Id: unprec_ovext_fermact_array_w.h,v 3.4 2006-10-19 16:01:29 edwards Exp $
/*! \file
 *  \brief Unpreconditioned extended-Overlap (5D) (Naryanan&Neuberger) action
 */

#ifndef __unprec_ovext_fermact_array_w_h__
#define __unprec_ovext_fermact_array_w_h__

#include "unprec_wilstype_fermact_w.h"
#include "handle.h"
#include "actions/ferm/fermstates/overlap_state.h"
#include "actions/ferm/linop/unprec_wilson_linop_w.h"
#include "io/enum_io/enum_coeffs_io.h"

#include "actions/ferm/fermacts/ovext_tuning_strategy.h"

namespace Chroma
{
  //! Name and registration
  namespace UnprecOvExtFermActArrayEnv
  {
    extern const std::string name;
    bool registerAll();
  }
  

  //! Params for NEFF
  struct UnprecOvExtFermActArrayParams
  {
    UnprecOvExtFermActArrayParams() {}
    UnprecOvExtFermActArrayParams(XMLReader& in, const std::string& path);
    
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
  void read(XMLReader& xml, const string& path, UnprecOvExtFermActArrayParams& param);
  void write(XMLWriter& xml, const string& path, const UnprecOvExtFermActArrayParams& param);


  //! Unpreconditioned Extended-Overlap (N&N) linear operator
  /*!
   * \ingroup fermacts
   *
   * This operator applies the extended version of the hermitian overlap operator
   *   Chi  =   ((1+Mass)/(1-Mass)*gamma_5 + B) . Psi
   *  where  B  is the continued fraction of the pole approx. to eps(H(m))
   */
  class UnprecOvExtFermActArray : public UnprecWilsonTypeFermAct5D<LatticeFermion, 
				  multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! General FermBC
    UnprecOvExtFermActArray(Handle< CreateFermState<T,P,Q> > cfs_, 
			    const UnprecOvExtFermActArrayParams& param_);

    //! Copy constructor
    UnprecOvExtFermActArray(const UnprecOvExtFermActArray& a) : 
      cfs(a.cfs), param(a.param), theTuningStrategy(a.theTuningStrategy) {}

    //! Length of DW flavor index/space
    int size() const {return getN5FromRatPolyDeg(param.RatPolyDeg);}

    //! Return the quark mass
    Real getQuarkMass() const {return param.Mass;}

    //! Produce a linear operator for this action
    UnprecLinearOperatorArray<T,P,Q>* linOp(Handle< FermState<T,P,Q> > state) const;

    //! Produce a Pauli-Villars linear operator for this action
    UnprecLinearOperatorArray<T,P,Q>* linOpPV(Handle< FermState<T,P,Q> > state) const
    {
      QDPIO::cerr << "Ovext::linOpPV not implemented" << endl;
      QDP_abort(1);
      return 0;
    }

    //! Produce a hermitian version of the linear operator
    LinearOperatorArray<T>* hermitianLinOp(Handle< FermState<T,P,Q> > state) const
      {
	QDPIO::cerr << "UnprecOvExtFermActArray::gamma5HermLinOp not implemented" << endl;
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
      Handle< const LinearOperator<T> >  lin(linOp4D(state,Real(0),invParam));
      QDPIO::cout << "NOt yet implemented" << endl;
      QDP_abort(1);
      return 0x0;
    }


    //! Compute quark propagator over base type
    SystemSolver<LatticeFermion>* qprop(Handle< FermState<T,P,Q> > state,
					const GroupXML_t& invParam) const;
    
    //! Destructor is automatic
    ~UnprecOvExtFermActArray() {}

  protected:
    //! Return the fermion create state for this action
    const CreateFermState<T,P,Q>& getCreateState() const {return *cfs;}

  private:
    //! Initializer
    int getN5FromRatPolyDeg(const int& RatPolyDeg) const;

    void init(int& Npoles, 
	      Real& coeffP, 
	      multi1d<Real>& resP,
	      multi1d<Real>& rootQ) const;

  private:
    // Hide partial constructor
    UnprecOvExtFermActArray() {}
    //! Hide =
    void operator=(const UnprecOvExtFermActArray& a) {}

  private:
    Handle< CreateFermState<T,P,Q> >  cfs;
    Handle< AbsOvExtTuningStrategy > theTuningStrategy;

    UnprecOvExtFermActArrayParams param;
  };

}



#endif
