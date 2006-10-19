// $Id: eoprec_zolo_nef_fermact_array_w.cc,v 3.1 2006-10-19 16:01:27 edwards Exp $
/*! \file
 *  \brief Unpreconditioned NEF fermion action
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/eoprec_zolo_nef_fermact_array_w.h"
#include "actions/ferm/fermacts/unprec_wilson_fermact_w.h"
#include "actions/ferm/linop/unprec_nef_linop_array_w.h"
#include "actions/ferm/linop/eoprec_nef_general_linop_array_w.h"

#include "actions/ferm/fermacts/zolotarev.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermstates/ferm_createstate_reader_w.h"

#include "actions/ferm/qprop/quarkprop4_w.h"
#include "actions/ferm/qprop/nef_quarkprop4_w.h"

namespace Chroma
{
  //! Hooks to register the class with the fermact factory
  namespace EvenOddPrecZoloNEFFermActArrayEnv
  {
    //! Callback function
    WilsonTypeFermAct5D<LatticeFermion,
			multi1d<LatticeColorMatrix>,
			multi1d<LatticeColorMatrix> >* createFermAct5D(XMLReader& xml_in,
								       const std::string& path)
    {
      return new EvenOddPrecZoloNEFFermActArray(CreateFermStateEnv::reader(xml_in, path), 
						EvenOddPrecZoloNEFFermActArrayParams(xml_in, path));
    }

    //! Callback function
    /*! Differs in return type */
    FermionAction<LatticeFermion,
		  multi1d<LatticeColorMatrix>,
		  multi1d<LatticeColorMatrix> >* createFermAct(XMLReader& xml_in,
							       const std::string& path)
    {
      return createFermAct5D(xml_in, path);
    }

    //! Name to be used
    const std::string name = "ZOLO_NEF";

    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::TheFermionActionFactory::Instance().registerObject(name, createFermAct);
	success &= Chroma::TheWilsonTypeFermAct5DFactory::Instance().registerObject(name, createFermAct5D);
	registered = true;
      }
      return success;
    }
  }


  //! Read parameters
  EvenOddPrecZoloNEFFermActArrayParams::EvenOddPrecZoloNEFFermActArrayParams(XMLReader& xml, 
									     const std::string& path)
  {
    XMLReader paramtop(xml, path);
    try {
      // Read the stuff for the action
      read(paramtop, "OverMass", OverMass);
      read(paramtop, "Mass", Mass);
      read(paramtop, "b5", b5);
      read(paramtop, "c5", c5);
      read(paramtop, "N5", N5);

      if( paramtop.count("ApproximationType") == 1 ) 
      { 
      	read(paramtop, "ApproximationType", approximation_type);
      }
      else 
      {
	// Default coeffs are unscaled tanh
	approximation_type = COEFF_TYPE_TANH_UNSCALED;
      }

      if (approximation_type == COEFF_TYPE_ZOLOTAREV)
      {
	read(paramtop, "ApproxMin", ApproxMin);
	read(paramtop, "ApproxMax", ApproxMax);
      }
      else
      {
	ApproxMin = ApproxMax = 0.0;
      }
    }
    catch(const string& e) { 
      QDPIO::cerr << "Caught Exception : " << e << endl;
    }
  }


  //! Read parameters
  void read(XMLReader& xml, const string& path, EvenOddPrecZoloNEFFermActArrayParams& param)
  {
    EvenOddPrecZoloNEFFermActArrayParams tmp(xml, path);
    param = tmp;
  }



  void EvenOddPrecZoloNEFFermActArray::initCoeffs(multi1d<Real>& b5_arr,
						  multi1d<Real>& c5_arr) const
  {
    b5_arr.resize(params.N5);
    c5_arr.resize(params.N5);

    Real approxMin = params.ApproxMin;
    Real approxMax = params.ApproxMax;

    Real epsilon;
    Real scale_fac;
    zolotarev_data *rdata;

    switch(params.approximation_type) 
    {
    case COEFF_TYPE_ZOLOTAREV:
      epsilon = approxMin / approxMax;
      QDPIO::cout << "Initing Linop with Zolotarev Coefficients: epsilon = " << epsilon << endl;
      rdata = zolotarev(toFloat(epsilon), params.N5, 0);
      scale_fac = approxMax;
      break;

    case COEFF_TYPE_TANH_UNSCALED:
      epsilon = approxMin;
      QDPIO::cout << "Initing Linop with Higham Rep tanh Coefficients" << endl;
      rdata = higham(toFloat(epsilon), params.N5);
      scale_fac = Real(1);
      break;

    default:
      // The map system should ensure that we never get here but 
      // just for style
      QDPIO::cerr << "Unknown coefficient type: " << params.approximation_type
		  << endl;
      QDP_abort(1);
    }

    if( rdata->n != params.N5 ) { 
      QDPIO::cerr << "Error:rdata->n != N5" << endl;
      QDP_abort(1);
    }

    multi1d<Real> gamma(params.N5);
    for(int i=0; i < params.N5; i++) { 
      gamma[i] = Real(rdata->gamma[i])*scale_fac;
    }

    Real maxerr=rdata->Delta;
    zolotarev_free(rdata);

    QDPIO::cout << "Initing Zolo NEF Linop: N5=" << params.N5 
		<< " b5+c5=" << params.b5+params.c5 
		<<  " b5-c5=" << params.b5-params.c5 
		<<  " epsilon=" << epsilon 
                <<  " scale=" << approxMax 
		<<  " maxerr=" << maxerr << endl;

    
    for(int i = 0; i < params.N5; i++) { 
      Real omega = Real(1)/gamma[i];

      b5_arr[i] = Real(0.5)*( (omega + Real(1))*params.b5 + (omega - Real(1))*params.c5 );
      c5_arr[i] = Real(0.5)*( (omega - Real(1))*params.b5 + (omega + Real(1))*params.c5 );
    
      QDPIO::cout << "i=" << i << " omega["<<i<<"]=" << omega
		  << " b5["<< i << "] ="<< b5_arr[i] 
		  << " c5["<< i << "] ="<< c5_arr[i] << endl;
      QDPIO::cout << "i=" << i 
		  << " b5["<<i<<"]+c5["<<i<<"]/omega["<<i<<"] =" 
		  << (b5_arr[i]+c5_arr[i])/omega
		  << " b5["<<i<<"]-c5["<<i<<"]=" << b5_arr[i]-c5_arr[i] << endl;
    }

  }


  //! Produce an even-odd preconditioned linear operator for this action with arbitrary quark mass
  EvenOddPrecDWLikeLinOpBaseArray<LatticeFermion,
				  multi1d<LatticeColorMatrix>,
				  multi1d<LatticeColorMatrix> >* 
  EvenOddPrecZoloNEFFermActArray::precLinOp(Handle< FermState<T,P,Q> > state,
					    const Real& m_q) const
  {
    multi1d<Real> b5_arr;
    multi1d<Real> c5_arr;

    initCoeffs(b5_arr,c5_arr);
    
    return new EvenOddPrecGenNEFDWLinOpArray(state,params.OverMass,b5_arr,c5_arr,m_q,params.N5);
  }


  //! Produce an unpreconditioned linear operator for this action with arbitrary quark mass
  UnprecDWLikeLinOpBaseArray<LatticeFermion,
			     multi1d<LatticeColorMatrix>,
			     multi1d<LatticeColorMatrix> >* 
  EvenOddPrecZoloNEFFermActArray::unprecLinOp(Handle< FermState<T,P,Q> > state,
					      const Real& m_q) const
  {
    multi1d<Real> b5_arr;
    multi1d<Real> c5_arr;

    initCoeffs(b5_arr,c5_arr);
    
    return new UnprecNEFDWLinOpArray(state,params.OverMass,b5_arr,c5_arr,m_q,params.N5);
  }


  // Given a complete propagator as a source, this does all the inversions needed
  void 
  EvenOddPrecZoloNEFFermActArray::quarkProp(LatticePropagator& q_sol, 
					    XMLWriter& xml_out,
					    const LatticePropagator& q_src,
					    int t_src, int j_decay,
					    Handle< FermState<T,P,Q> > state,
					    const GroupXML_t& invParam,
					    QuarkSpinType quarkSpinType,
					    bool obsvP,
					    int& ncg_had) const
  {
    if (obsvP && (quarkSpinType == QUARK_SPIN_TYPE_FULL))
      nef_quarkProp4(q_sol, xml_out, q_src, t_src, j_decay, *this, state, invParam, ncg_had);
    else
    {
      Handle< SystemSolver<LatticeFermion> > qprop(this->qprop(state,invParam));
      quarkProp4(q_sol, xml_out, q_src, qprop, quarkSpinType, ncg_had);
    }
  }

}
