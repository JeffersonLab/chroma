// $Id: prec_zolo_nef_fermact_array_w.cc,v 1.10 2004-12-24 04:23:20 edwards Exp $
/*! \file
 *  \brief Unpreconditioned NEF fermion action
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/prec_zolo_nef_fermact_array_w.h"
#include "actions/ferm/fermacts/unprec_wilson_fermact_w.h"
#include "actions/ferm/linop/unprec_nef_linop_array_w.h"
#include "actions/ferm/linop/prec_nef_general_linop_array_w.h"

#include "actions/ferm/fermacts/zolotarev.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermbcs/fermbcs_w.h"

using namespace Chroma;

namespace Chroma
{
  //! Hooks to register the class with the fermact factory
  namespace EvenOddPrecZoloNEFFermActArrayEnv
  {
    //! Callback function
    WilsonTypeFermAct< multi1d<LatticeFermion> >* createFermAct(XMLReader& xml_in,
								const std::string& path)
    {
      return new EvenOddPrecZoloNEFFermActArray(WilsonTypeFermBCArrayEnv::reader(xml_in, path), 
						EvenOddPrecZoloNEFFermActArrayParams(xml_in, path));
    }

    //! Callback function
    /*! Differs in return type */
    EvenOddPrecDWFermActBaseArray<LatticeFermion>* createDWFermAct(XMLReader& xml_in,
								   const std::string& path)
    {
      return new EvenOddPrecZoloNEFFermActArray(WilsonTypeFermBCArrayEnv::reader(xml_in, path), 
						EvenOddPrecZoloNEFFermActArrayParams(xml_in, path));
    }

    //! Name to be used
    const std::string name = "ZOLO_NEF";

    //! Register the Wilson fermact
    const bool registered = Chroma::TheWilsonTypeFermActArrayFactory::Instance().registerObject(name, createFermAct)
      & Chroma::TheEvenOddPrecDWFermActBaseArrayFactory::Instance().registerObject(name, createDWFermAct); 
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
      read(paramtop, "ApproximationType", approximation_type);
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



  //! Check stuff
  void EvenOddPrecZoloNEFFermActArray::init()
  {
  }

  void EvenOddPrecZoloNEFFermActArray::initCoeffs(multi1d<Real>& b5_arr,
						  multi1d<Real>& c5_arr,
						  Handle<const ConnectState>& state) const
  {
    b5_arr.resize(N5);
    c5_arr.resize(N5);

    Real approxMin;
    Real approxMax;
    try {
      const OverlapConnectState& ov_state=dynamic_cast<const OverlapConnectState&>(*state);
      approxMin = ov_state.getApproxMin();
      approxMax = ov_state.getApproxMax();

      // You can do e-value stuff here later
    }
    catch( bad_cast ) {
      QDPIO::cerr << "Failed to cast ConnectState to OverlapConnectState" <<endl;
      QDPIO::cerr << " in EvenOddPrecNEFFermActArray::linOp() " << endl;
      QDP_abort(1);
    }


    Real epsilon;
    Real scale_fact;
    zolotarev_data *rdata;

    switch( approximation_type ) { 
    case COEFF_TYPE_ZOLOTAREV: 
    {
      epsilon = approxMin / approxMax;
      QDPIO::cout << "Epsilon= " << epsilon << endl << flush;
      rdata=zolotarev(toFloat(epsilon), N5, 0);
      scale_fact = approxMax;
      break;
    }
    case COEFF_TYPE_TANH_UNSCALED:
    {
      epsilon = approxMin;
      scale_fact=Real(1);
      rdata=higham(toFloat(epsilon), N5);
      break;
    }
    default:
    {
      QDPIO::cout << "Unsupported Coeff Type : " << approximation_type << endl;
      QDP_abort(1);
      break;
    }
    };

    if( rdata->n != N5 ) { 
      QDPIO::cerr << "Error:rdata->n != N5" << endl;
      QDP_abort(1);
    }

    multi1d<Real> gamma(N5);
    for(int i=0; i < N5; i++) { 
      gamma[i] = Real(rdata->gamma[i])*scale_fact;
    }

    Real maxerr=rdata->Delta;
    zolotarev_free(rdata);

    QDPIO::cout << "Initing Zolo NEF Linop: N5=" << N5 
		<< " b5+c5=" << b5+c5 
		<<  " b5-c5=" << b5-c5 
		<<  " epsilon=" << epsilon 
                <<  " scale=" << approxMax 
		<<  " maxerr=" << maxerr << endl;

    
    for(int i = 0; i < N5; i++) { 
      Real omega = Real(1)/gamma[i];

      b5_arr[i] = Real(0.5)*( (omega + Real(1))*b5 + (omega - Real(1))*c5 );
      c5_arr[i] = Real(0.5)*( (omega - Real(1))*b5 + (omega + Real(1))*c5 );
    
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
  const EvenOddPrecDWLinOpBaseArray<LatticeFermion>* 
  EvenOddPrecZoloNEFFermActArray::precLinOp(Handle<const ConnectState> state,
					    const Real& m_q) const
  {
    multi1d<Real> b5_arr;
    multi1d<Real> c5_arr;

    // Cast the state up to an overlap state
    initCoeffs(b5_arr,c5_arr,state);
    
    return new EvenOddPrecGenNEFDWLinOpArray(state->getLinks(),OverMass,b5_arr,c5_arr,m_q,N5);
  }

  //! Produce an unpreconditioned linear operator for this action with arbitrary quark mass
  const UnprecDWLinOpBaseArray<LatticeFermion>* 
  EvenOddPrecZoloNEFFermActArray::unprecLinOp(Handle<const ConnectState> state,
					      const Real& m_q) const
  {
    multi1d<Real> b5_arr;
    multi1d<Real> c5_arr;

    // Cast the state up to an overlap state
    initCoeffs(b5_arr,c5_arr,state);
    
    return new UnprecNEFDWLinOpArray(state->getLinks(),OverMass,b5_arr,c5_arr,m_q,N5);
  }

  //! Create a ConnectState with just the gauge fields
  const OverlapConnectState*
  EvenOddPrecZoloNEFFermActArray::createState(const multi1d<LatticeColorMatrix>& u_) const
  {
    const OverlapConnectState *ret_val;
    try { 
      ret_val = 
	OverlapConnectStateEnv::createOverlapState(u_, 
						   getFermBC()
	  );
    } 
    catch(const string& e) { 
      QDPIO::cerr << "Caught Exception: " << e << endl;
      QDP_abort(1);
    }
    
    return ret_val;
  }
  
  //! Create a ConnectState with just the gauge fields, and a lower
  //  approximation bound
  const OverlapConnectState*
  EvenOddPrecZoloNEFFermActArray::createState(const multi1d<LatticeColorMatrix>& u_,
					      const Real& approxMin_) const 
  {
    const OverlapConnectState *ret_val;
    try { 
      ret_val = 
	OverlapConnectStateEnv::createOverlapState(u_, 
						   getFermBC(),
						   approxMin_
	  );
    } 
    catch(const string& e) { 
      QDPIO::cerr << "Caught Exception: " << e << endl;
      QDP_abort(1);
    }

    return ret_val;
  }

  //! Create a connect State with just approximation range bounds
  const OverlapConnectState*
  EvenOddPrecZoloNEFFermActArray::createState(const multi1d<LatticeColorMatrix>& u_,
					      const Real& approxMin_,
					      const Real& approxMax_) const
  {
    const OverlapConnectState *ret_val;
    try { 
      ret_val = 
	OverlapConnectStateEnv::createOverlapState(u_, 
						   getFermBC(),
						   approxMin_,
						   approxMax_
	  );
    } 
    catch(const string& e) { 
      QDPIO::cerr << "Caught Exception: " << e << endl;
      QDP_abort(1);
    }

    return ret_val;
  }
  
  //! Create OverlapConnectState with eigenvalues/vectors
  const OverlapConnectState*
  EvenOddPrecZoloNEFFermActArray::createState(const multi1d<LatticeColorMatrix>& u_,
					      const multi1d<Real>& lambda_lo_, 
					      const multi1d<LatticeFermion>& evecs_lo_,
					      const Real& lambda_hi_) const
  {
    const OverlapConnectState *ret_val;
    try { 
      ret_val = 
	OverlapConnectStateEnv::createOverlapState(u_, 
						   getFermBC(),
						   lambda_lo_, 
						   evecs_lo_, 
						   lambda_hi_);
      
    } 
    catch(const string& e) { 
      QDPIO::cerr << "Caught Exception: " << e << endl;
      QDP_abort(1);
    }
    
    return ret_val;
  }    
  
  //! Create OverlapConnectState from XML
  const OverlapConnectState*
  EvenOddPrecZoloNEFFermActArray::createState(const multi1d<LatticeColorMatrix>& u_,
					      XMLReader& state_info_xml,
					      const string& state_info_path) const
  {
    multi1d<LatticeColorMatrix> u_tmp = u_;
    
    // HACK UP A LINEAR OPERATOR TO CHECK EIGENVALUES/VECTORS WITH
    getFermBC().modifyU(u_tmp);
    Handle< const ConnectState > state_aux = new SimpleConnectState(u_tmp);
    Handle< FermBC<LatticeFermion> > fbc_aux = new PeriodicFermBC<LatticeFermion>;
    Real aux_over_mass = -OverMass;
    UnprecWilsonFermAct S_aux( fbc_aux, aux_over_mass );
    
    Handle< const LinearOperator<LatticeFermion> > Maux = 
      S_aux.gamma5HermLinOp(state_aux);
    
    const OverlapConnectState *ret_val;
    
    try {
      ret_val = 
	OverlapConnectStateEnv::createOverlapState(
	  u_,
	  getFermBC(),
	  state_info_xml,
	  state_info_path,
	  *Maux);
    }
    catch(const string& e) { 
      QDPIO::cerr << "Caught Exception: " << e << endl;
      QDP_abort(1);
    }
    
    return ret_val;
  }

  //! Create OverlapConnectState from XML
  const OverlapConnectState*
  EvenOddPrecZoloNEFFermActArray::createState(const multi1d<LatticeColorMatrix>& u_,
					      const OverlapStateInfo& state_info) const
  {

    multi1d<LatticeColorMatrix> u_tmp = u_;
    getFermBC().modifyU(u_tmp);
    Handle< const ConnectState > state_aux = new SimpleConnectState(u_tmp);
    Handle< FermBC<LatticeFermion> > fbc_aux = new PeriodicFermBC<LatticeFermion>;

    Real aux_over_mass = -OverMass;
    UnprecWilsonFermAct S_aux( fbc_aux, aux_over_mass );
    
    Handle< const LinearOperator<LatticeFermion> > Maux = 
      S_aux.gamma5HermLinOp(state_aux);
    
    const OverlapConnectState* ret_val;    
    try {
      ret_val = 
	OverlapConnectStateEnv::createOverlapState(
	  u_,
	  getFermBC(),
	  state_info,
	  *Maux);
    }
    catch(const string& e) { 
      QDPIO::cerr << "Caught Exception: " << e << endl;
      QDP_abort(1);
    }
    
    return ret_val;

  }


}
