// $Id: unprec_zolo_nef_fermact_array_w.cc,v 1.15 2005-05-28 22:37:42 edwards Exp $
/*! \file
 *  \brief Unpreconditioned NEF fermion action
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/unprec_zolo_nef_fermact_array_w.h"
#include "actions/ferm/fermacts/unprec_wilson_fermact_w.h"
#include "actions/ferm/linop/unprec_nef_linop_array_w.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermbcs/fermbcs_w.h"

#include "actions/ferm/fermacts/zolotarev.h"

#include "actions/ferm/qprop/quarkprop4_w.h"
#include "actions/ferm/qprop/nef_quarkprop4_w.h"

namespace Chroma
{
  //! Hooks to register the class with the fermact factory
  namespace UnprecZoloNEFFermActArrayEnv
  {
    //! Callback function
    WilsonTypeFermAct5D<LatticeFermion,multi1d<LatticeColorMatrix> >* createFermAct5D(XMLReader& xml_in,
										      const std::string& path)
    {
      return new UnprecZoloNEFFermActArray(WilsonTypeFermBCArrayEnv::reader(xml_in, path), 
					   UnprecZoloNEFFermActArrayParams(xml_in, path));
    }

    //! Callback function
    /*! Differs in return type */
    FermionAction<LatticeFermion>* createFermAct(XMLReader& xml_in,
						 const std::string& path)
    {
      return createFermAct5D(xml_in, path);
    }

    //! Name to be used
    const std::string name = "UNPRECONDITIONED_ZOLO_NEF";

    //! Register all the factories
    bool registerAll()
    {
      return Chroma::TheFermionActionFactory::Instance().registerObject(name, createFermAct)
	   & Chroma::TheWilsonTypeFermAct5DFactory::Instance().registerObject(name, createFermAct5D);
    }

    //! Register the fermact
    const bool registered = registerAll();
  }


  //! Read parameters
  UnprecZoloNEFFermActArrayParams::UnprecZoloNEFFermActArrayParams(XMLReader& xml, 
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
  void read(XMLReader& xml, const string& path, UnprecZoloNEFFermActArrayParams& param)
  {
    UnprecZoloNEFFermActArrayParams tmp(xml, path);
    param = tmp;
  }



  void UnprecZoloNEFFermActArray::initCoeffs(multi1d<Real>& b5_arr,
					     multi1d<Real>& c5_arr,
					     Handle<const ConnectState> state) const
  {
    b5_arr.resize(params.N5);
    c5_arr.resize(params.N5);

    Real approxMin;
    Real approxMax;
    try {
      const OverlapConnectState& ov_state=dynamic_cast<const OverlapConnectState&>(*state);
      approxMin = (ov_state.getEigVal().size() != 0) ? ov_state.getApproxMin() : params.ApproxMin;
      approxMax = (ov_state.getEigVal().size() != 0) ? ov_state.getApproxMax() : params.ApproxMax;

      // You can do e-value stuff here later
    }
    catch( bad_cast ) {
      QDPIO::cerr << "Failed to cast ConnectState to OverlapConnectState" <<endl;
      QDPIO::cerr << " in UnprecNEFFermActArray::linOp() " << endl;
      QDP_abort(1);
    }

    zolotarev_data *rdata;
    Real epsilon;
    Real scale_fac;

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

    Real maxerr = rdata->Delta;

    multi1d<Real> gamma(params.N5);
    for(int i=0; i < params.N5; i++) { 
      gamma[i] = Real(rdata->gamma[i])*scale_fac;
    }

    zolotarev_free(rdata);

    QDPIO::cout << "Initing Zolo NEF Linop: N5=" << params.N5 
		<< " b5+c5=" << params.b5+params.c5 
		<<  " b5-c5=" << params.b5-params.c5 
		<<  " epsilon=" << epsilon
                <<  " scale=" << approxMax 
		<<  " Maxerr=" << maxerr << endl;

    for(int i=0; i < params.N5; i++) { 
      QDPIO::cout << "gamma[" << i << "] = " << gamma[i] << endl;
    }
    
    for(int i = 0; i < params.N5; i++) { 
      Real omega = Real(1)/gamma[i];

      b5_arr[i] = Real(0.5)*( (omega + Real(1))*params.b5 + (omega - Real(1))*params.c5 );
      c5_arr[i] = Real(0.5)*( (omega - Real(1))*params.b5 + (omega + Real(1))*params.c5 );

      QDPIO::cout << "i=" << i << " omega["<<i<<"]=" << omega
		  << " b5["<< i << "] ="<< b5_arr[i] 
		  << " c5["<< i << "] ="<< c5_arr[i] << endl;
      QDPIO::cout << "i=" << i << " b5["<<i<<"]+c5["<<i<<"]/omega["<<i<<"] =" 
		  << (b5_arr[i]+c5_arr[i])/omega
		  << " b5["<<i<<"]-c5["<<i<<"]=" << b5_arr[i]-c5_arr[i] << endl;
    }
  }


  //! Produce an unpreconditioned linear operator for this action with arbitrary quark mass
  const UnprecDWLikeLinOpBaseArray<LatticeFermion,multi1d<LatticeColorMatrix> >* 
  UnprecZoloNEFFermActArray::unprecLinOp(Handle<const ConnectState> state, 
					 const Real& m_q) const
  {
    multi1d<Real> b5_arr;
    multi1d<Real> c5_arr;
    
    // Cast the state up to an overlap state
    initCoeffs(b5_arr,c5_arr,state);
    
    return new UnprecNEFDWLinOpArray(state->getLinks(),params.OverMass,b5_arr,c5_arr,m_q,params.N5);
  }


  //! Create a ConnectState with just the gauge fields
  const OverlapConnectState*
  UnprecZoloNEFFermActArray::createState(const multi1d<LatticeColorMatrix>& u_) const
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
  UnprecZoloNEFFermActArray::createState(const multi1d<LatticeColorMatrix>& u_,
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
  UnprecZoloNEFFermActArray::createState(const multi1d<LatticeColorMatrix>& u_,
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
  UnprecZoloNEFFermActArray::createState(const multi1d<LatticeColorMatrix>& u_,
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
  UnprecZoloNEFFermActArray::createState(const multi1d<LatticeColorMatrix>& u_,
				      XMLReader& state_info_xml,
				      const string& state_info_path) const
  {
    multi1d<LatticeColorMatrix> u_tmp = u_;
    
    // HACK UP A LINEAR OPERATOR TO CHECK EIGENVALUES/VECTORS WITH
    getFermBC().modifyU(u_tmp);
    Handle< const ConnectState > state_aux = new SimpleConnectState(u_tmp);
    Handle< FermBC<LatticeFermion> > fbc_aux = new PeriodicFermBC<LatticeFermion>;
    Real aux_over_mass = -params.OverMass;
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
  UnprecZoloNEFFermActArray::createState(const multi1d<LatticeColorMatrix>& u_,
				      const OverlapStateInfo& state_info) const
  {

    multi1d<LatticeColorMatrix> u_tmp = u_;
    getFermBC().modifyU(u_tmp);
    Handle< const ConnectState > state_aux = new SimpleConnectState(u_tmp);
    Handle< FermBC<LatticeFermion> > fbc_aux = new PeriodicFermBC<LatticeFermion>;

    Real aux_over_mass = -params.OverMass;
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

  
  // Given a complete propagator as a source, this does all the inversions needed
  void 
  UnprecZoloNEFFermActArray::quarkProp(LatticePropagator& q_sol, 
				       XMLWriter& xml_out,
				       const LatticePropagator& q_src,
				       int t_src, int j_decay,
				       Handle<const ConnectState> state,
				       const InvertParam_t& invParam,
				       bool nonRelProp,
				       bool obsvP,
				       int& ncg_had)
  {
    if (obsvP)
      nef_quarkProp4(q_sol, xml_out, q_src, t_src, j_decay, *this, state, invParam, ncg_had);
    else
    {
      Handle< const SystemSolver<LatticeFermion> > qprop(this->qprop(state,invParam));
      quarkProp4(q_sol, xml_out, q_src, qprop, nonRelProp, ncg_had);
    }
  }


}
