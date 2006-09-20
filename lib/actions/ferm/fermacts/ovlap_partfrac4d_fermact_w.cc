// $Id: ovlap_partfrac4d_fermact_w.cc,v 3.2 2006-09-20 20:27:59 edwards Exp $
/*! \file
 *  \brief 4D Zolotarev variant of Overlap-Dirac operator
 */
#include <iostream>
#include <sstream>
#include <chromabase.h>
#include <linearop.h>
#include <string>
#include "io/enum_io/enum_io.h"  // Read/Write OverlapInnerSolver Type
// #include "io/overlap_state_info.h"

// For creating auxiliary actions
#include "actions/ferm/fermacts/unprec_wilson_fermact_w.h"

// My defs
#include "actions/ferm/fermacts/ovlap_partfrac4d_fermact_w.h"

// the zolo coeffs
#include "actions/ferm/fermacts/zolotarev.h"

// Linops I can make
#include "actions/ferm/linop/lovlapms_w.h"
#include "actions/ferm/linop/lovlap_double_pass_w.h"
#include "actions/ferm/linop/lovddag_w.h"
#include "actions/ferm/linop/lovddag_double_pass_w.h"
#include "actions/ferm/linop/lg5eps_w.h"
#include "actions/ferm/linop/lg5eps_double_pass_w.h"
#include "meas/eig/ischiral_w.h"


#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermstates/simple_fermstate.h"
#include "actions/ferm/fermstates/ferm_createstate_reader_w.h"
#include "actions/ferm/fermbcs/fermbcs_reader_w.h"

namespace Chroma
{

  //! Hooks to register the class with the fermact factory
  namespace OvlapPartFrac4DFermActEnv
  {
    //! Callback function
    WilsonTypeFermAct<LatticeFermion, 
		      multi1d<LatticeColorMatrix>,
		      multi1d<LatticeColorMatrix> >* createFermAct4D(XMLReader& xml_in,
								     const std::string& path)
    {
      return new OvlapPartFrac4DFermAct(WilsonTypeFermBCEnv::reader(xml_in, path), 
					OvlapPartFrac4DFermActParams(xml_in, path));
    }

    //! Callback function
    /*! Differs in return type */
    FermionAction<LatticeFermion,
		  multi1d<LatticeColorMatrix>,
		  multi1d<LatticeColorMatrix> >* createFermAct(XMLReader& xml_in,
							       const std::string& path)
    {
      return createFermAct4D(xml_in, path);
    }

    //! Name to be used
    const std::string name ="OVERLAP_PARTIAL_FRACTION_4D";

    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::TheFermionActionFactory::Instance().registerObject(name, createFermAct);
	success &= Chroma::TheWilsonTypeFermActFactory::Instance().registerObject(name, createFermAct4D);
	registered = true;
      }
      return success;
    }
  }



  OvlapPartFrac4DFermActParams::OvlapPartFrac4DFermActParams(XMLReader& xml, const std::string& path) : ReorthFreqInner(10), inner_solver_type(OVERLAP_INNER_CG_SINGLE_PASS)
  {
    XMLReader in(xml, path);

    try 
    { 
      if(in.count("AuxFermAct") == 1 )
      { 
	XMLReader xml_tmp(in, "AuxFermAct");
	std::ostringstream os;
	xml_tmp.print(os);
	AuxFermAct = os.str();
      }
      else
      {
	throw std::string("No auxilliary action");
      }
  	
      read(in, "Mass", Mass);
      read(in, "RatPolyDeg", RatPolyDeg);

      if(in.count("RatPolyDegPrecond") == 1 ) { 
	read(in, "RatPolyDegPrecond", RatPolyDegPrecond);
      }
      else {
	RatPolyDegPrecond = RatPolyDeg;
      }

      if( in.count("ApproximationType") == 1 ) { 
	read(in, "ApproximationType", approximation_type);
      }
      else { 
	// Default coeffs are Zolotarev
	approximation_type = COEFF_TYPE_ZOLOTAREV;
      }

      read(in, "ApproxMin", approxMin);
      read(in, "ApproxMax", approxMax);

      read(in, "InnerSolve/MaxCG", invParamInner.MaxCG);
      read(in, "InnerSolve/RsdCG", invParamInner.RsdCG);
      if( in.count("InnerSolve/ReorthFreq") == 1 ) {
	read(in, "InnerSolve/ReorthFreq", ReorthFreqInner);
      }
      else {
	// ReorthFreqInner = 10; // Some default
	// This is now set automagically -- constructor initialisation
      }

      if( in.count("InnerSolve/SolverType") == 1 ) { 
	read(in, "InnerSolve/SolverType", inner_solver_type);
      }
      else {
	// inner_solver_type = OVERLAP_INNER_CG_SINGLE_PASS;
	// This is now set automagically -- constructor initialisation
      }

      read(in, "IsChiral", isChiralP);
    }
    catch( const string &e ) {
      QDPIO::cerr << "Caught Exception reading Zolo4D Fermact params: " << e << endl;
      QDP_abort(1);
    }
  }

  void write(XMLWriter& xml_out, const string& path, const OvlapPartFrac4DFermActParams& p)
  {
    if ( path != "." ) { 
      push( xml_out, path);
    }
  
    xml_out << p.AuxFermAct;
    write(xml_out, "Mass", p.Mass);
    write(xml_out, "RatPolyDeg", p.RatPolyDeg);
    write(xml_out, "RatPolyDegPrecond", p.RatPolyDegPrecond);
    push(xml_out, "InnerSolve");
    write(xml_out, "MaxCG", p.invParamInner.MaxCG);
    write(xml_out, "RsdCG", p.invParamInner.RsdCG);
    write(xml_out, "ReorthFreq", p.ReorthFreqInner);
    write(xml_out, "SolverType", p.inner_solver_type);
    write(xml_out, "ApproximationType", p.approximation_type);
    write(xml_out, "ApproxMin", p.approxMin);
    write(xml_out, "ApproxMax", p.approxMax);
    write(xml_out, "IsChiral", p.isChiralP); 
    pop(xml_out);

//    write(xml_out, "StateInfo", p.state_info);

    if( path != "." ) { 
      pop(xml_out);
    }
  }

  //! Read parameters
  void read(XMLReader& xml, const string& path, OvlapPartFrac4DFermActParams& param)
  {
    OvlapPartFrac4DFermActParams tmp(xml, path);
    param = tmp;
  }



  //! Constructor
  OvlapPartFrac4DFermAct::OvlapPartFrac4DFermAct(Handle< FermBC<T,P,Q> > fbc_, 
						 const OvlapPartFrac4DFermActParams& params_) : 
    fbc(fbc_), params(params_)
  {
    QDPIO::cout << "Constructing OvlapPartFrac4D FermAct from params" << endl;
    std::istringstream  xml_s(params.AuxFermAct);
    XMLReader  fermacttop(xml_s);
    const string fermact_path = "/AuxFermAct";

    // Sanity check
    if (fbc.operator->() == 0)
    {
      QDPIO::cerr << OvlapPartFrac4DFermActEnv::name << ": error: fbc is null" << endl;
      QDP_abort(1);
    }

    // Fake a creator. This should be cleaned up
    Handle< CreateFermState<T,P,Q> > cfs_(new CreateSimpleFermState<T,P,Q>(fbc));
    cfs = cfs_;

    struct UnprecCastFailure {
      UnprecCastFailure(std::string e) : auxfermact(e) {};
      const string auxfermact;
    };

    try
    {
      string auxfermact;
      read(fermacttop, fermact_path + "/FermAct", auxfermact);
      QDPIO::cout << "AuxFermAct: " << auxfermact << endl;

      read(fermacttop, fermact_path + "/Mass", params.AuxMass);
      QDPIO::cout << "AuxFermAct Mass: " << params.AuxMass << endl;
      // Generic Wilson-Type stuff
      FermionAction<T,P,Q>* S_f =
	TheFermionActionFactory::Instance().createObject(auxfermact,
							 fermacttop,
							 fermact_path);

      UnprecWilsonTypeFermAct<T,P,Q>* S_aux = 
	dynamic_cast<UnprecWilsonTypeFermAct<T,P,Q>*>(S_f);

      // Dumbass User specifies something that is not UnpreWilsonTypeFermAct
      // dynamic_cast MUST be checked for 0
      if( S_aux == 0 ) throw UnprecCastFailure(auxfermact);
     

      // Drop AuxFermAct into a Handle immediately.
      // This should free things up at the end
      Handle<UnprecWilsonTypeFermAct<T,P,Q> >  S_w(S_aux);
      Mact = S_w;
    }
    catch( const UnprecCastFailure& e) {

      // Breakage Scenario
      QDPIO::cerr << "Unable to upcast auxiliary fermion action to "
		  << "UnprecWilsonTypeFermAct " << endl;
      QDPIO::cerr << OvlapPartFrac4DFermActEnv::name << " does not support even-odd preconditioned "
		  << "auxiliary FermActs" << endl;
      QDPIO::cerr << "You passed : " << endl;
      QDPIO::cerr << e.auxfermact << endl;
      QDP_abort(1);
    }
    catch (const std::exception& e) {
      // General breakage Scenario
      QDPIO::cerr << "Error reading data: " << e.what() << endl;
      throw;
    }

  }


  //! Creation routine
  /*! */
  void 
  OvlapPartFrac4DFermAct::init(int& numroot, 
			       Real& coeffP, 
			       multi1d<Real>& resP,
			       multi1d<Real>& rootQ, 
			       int& NEig, 
			       multi1d<Real>& EigValFunc,
			       const EigenConnectState& state ) const
  {
    /* A scale factor which should bring the spectrum of the hermitian
       Wilson Dirac operator H into |H| < 1. */
    Real scale_fac;
  
    /* Contains all the data necessary for Zolotarev partial fraction */
    /* -------------------------------------------------------------- */
    zolotarev_data *rdata ;
    /* The lower (positive) interval bound for the approximation 
       interval [-1,-eps] U [eps,1] */

    Real eps;
    /* The type of the approximation R(x): 
       type = 0 -> R(x) = 0        at x = 0 
       type = 1 -> R(x) = infinity at x = 0 */

    int type;
    /* The maximal error of the approximation in the interval 
       [-1,-eps] U [eps,1]*/

    Real maxerr;

    // The residual for the solutions of the multi-shift linear system
    // I put this in the class constructor 

    //  RsdCGinner = 1.0e-7;  // Hardwired the accuracy


    /* Hermitian 4D overlap operator 1/2 ( 1 + Mass + (1 - Mass) gamma5 * sgn(H)) 
       using a partial fraction expansion of the optimal rational function
       approximation to sgn. Here, H = 1/2 * gamma_5 * (1/kappa - D'). 
       The coefficients are computed by Zolotarev's formula. */

    int NEigVal = state.getNEig();

    /* The operator gamma_5 * M with the M constructed here has its eigenvalues
       in the range m/(m + Nd) <= |gamma_5 * M| <= (m + 2*Nd)/(m + Nd) (in the 
       free case) where here m is arbitrary.
       So if we multiply M by a factor scale_fac = (m + Nd)/(m + 2*Nd) we have
       |gamma_5 * M| <= 1. */

    // REmove this evil dublication later please.
    NEig = NEigVal;

    switch(params.approximation_type) { 
    case COEFF_TYPE_ZOLOTAREV:
      scale_fac = Real(1) / params.approxMax;
      eps = params.approxMin * scale_fac;

      QDPIO::cout << "Initing Linop with Zolotarev Coefficients" << endl;
      QDPIO::cout << "  MaxCGInner =  " << params.invParamInner.MaxCG << endl;
      QDPIO::cout << "  RsdCGInner =  " << params.invParamInner.RsdCG << endl;
      QDPIO::cout << "  NEigVal    =  " << NEigVal << endl;
      
      /* Below, when we fill in the coefficents for the partial fraction, 
	 we include this factor, say t, appropriately, i.e.
	 R(x) = alpha[da] * t * x + sum(alpha[j] * t * x / (t^2 * x^2 - ap[j]), 
	 j = 0 .. da-1)
	 = (alpha[da] + + sum(alpha[j] / (x^2 - ap[j] / t^2) ) / t^2 ) * t * x 
      */
      
      /* ZOLOTAREV_4D uses Zolotarev's formula for the coefficients. 
	 The coefficents produced are for an optimal uniform approximation
	 to the sign-function in the interval [-1,-eps] U [eps,1] and of order n. 
	 type can be set to 0 or 1 corresponding to an approximation which is 
	 is zero or infinite at x = 0, respectively. 
	 Here we are interested in the partial fraction form 
	 
	 R(x) = alpha[da] * x + sum(alpha[j] * x / (x^2 - ap[j]), j = 0 .. da-1) 
	 
	 where da = dd for type 0 and da = dd + 1 with ap[dd] = 0 for type 1. 
      */
      type = 0;
      rdata = zolotarev(toFloat(eps), params.RatPolyDeg, type);
      if( rdata == 0x0 ) { 
	QDPIO::cerr << "Failed to get Zolo Coeffs" << endl;
        QDP_abort(1);
      } 
      break;

    case COEFF_TYPE_TANH:
      scale_fac = Real(1) / params.approxMax;
      eps = params.approxMin * scale_fac;

      QDPIO::cout << "Initing Linop with Higham Rep tanh Coefficients" << endl;
      QDPIO::cout << "  MaxCGInner =  " << params.invParamInner.MaxCG << endl;
      QDPIO::cout << "  RsdCGInner =  " << params.invParamInner.RsdCG << endl;
      QDPIO::cout << "  NEigVal    =  " << NEigVal << endl;
      
      /* Below, when we fill in the coefficents for the partial fraction, 
	 we include this factor, say t, appropriately, i.e.
	 R(x) = alpha[da] * t * x + sum(alpha[j] * t * x / (t^2 * x^2 - ap[j]), 
	 j = 0 .. da-1)
	 = (alpha[da] + + sum(alpha[j] / (x^2 - ap[j] / t^2) ) / t^2 ) * t * x 
      */
      
      /*  use the tanh formula (Higham Rep) for the coefficients. 
	 The coefficents produced are for the tanh approximation
	 to the sign-function in the interval [-1,-eps] U [eps,1] and of order n.	 R(x) = alpha[da] * x + sum(alpha[j] * x / (x^2 - ap[j]), j = 0 .. da-1) 
	 where da = dd for type 0 and da = dd + 1 with ap[dd] = 0 for type 1. 
      */
      rdata = higham(toFloat(eps), params.RatPolyDeg);
      break;

    case COEFF_TYPE_TANH_UNSCALED:
      scale_fac = Real(1) ;
      eps = params.approxMin;

      QDPIO::cout << "Initing Linop with Unscaled Higham Rep tanh Coefficients" << endl;
      QDPIO::cout << "  MaxCGInner =  " << params.invParamInner.MaxCG << endl;
      QDPIO::cout << "  RsdCGInner =  " << params.invParamInner.RsdCG << endl;
      QDPIO::cout << "  NEigVal    =  " << NEigVal << endl;
      
      /* Below, when we fill in the coefficents for the partial fraction, 
	 we include this factor, say t, appropriately, i.e.
	 R(x) = alpha[da] * t * x + sum(alpha[j] * t * x / (t^2 * x^2 - ap[j]), 
	 j = 0 .. da-1)
	 = (alpha[da] + + sum(alpha[j] / (x^2 - ap[j] / t^2) ) / t^2 ) * t * x 
      */
      
      /*  use the tanh formula (Higham Rep) for the coefficients. 
	 The coefficents produced are for the tanh approximation
	 to the sign-function in the interval [-1,-eps] U [eps,1] and of order n.	 R(x) = alpha[da] * x + sum(alpha[j] * x / (x^2 - ap[j]), j = 0 .. da-1) 
	 where da = dd for type 0 and da = dd + 1 with ap[dd] = 0 for type 1. 
      */
      rdata = higham(toFloat(eps), params.RatPolyDeg);
      break;

    default:
      // The map system should ensure that we never get here but 
      // just for style
      QDPIO::cerr << "Unknown coefficient type: " << params.approximation_type
		  << endl;
      QDP_abort(1);
    }
    
    maxerr = (Real)(rdata -> Delta);
    QDPIO::cout << "Maxerr " << maxerr << flush << endl; 
      /*
	push(my_writer, "ZolotarevApprox");
	write(my_writer, "eps", eps);
	write(my_writer, "scale_fac", scale_fac);
	write(my_writer, "RatPolyDeg", params.RatPolyDeg);
	write(my_writer, "type", type);
	write(my_writer, "maxerr", maxerr);
	write(my_writer, "InnerSolverType", params.inner_solver_type);
	pop(my_writer);
      */

      /* The number of residuals and poles */
      /* Allocate the roots and residua */
    numroot = rdata -> dd;
    /* The roots, i.e., the shifts in the partial fraction expansion */
    rootQ.resize(numroot);
    /* The residuals in the partial fraction expansion */
    resP.resize(numroot);
 
    /* Fill in alpha[0] = alpha[da] if it is not zero*/
    coeffP = Real(0);
    coeffP = rdata -> alpha[rdata -> da - 1];
    /* The coefficients from the partial fraction.
       Here, we write them out for the sake of bookkeeping. */
    for(int n=0; n < numroot; ++n) {
      resP[n] = rdata -> alpha[n];
      rootQ[n] = rdata -> ap[n];
      rootQ[n] = -rootQ[n];
    }
    /*
    push(my_writer,"ZolotarevPartFrac");
    write(my_writer, "scale_fac", scale_fac);
    write(my_writer, "coeffP", coeffP);
    write(my_writer, "resP", resP);
    write(my_writer, "rootQ", rootQ);
    pop(my_writer);
    */

    /* Now fill in the coefficients for real, i.e., taking the rescaling
       into account */
    /* Fill in alpha[0] = alpha[da] if it is not zero*/
    coeffP = rdata -> alpha[rdata -> da - 1] * scale_fac;
    /* Fill in the coefficients for the roots and the residua */
    /* Make sure that the smallest shift is in the last value rootQ(numroot-1)*/
    Real t = Real(1) / (scale_fac * scale_fac);
    for(int n=0; n < numroot; ++n) {
    
      resP[n] = rdata -> alpha[n] / scale_fac;
      rootQ[n] = rdata -> ap[n];
      rootQ[n] = -(rootQ[n] * t);
    }
  
    /* Write them out into the namelist */
    /*
    push(my_writer,"ZolotarevPartFracResc");
    write(my_writer, "scale_fac", scale_fac);
    write(my_writer, "coeffP", coeffP);
    write(my_writer, "resP", resP);
    write(my_writer, "rootQ", rootQ);
    pop(my_writer);
    */
  
    QDPIO::cout << "PartFracApprox 4d n=" << params.RatPolyDeg << " scale=" << scale_fac
		<< " coeff=" << coeffP << " Nwils= " << NEigVal <<" Mass="
		<< params.Mass << " Rsd=" << params.invParamInner.RsdCG << endl;
  
    QDPIO::cout << "Approximation on [-1,-eps] U [eps,1] with eps = " << eps <<
      endl;
    QDPIO::cout << "Maximum error |R(x) - sqn(x)| <= " << maxerr << endl;
  
    switch( params.approximation_type) {
    case COEFF_TYPE_ZOLOTAREV:
      QDPIO::cout << "Coefficients from Zolotarev" << endl;

      if(type == 0) {
	QDPIO::cout << "Approximation type " << type << " with R(0) = 0"
		    << endl;
      }
      else {
	QDPIO::cout << "Approximation type " << type << " with R(0) =  infinity"                    << endl;
      }

      break;
    case COEFF_TYPE_TANH:
      QDPIO::cout << "Coefficients from Higham Tanh representation" << endl;
      break;

    case COEFF_TYPE_TANH_UNSCALED:
      QDPIO::cout << "Coefficients from Unscaled Higham Tanh representation" << endl;
      break;

    default:
      QDPIO::cerr << "Unknown coefficient type " << params.approximation_type 
		  << endl;
      break;
    }


    switch(params.inner_solver_type) { 
    case OVERLAP_INNER_CG_SINGLE_PASS:
      QDPIO::cout << "Using Single Pass Inner Solver" << endl;
      break;
    case OVERLAP_INNER_CG_DOUBLE_PASS:
      QDPIO::cout << "Using Neuberger/Chu Double Pass Inner Solver" << endl;
      break;
    default:
      QDPIO::cerr << "Unknown inner solver type " << endl;
      QDP_abort(1);
    }

    QDPIO::cout << "Number of poles= " << numroot << endl;
    QDPIO::cout << "Overall Factor=  " << coeffP << endl;
    QDPIO::cout << "Numerator coefficients:" << endl;
    for(int n=0; n < numroot; n++) { 
      QDPIO::cout <<"  resP[" << n << "]= " << resP[n] << endl;
    }
    QDPIO::cout << "Denominator roots: " << endl;
    for(int n=0; n < numroot; n++) { 
      QDPIO::cout <<"  rootQ[" << n<< "]= " << rootQ[n] << endl;
    }
 
    /* We will also compute the 'function' of the eigenvalues */
    /* for the Wilson vectors to be projected out. */
    if (NEig > 0)
    {
      for(int i = 0; i < NEigVal; i++)
      {
	if (toBool(state.getEvalues()[i] > 0.0))
	  EigValFunc[i] = 1.0;
	else if (toBool(state.getEvalues()[i] < 0.0))
	  EigValFunc[i] = -1.0;
	else
	  EigValFunc[i] = 0.0;
      }
    }


    // Free the arrays allocate by Tony's zolo
    zolotarev_free(rdata);
  }

  void 
  OvlapPartFrac4DFermAct::initPrec(int& numroot, 
			       Real& coeffP, 
			       multi1d<Real>& resP,
			       multi1d<Real>& rootQ, 
			       int& NEig, 
			       multi1d<Real>& EigValFunc,
			       const EigenConnectState& state ) const
  {
    /* A scale factor which should bring the spectrum of the hermitian
       Wilson Dirac operator H into |H| < 1. */
    Real scale_fac;
    XMLBufferWriter my_writer;
  
    /* Contains all the data necessary for Zolotarev partial fraction */
    /* -------------------------------------------------------------- */
    zolotarev_data *rdata ;
    /* The lower (positive) interval bound for the approximation 
       interval [-1,-eps] U [eps,1] */

    Real eps;
    /* The type of the approximation R(x): 
       type = 0 -> R(x) = 0        at x = 0 
       type = 1 -> R(x) = infinity at x = 0 */

    int type;
    /* The maximal error of the approximation in the interval 
       [-1,-eps] U [eps,1]*/

    Real maxerr;

    // The residual for the solutions of the multi-shift linear system
    // I put this in the class constructor 

    //  RsdCGinner = 1.0e-7;  // Hardwired the accuracy


    /* Hermitian 4D overlap operator 1/2 ( 1 + Mass + (1 - Mass) gamma5 * sgn(H)) 
       using a partial fraction expansion of the optimal rational function
       approximation to sgn. Here, H = 1/2 * gamma_5 * (1/kappa - D'). 
       The coefficients are computed by Zolotarev's formula. */

    int NEigVal = state.getNEig();

    /* The operator gamma_5 * M with the M constructed here has its eigenvalues
       in the range m/(m + Nd) <= |gamma_5 * M| <= (m + 2*Nd)/(m + Nd) (in the 
       free case) where here m is arbitrary.
       So if we multiply M by a factor scale_fac = (m + Nd)/(m + 2*Nd) we have
       |gamma_5 * M| <= 1. */
    NEig = NEigVal;


    switch(params.approximation_type) { 
    case COEFF_TYPE_ZOLOTAREV:
      scale_fac = Real(1) / params.approxMax;
      eps = params.approxMin * scale_fac;
      
      QDPIO::cout << "Initing Linop with Zolotarev Coefficients" << endl;
      
      QDPIO::cout << "  MaxCGInner =  " << params.invParamInner.MaxCG << endl;
      QDPIO::cout << "  RsdCGInner =  " << params.invParamInner.RsdCG << endl;
      QDPIO::cout << "  NEigVal    =  " << NEigVal << endl;
      
      /* Below, when we fill in the coefficents for the partial fraction, 
	 we include this factor, say t, appropriately, i.e.
	 R(x) = alpha[da] * t * x + sum(alpha[j] * t * x / (t^2 * x^2 - ap[j]), 
	 j = 0 .. da-1)
	 = (alpha[da] + + sum(alpha[j] / (x^2 - ap[j] / t^2) ) / t^2 ) * t * x 
      */
      
      /* ZOLOTAREV_4D uses Zolotarev's formula for the coefficients. 
	 The coefficents produced are for an optimal uniform approximation
	 to the sign-function in the interval [-1,-eps] U [eps,1] and of order n. 
	 type can be set to 0 or 1 corresponding to an approximation which is 
	 is zero or infinite at x = 0, respectively. 
	 Here we are interested in the partial fraction form 
	 
	 R(x) = alpha[da] * x + sum(alpha[j] * x / (x^2 - ap[j]), j = 0 .. da-1) 
	 
	 where da = dd for type 0 and da = dd + 1 with ap[dd] = 0 for type 1. 
      */
      type = 0;
      rdata = zolotarev(toFloat(eps), params.RatPolyDeg, type);
      maxerr = (Real)(rdata -> Delta);
      break;
    case COEFF_TYPE_TANH:
      scale_fac = Real(1) / params.approxMax;
      eps = params.approxMin * scale_fac;

      QDPIO::cout << "Initing Linop with Higham Rep tanh Coefficients" << endl;
      
      QDPIO::cout << "  MaxCGInner =  " << params.invParamInner.MaxCG << endl;
      QDPIO::cout << "  RsdCGInner =  " << params.invParamInner.RsdCG << endl;
      QDPIO::cout << "  NEigVal    =  " << NEigVal << endl;
      
      /* Below, when we fill in the coefficents for the partial fraction, 
	 we include this factor, say t, appropriately, i.e.
	 R(x) = alpha[da] * t * x + sum(alpha[j] * t * x / (t^2 * x^2 - ap[j]), 
	 j = 0 .. da-1)
	 = (alpha[da] + + sum(alpha[j] / (x^2 - ap[j] / t^2) ) / t^2 ) * t * x 
      */
      
      /*  use the tanh formula (Higham Rep) for the coefficients. 
	 The coefficents produced are for the tanh approximation
	 to the sign-function in the interval [-1,-eps] U [eps,1] and of order n.	 R(x) = alpha[da] * x + sum(alpha[j] * x / (x^2 - ap[j]), j = 0 .. da-1) 
	 where da = dd for type 0 and da = dd + 1 with ap[dd] = 0 for type 1. 
      */
      rdata = higham(toFloat(eps), params.RatPolyDeg);
      maxerr = (Real)(rdata -> Delta);
      break;
    case COEFF_TYPE_TANH_UNSCALED:
      scale_fac = Real(1) ;
      eps = params.approxMin;
      
      QDPIO::cout << "Initing Preconditioning Linop with Unscaled Higham Rep tanh Coefficients" << endl;
      QDPIO::cout << "  MaxCGInner =  " << params.invParamInner.MaxCG << endl;
      QDPIO::cout << "  RsdCGInner =  " << params.invParamInner.RsdCG << endl;
      QDPIO::cout << "  NEigVal    =  " << NEigVal << endl;
      
      /* Below, when we fill in the coefficents for the partial fraction, 
	 we include this factor, say t, appropriately, i.e.
	 R(x) = alpha[da] * t * x + sum(alpha[j] * t * x / (t^2 * x^2 - ap[j]), 
	 j = 0 .. da-1)
	 = (alpha[da] + + sum(alpha[j] / (x^2 - ap[j] / t^2) ) / t^2 ) * t * x 
      */
      
      /*  use the tanh formula (Higham Rep) for the coefficients. 
	 The coefficents produced are for the tanh approximation
	 to the sign-function in the interval [-1,-eps] U [eps,1] and of order n.	 R(x) = alpha[da] * x + sum(alpha[j] * x / (x^2 - ap[j]), j = 0 .. da-1) 
	 where da = dd for type 0 and da = dd + 1 with ap[dd] = 0 for type 1. 
      */
      rdata = higham(toFloat(eps), params.RatPolyDeg);


    default:
      // The map system should ensure that we never get here but 
      // just for style
      QDPIO::cerr << "Unknown coefficient type: " << params.approximation_type
		  << endl;
      QDP_abort(1);
    }
  
    /* The number of residuals and poles */
    /* Allocate the roots and residua */
    numroot = rdata -> dd;
    /* The roots, i.e., the shifts in the partial fraction expansion */
    rootQ.resize(numroot);
  
    /* The residuals in the partial fraction expansion */
    resP.resize(numroot);
  
    /* Fill in alpha[0] = alpha[da] if it is not zero*/
    coeffP = 0;
    coeffP = rdata -> alpha[rdata -> da - 1];
    /* The coefficients from the partial fraction.
       Here, we write them out for the sake of bookkeeping. */
    for(int n=0; n < numroot; ++n) {
      resP[n] = rdata -> alpha[n];
      rootQ[n] = rdata -> ap[n];
      rootQ[n] = -rootQ[n];
    }
  
    /*
    push(my_writer,"ZolotarevPreconditionerPartFrac");
    write(my_writer, "scale_fac", scale_fac);
    write(my_writer, "coeffP", coeffP);
    write(my_writer, "resP", resP);
    write(my_writer, "rootQ", rootQ);
    pop(my_writer);
    */

    /* Now fill in the coefficients for real, i.e., taking the rescaling
       into account */
    /* Fill in alpha[0] = alpha[da] if it is not zero*/
    coeffP = rdata -> alpha[rdata -> da - 1] * scale_fac;
    /* Fill in the coefficients for the roots and the residua */
    /* Make sure that the smallest shift is in the last value rootQ(numroot-1)*/
    Real t = Real(1) / (scale_fac * scale_fac);
    for(int n=0; n < numroot; ++n) {
    
      resP[n] = rdata -> alpha[n] / scale_fac;
      rootQ[n] = rdata -> ap[n];
      rootQ[n] = -(rootQ[n] * t);
    }
  
  
    /* Write them out into the namelist */
    /*
    push(my_writer,"ZolotarevPreconditionerPartFracResc");
    write(my_writer, "scale_fac", scale_fac);
    write(my_writer, "coeffP", coeffP);
    write(my_writer, "resP", resP);
    write(my_writer, "rootQ", rootQ);
    pop(my_writer);
    */

  
    QDPIO::cout << "PartFrac Preconditioner 4d n=" << params.RatPolyDegPrecond << " scale=" << scale_fac
		<< " coeff=" << coeffP << " Nwils= " << NEigVal <<" Mass="
		<< params.Mass << " Rsd=" << params.invParamInner.RsdCG << endl;
  
    QDPIO::cout << "Approximation on [-1,-eps] U [eps,1] with eps = " << eps <<
      endl;
    QDPIO::cout << "Maximum error |R(x) - sqn(x)| <= " << maxerr << endl;

    switch( params.approximation_type) {
    case COEFF_TYPE_ZOLOTAREV:
      QDPIO::cout << "Coefficients from Zolotarev" << endl;

      if(type == 0) {
	QDPIO::cout << "Approximation type " << type << " with R(0) = 0"
		    << endl;
      }
      else {
	QDPIO::cout << "Approximation type " << type << " with R(0) =  infinity"                    << endl;
      }
      break;
    case COEFF_TYPE_TANH:

      QDPIO::cout << "Coefficients from Higham Tanh representation" << endl;
      break;

    case COEFF_TYPE_TANH_UNSCALED:

      QDPIO::cout << "Coefficients from Higham Tanh representation" << endl;
      break;
    default:
      QDPIO::cerr << "Unknown coefficient type " << params.approximation_type 
		  << endl;
      break;
    }
    QDPIO::cout << flush;  

    switch(params.inner_solver_type) { 
    case OVERLAP_INNER_CG_SINGLE_PASS:
      QDPIO::cout << "Using Single Pass Inner Solver" << endl;
      break;
    case OVERLAP_INNER_CG_DOUBLE_PASS:
      QDPIO::cout << "Using Neuberger/Chu Double Pass Inner Solver" << endl;
      break;
    default:
      QDPIO::cerr << "Unknown inner solver type " << endl;
      QDP_abort(1);
    }

    QDPIO::cout << "Number of poles= " << numroot << endl;
    QDPIO::cout << "Overall Factor=  " << coeffP << endl;
    QDPIO::cout << "Numerator coefficients:" << endl;
    for(int n=0; n < numroot; n++) { 
      QDPIO::cout <<"  resP[" << n << "]= " << resP[n] << endl;
    }
    QDPIO::cout << "Denominator roots: " << endl;
    for(int n=0; n < numroot; n++) { 
      QDPIO::cout <<"  rootQ[" << n<< "]= " << rootQ[n] << endl;
    }
  
      
    /* We will also compute the 'function' of the eigenvalues */
    /* for the Wilson vectors to be projected out. */
    if (NEig > 0)
    {
      for(int i = 0; i < NEigVal; i++)
      {
	if (toBool(state.getEvalues()[i] > 0.0))
	  EigValFunc[i] = 1.0;
	else if (toBool(state.getEvalues()[i] < 0.0))
	  EigValFunc[i] = -1.0;
	else
	  EigValFunc[i] = 0.0;
      }
    }

    // Free the arrays allocate by Tony's zolo
    zolotarev_free(rdata);

    QDPIO::cout << "Leaving Init!" << endl << flush;

  }

  //! Produce a linear operator for this action
  /*!
   * The operator acts on the entire lattice
   *
   * \param state_	 gauge field state  	 (Read)
   * \param m_q	         mass for this operator	 (Read)
   */
  UnprecLinearOperator<LatticeFermion,
		       multi1d<LatticeColorMatrix>,
		       multi1d<LatticeColorMatrix> >* 
  OvlapPartFrac4DFermAct::unprecLinOp(Handle< FermState<T,P,Q> > state_, const Real& m_q) const
  {
    START_CODE();

    try {
      const EigenConnectState& state = dynamic_cast<EigenConnectState&>(*state_);
          
      int NEigVal = state.getNEig();
      
      /* The actual number of eigenvectors to project out.
	 The highest of the valid low eigenmodes is not
	 projected out. So we will put NEig = NEigVal - 1 */  
      int NEig = NEigVal;
      
      /* The number of residuals and poles */
      int numroot;
      
      /* The roots, i.e., the shifts in the partial fraction expansion */
      multi1d<Real> rootQ;
      
      /* The residuals in the partial fraction expansion */
      multi1d<Real> resP;
      
      /* This will be our alpha(0) which can be 0 depending on type */
      /* an even- or oddness of RatPolyDeg*/
      Real coeffP; 
      
      /* Array of values of the sign function evaluated on the eigenvectors of H */
      multi1d<Real> EigValFunc(NEigVal);
      
      // Common initialization
      init(numroot, coeffP, resP, rootQ, NEig, EigValFunc, state);
      
      
      /* Finally construct and pack the operator */
      /* This is the operator of the form (1/2)*[(1+mu) + (1-mu)*gamma_5*eps] */
      switch( params.inner_solver_type ) {
      case OVERLAP_INNER_CG_SINGLE_PASS:
	return new lovlapms(*Mact, state_, m_q,
			    numroot, coeffP, resP, rootQ, 
			    NEig, EigValFunc, state.getEvectors(),
			    params.invParamInner.MaxCG, 
			    params.invParamInner.RsdCG, 
			    params.ReorthFreqInner);
	break;
      case OVERLAP_INNER_CG_DOUBLE_PASS:
	return new lovlap_double_pass(*Mact, state_, m_q,
				      numroot, coeffP, resP, rootQ, 
				      NEig, EigValFunc, state.getEvectors(),
				      params.invParamInner.MaxCG, 
				      params.invParamInner.RsdCG, 
				      params.ReorthFreqInner);
	break;
      default:
	QDPIO::cerr << "Unknown OverlapInnerSolverType " << params.inner_solver_type << flush << endl;
	QDP_abort(1);
      }
      
    }
    catch(bad_cast) { 
      QDPIO::cerr << "OverlapPartFrac4DFermAct::unprecLinOp: "
		  << " Failed to downcast ConnectState to OverlapConnectState"
		  << endl;
      QDP_abort(1);
    }

    END_CODE();

    QDPIO::cout << "DANGER!!! About to return zero! " << endl;
    return 0;
  }

  //! Produce a linear operator for this action
  /*!
   * The operator acts on the entire lattice
   *
   * \param state_	 gauge field state  	 (Read)
   */
  LinearOperator<LatticeFermion>* 
  OvlapPartFrac4DFermAct::linOpPrecondition(Handle< FermState<T,P,Q> > state_) const
  {
    START_CODE();
    try {
      const EigenConnectState& state = dynamic_cast<const EigenConnectState&>(*state_);

      int NEigVal = state.getNEig();
      
      /* The actual number of eigenvectors to project out.
	 The highest of the valid low eigenmodes is not
	 projected out. So we will put NEig = NEigVal - 1 */  
      int NEig=NEigVal;
      
      /* The number of residuals and poles */
      int numroot;
      
      /* The roots, i.e., the shifts in the partial fraction expansion */
      multi1d<Real> rootQ;
      
      /* The residuals in the partial fraction expansion */
      multi1d<Real> resP;
      
      /* This will be our alpha(0) which can be 0 depending on type */
      /* an even- or oddness of RatPolyDeg*/
      Real coeffP; 
      
      /* Array of values of the sign function evaluated on the eigenvectors of H */
      multi1d<Real> EigValFunc(NEigVal);
      
      // Common initialization
      initPrec(numroot, coeffP, resP, rootQ, NEig, EigValFunc, state);
      
      
      /* Finally construct and pack the operator */
      /* This is the operator of the form (1/2)*[(1+mu) + (1-mu)*gamma_5*eps] */
      switch( params.inner_solver_type ) {
      case OVERLAP_INNER_CG_SINGLE_PASS:
	return new lovlapms(*Mact, state_, params.Mass,
			    numroot, coeffP, resP, rootQ, 
			    NEig, EigValFunc, state.getEvectors(),
			    params.invParamInner.MaxCG, params.invParamInner.RsdCG, params.ReorthFreqInner);
	break;
      case OVERLAP_INNER_CG_DOUBLE_PASS:
	return new lovlap_double_pass(*Mact, state_, params.Mass,
				      numroot, coeffP, resP, rootQ, 
				      NEig, EigValFunc, state.getEvectors(),
				      params.invParamInner.MaxCG, params.invParamInner.RsdCG, params.ReorthFreqInner);
	break;
      default:
	QDPIO::cerr << "Unknown OverlapInnerSolverType " << params.inner_solver_type << endl;
	QDP_abort(1);
      }
  
    }
    catch(bad_cast) { 
      QDPIO::cerr << "OverlapPartFrac4DFermAct::linOpPrecondition: "
		  << " Failed to downcast ConnectState to OverlapConnectState"
		  << endl;
      QDP_abort(1);
    }

    END_CODE();

    return 0;
  }

  //! Produce a linear operator for gamma5 epsilon(H) psi
  /*!
   * The operator acts on the entire lattice
   *
   * \param state_	 gauge field state  	 (Read)
   */
  LinearOperator<LatticeFermion>* 
  OvlapPartFrac4DFermAct::lgamma5epsH(Handle< FermState<T,P,Q> > state_) const
  {
    START_CODE();

    try {
      const EigenConnectState& state = dynamic_cast<const EigenConnectState&>(*state_);

      int NEigVal = state.getNEig();
      
      /* The actual number of eigenvectors to project out.
       The highest of the valid low eigenmodes is not
       projected out. */  
      int NEig = NEigVal;
      
      /* The number of residuals and poles */
      int numroot;
      
      /* The roots, i.e., the shifts in the partial fraction expansion */
      multi1d<Real> rootQ;
      
      /* The residuals in the partial fraction expansion */
      multi1d<Real> resP;
      
      /* This will be our alpha(0) which can be 0 depending on type */
      /* an even- or oddness of RatPolyDeg*/
      Real coeffP; 
      
      /* Array of values of the sign function evaluated on the eigenvectors of H */
      multi1d<Real> EigValFunc(NEigVal);
      
      // Common initialization
      init(numroot, coeffP, resP, rootQ, NEig, EigValFunc, state);
      
      
      /* Finally construct and pack the operator */
      /* This is the operator of the form (1/2)*[(1+mu) + (1-mu)*gamma_5*eps] */
      switch( params.inner_solver_type ) { 
      case OVERLAP_INNER_CG_SINGLE_PASS:
	return new lg5eps(*Mact, state_,
			  numroot, coeffP, resP, rootQ, 
			  NEig, EigValFunc, state.getEvectors(),
			  params.invParamInner.MaxCG, params.invParamInner.RsdCG, params.ReorthFreqInner);
	break;
      case OVERLAP_INNER_CG_DOUBLE_PASS:
	return new lg5eps_double_pass(*Mact, state_,
				      numroot, coeffP, resP, rootQ, 
				      NEig, EigValFunc, state.getEvectors(),
				      params.invParamInner.MaxCG, params.invParamInner.RsdCG, params.ReorthFreqInner);
	break;
      default:
	QDPIO::cerr << "Unknown OverlapInnerSolverType " << params.inner_solver_type << endl;
	QDP_abort(1);
      }
    }
    catch(bad_cast) { 
      QDPIO::cerr << "OverlapPartFrac4DFermAct::lgamma5epsH: "
		  << " Failed to downcast ConnectState to OverlapConnectState"
		  << endl;
      QDP_abort(1);
    }
      
    END_CODE();
    
    return 0;
  }
  
  //! Produce a linear operator for gamma5 epsilon(H) psi
  /*!
   * The operator acts on the entire lattice
   *
   * \param state_	 gauge field state  	 (Read)
   */
  LinearOperator<LatticeFermion>* 
  OvlapPartFrac4DFermAct::lgamma5epsHPrecondition(Handle< FermState<T,P,Q> > state_) const
  {
    START_CODE();
  
    const EigenConnectState& state = dynamic_cast<EigenConnectState&>(*state_);

    int NEigVal = state.getNEig();

    /* The actual number of eigenvectors to project out.
       The highest of the valid low eigenmodes is not
       projected out. So we will put NEig = NEigVal - 1 */  
    int NEig=NEigVal;

    /* The number of residuals and poles */
    int numroot;

    /* The roots, i.e., the shifts in the partial fraction expansion */
    multi1d<Real> rootQ;

    /* The residuals in the partial fraction expansion */
    multi1d<Real> resP;

    /* This will be our alpha(0) which can be 0 depending on type */
    /* an even- or oddness of RatPolyDeg*/
    Real coeffP; 

    /* Array of values of the sign function evaluated on the eigenvectors of H */
    multi1d<Real> EigValFunc(NEigVal);

    // Common initialization
    initPrec(numroot, coeffP, resP, rootQ, NEig, EigValFunc, state);

  
    /* Finally construct and pack the operator */
    /* This is the operator of the form (1/2)*[(1+mu) + (1-mu)*gamma_5*eps] */
    switch( params.inner_solver_type ) { 
    case OVERLAP_INNER_CG_SINGLE_PASS:
      return new lg5eps(*Mact, state_,
			numroot, coeffP, resP, rootQ, 
			NEig, EigValFunc, state.getEvectors(),
			params.invParamInner.MaxCG, params.invParamInner.RsdCG, params.ReorthFreqInner);
      break;
    case OVERLAP_INNER_CG_DOUBLE_PASS:
      return new lg5eps_double_pass(*Mact, state_,
				    numroot, coeffP, resP, rootQ, 
				    NEig, EigValFunc, state.getEvectors(),
				    params.invParamInner.MaxCG, params.invParamInner.RsdCG, params.ReorthFreqInner);
      break;
    default:
      QDPIO::cerr << "Unknown OverlapInnerSolverType " << params.inner_solver_type << endl;
      QDP_abort(1);
    }
  
    END_CODE();

    return 0;
  }

  //! Produce a conventional MdagM operator for this action
  /*!
   * This is the operator which you can always use
   * The operator acts on the entire lattice
   *
   * \param state_     gauge field state       (Read)
   */
  DiffLinearOperator<LatticeFermion,
		     multi1d<LatticeColorMatrix>,
		     multi1d<LatticeColorMatrix> >* 
  OvlapPartFrac4DFermAct::lMdagM(Handle< FermState<T,P,Q> > state_) const
  {
    // linOp news the linear operator and gives back pointer, 
    // We call lmdagm with this pointer.
    // lmdagm is the only owner
    // No need to grab linOp with handle at this stage.
    return new DiffMdagMLinOp<T,P,Q>( linOp(state_) );
  }
 
  //! Produce a linear operator for this action
  /*!
   * The operator acts on the entire lattice
   *
   * \param state_	 gauge field state  	 (Read)
   */
  DiffLinearOperator<LatticeFermion,
		     multi1d<LatticeColorMatrix>,
		     multi1d<LatticeColorMatrix> >* 
  OvlapPartFrac4DFermAct::lMdagM(Handle< FermState<T,P,Q> > state_, const Chirality& ichiral) const
  {

    // If chirality is none, return traditional MdagM
    if ( ichiral == CH_NONE ) {
      return lMdagM(state_);
    }
    else { 
      const EigenConnectState& state = dynamic_cast<EigenConnectState&>(*state_);
    
      int NEigVal = state.getNEig();

      /* The actual number of eigenvectors to project out.
	 The highest of the valid low eigenmodes is not
	 projected out. So we will put NEig = NEigVal - 1 */  
      int NEig=NEigVal;
    
      /* The number of residuals and poles */
      int numroot;

      /* The roots, i.e., the shifts in the partial fraction expansion */
      multi1d<Real> rootQ;
    
      /* The residuals in the partial fraction expansion */
      multi1d<Real> resP;
    
      /* This will be our alpha(0) which can be 0 depending on type */
      /* an even- or oddness of RatPolyDeg*/
      Real coeffP; 
    
      /* Array of values of the sign function evaluated on the eigenvectors of H */
      multi1d<Real> EigValFunc(NEigVal);
    
      // Common initialization
      init(numroot, coeffP, resP, rootQ, NEig, EigValFunc, state);

  
      // Finally construct and pack the operator 
      // This is the operator of the form (1/2)*[(1+mu) + (1-mu)*gamma_5*eps]
      switch( params.inner_solver_type ) { 
      case OVERLAP_INNER_CG_SINGLE_PASS:
	return new lovddag(*Mact, state_, params.Mass,
			   numroot, coeffP, resP, rootQ, 
			   NEig, EigValFunc, state.getEvectors(),
			   params.invParamInner.MaxCG, params.invParamInner.RsdCG, params.ReorthFreqInner, ichiral);
	break;
      case OVERLAP_INNER_CG_DOUBLE_PASS:
	return new lovddag_double_pass(*Mact, state_, params.Mass,
				       numroot, coeffP, resP, rootQ, 
				       NEig, EigValFunc, state.getEvectors(),
				       params.invParamInner.MaxCG, params.invParamInner.RsdCG, params.ReorthFreqInner, ichiral);
	break;
      default:
	QDPIO::cerr << "Unknown OverlapInnerSolverType " << params.inner_solver_type << endl;
	QDP_abort(1);
      }

    }
    END_CODE();

    return 0;
  }


  //! Produce a linear operator for this action
  /*!
   * The operator acts on the entire lattice
   *
   * \param state_	 gauge field state  	 (Read)
   */
  LinearOperator<LatticeFermion>* 
  OvlapPartFrac4DFermAct::lMdagMPrecondition(Handle< FermState<T,P,Q> > state_, const Chirality& ichiral) const
  {

    // If chirality is none, return traditional MdagM
    if ( ichiral == CH_NONE ) {
      return lMdagM(state_);
    }
    else { 
      const EigenConnectState& state = dynamic_cast<EigenConnectState&>(*state_);
    
      int NEigVal = state.getNEig();

      /* The actual number of eigenvectors to project out.
	 The highest of the valid low eigenmodes is not
	 projected out. So we will put NEig = NEigVal - 1 */  
      int NEig=NEigVal;
    
      /* The number of residuals and poles */
      int numroot;

      /* The roots, i.e., the shifts in the partial fraction expansion */
      multi1d<Real> rootQ;
    
      /* The residuals in the partial fraction expansion */
      multi1d<Real> resP;
    
      /* This will be our alpha(0) which can be 0 depending on type */
      /* an even- or oddness of RatPolyDeg*/
      Real coeffP; 
    
      /* Array of values of the sign function evaluated on the eigenvectors of H */
      multi1d<Real> EigValFunc(NEigVal);
    
      // Common initialization
      initPrec(numroot, coeffP, resP, rootQ, NEig, EigValFunc, state);

  
      // Finally construct and pack the operator 
      // This is the operator of the form (1/2)*[(1+mu) + (1-mu)*gamma_5*eps]
      switch( params.inner_solver_type ) { 
      case OVERLAP_INNER_CG_SINGLE_PASS:
	return new lovddag(*Mact, state_, params.Mass,
			   numroot, coeffP, resP, rootQ, 
			   NEig, EigValFunc, state.getEvectors(),
			   params.invParamInner.MaxCG, params.invParamInner.RsdCG, params.ReorthFreqInner, ichiral);
	break;
      case OVERLAP_INNER_CG_DOUBLE_PASS:
	return new lovddag_double_pass(*Mact, state_, params.Mass,
				       numroot, coeffP, resP, rootQ, 
				       NEig, EigValFunc, state.getEvectors(),
				       params.invParamInner.MaxCG, params.invParamInner.RsdCG, params.ReorthFreqInner, ichiral);
	break;
      default:
	QDPIO::cerr << "Unknown OverlapInnerSolverType " << params.inner_solver_type << endl;
	QDP_abort(1);
      }

    }
    END_CODE();

    return 0;
  }

  //! Produce a conventional MdagM operator for this action
  /*!
   * This is the operator which you can always use
   * The operator acts on the entire lattice
   *
   * \param state_	 gauge field state  	 (Read)
   */
  LinearOperator<LatticeFermion>* 
  OvlapPartFrac4DFermAct::lMdagMPrecondition(Handle< FermState<T,P,Q> > state_) const
  {
    // linOp news the linear operator and gives back pointer, 
    // We call lmdagm with this pointer.
    // lmdagm is the only owner
    // No need to grab linOp with handle at this stage.
    return new MdagMLinOp<LatticeFermion>( linOpPrecondition(state_) );
  }

  //! Create a ConnectState with just the gauge fields
  EigenConnectState*
  OvlapPartFrac4DFermAct::createState(const multi1d<LatticeColorMatrix>& u_) const
  {
    return new EigenConnectState(fbc, u_);
  }


  //! Create OverlapConnectState from XML
  EigenConnectState*
  OvlapPartFrac4DFermAct::createState(const multi1d<LatticeColorMatrix>& u_,
				      XMLReader& state_info_xml,
				      const string& state_info_path) const
  {
    XMLFileWriter test("./test");
    test << state_info_xml;

    //QDPIO::cout << "Creating State from XML: " << reader_contents.str() << endl << flush;
 
    std::string eigen_info_id;
    if( state_info_xml.count("/StateInfo/eigen_info_id") == 1 ) {
     
      read(state_info_xml, "/StateInfo/eigen_info_id", eigen_info_id);
      QDPIO::cout << "Using eigen_info_id: " << eigen_info_id << endl;

      EigenConnectState *ret_val = new EigenConnectState(fbc, u_, eigen_info_id);


      // Check the low evs
      Handle< FermState<T,P,Q>  > state_aux(new SimpleFermState<T,P,Q>(fbc, u_));
      Handle< LinearOperator<T> > Maux = Mact->hermitianLinOp(state_aux);

      for(int vec = 0; vec < ret_val->getNEig(); vec++) { 

	LatticeFermion lambda_e, H_e;
	Real lambda = ret_val->getEvalues()[vec];
	LatticeFermion ev = ret_val->getEvectors()[vec];

	(*Maux)(H_e, ev, PLUS);
	lambda_e = lambda * ev;

	LatticeFermion diff = H_e - lambda_e;
	Double norm_diff = sqrt(norm2(diff));
	Double norm_diff_rel = norm_diff/fabs(lambda);

	QDPIO::cout << "EV Check["<< vec<<"]: lambda="<< lambda <<" resid="<<norm_diff<<" rel. resid="<< norm_diff_rel << endl;

      }
      return ret_val;

    }
    else {
      QDPIO::cout << "No StateInfo Found: " << endl;

      return new EigenConnectState(fbc, u_);
    }
    

  }


}
