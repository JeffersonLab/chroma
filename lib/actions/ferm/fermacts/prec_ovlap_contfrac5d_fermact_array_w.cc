// $Id: prec_ovlap_contfrac5d_fermact_array_w.cc,v 1.7 2004-12-24 04:23:20 edwards Exp $
/*! \file
 *  \brief Unpreconditioned extended-Overlap (5D) (Naryanan&Neuberger) action
 */

#include "chromabase.h"

#include "actions/ferm/fermacts/unprec_wilson_fermact_w.h"
#include "actions/ferm/fermacts/prec_ovlap_contfrac5d_fermact_array_w.h"
#include "actions/ferm/linop/prec_ovlap_contfrac5d_linop_array_w.h"
#include "actions/ferm/linop/lmdagm.h"
#include "actions/ferm/invert/invcg2_array.h"
#include "zolotarev.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermbcs/fermbcs_w.h"

#include "io/enum_io/enum_io.h"
#include "io/overlap_state_info.h"

using namespace QDP;
using namespace Chroma;
namespace Chroma
{

  //! Hooks to register the class with the fermact factory
  namespace EvenOddPrecOvlapContFrac5DFermActArrayEnv
  {
    //! Callback function
    WilsonTypeFermAct< multi1d<LatticeFermion> >* createFermAct(XMLReader& xml_in,
								const std::string& path)
    {
      return new EvenOddPrecOvlapContFrac5DFermActArray(WilsonTypeFermBCArrayEnv::reader(xml_in, path), 
							EvenOddPrecOvlapContFrac5DFermActParams(xml_in, path));
    }
    
    //! Name to be used
    const std::string name = "OVERLAP_CONTINUED_FRACTION_5D";
    
    //! Register the Wilson fermact
    const bool registered = TheWilsonTypeFermActArrayFactory::Instance().registerObject(name, createFermAct);
  } // End Namespace EvenOddPrecOvlapContFrac5DFermActArrayEnv

  
  //! Read XML
  EvenOddPrecOvlapContFrac5DFermActParams::EvenOddPrecOvlapContFrac5DFermActParams(XMLReader& xml, const std::string& path)
  {
    XMLReader in(xml, path);
    
    try {
      
      read(in, "Mass", Mass);
      read(in, "RatPolyDeg", RatPolyDeg);
      
      if( in.count("ApproximationType") == 1 ) { 
	read(in, "ApproximationType", approximation_type);
      }
      else { 
	// Default coeffs are Zolotarev
	approximation_type = COEFF_TYPE_ZOLOTAREV;
      }
      read(in, "OverMass", OverMass);
    }
    catch( const string &e ) {
      QDPIO::cerr << "Caught Exception reading Zolo5D Fermact params: " << e << endl;
      QDP_abort(1);
    }
  }

  void read(XMLReader& xml_in, const string& path,
	    EvenOddPrecOvlapContFrac5DFermActParams& param) 
  {
    
    EvenOddPrecOvlapContFrac5DFermActParams tmp(xml_in, path);
    param = tmp;
  }
  
  void write(XMLWriter& xml_out, const string& path, const EvenOddPrecOvlapContFrac5DFermActParams& p)
  {
    if ( path != "." ) { 
      push( xml_out, path);
    }
    
    write(xml_out, "Mass", p.Mass);
    write(xml_out, "OverMass", p.OverMass);
    write(xml_out, "RatPolyDeg", p.RatPolyDeg);
    write(xml_out, "ApproximationType", p.approximation_type);
    
    pop(xml_out);
    
    
    if( path != "." ) { 
      pop(xml_out);
    }
  }

  
  
  // Construct the action out of a parameter structure
  EvenOddPrecOvlapContFrac5DFermActArray::EvenOddPrecOvlapContFrac5DFermActArray(Handle< FermBC< multi1d< LatticeFermion> > > fbc_a_, 
 const EvenOddPrecOvlapContFrac5DFermActParams& params_) :
    fbc(fbc_a_), params(params_) 
  {

    
    // WHAT IS BELOW ONLY WORKS FOR TYPE=0 approximations
    // which is what we use. Forget TYPE=1
    // the Tanh approximation (Higham) is of type TYPE=0
    // We have two cases.
    bool isEvenRatPolyDeg = ( params.RatPolyDeg % 2 == 0);

    if( isEvenRatPolyDeg ) { 
      N5 = params.RatPolyDeg+1;
      isLastZeroP = true;     
    }
    else { 
      N5 = params.RatPolyDeg;
      isLastZeroP = false;
    }
  }


  void
  EvenOddPrecOvlapContFrac5DFermActArray::init(Real& scale_fac,
					  multi1d<Real>& alpha,
					  multi1d<Real>& beta,
					  const OverlapConnectState& state) const
  {
  
    int type = 0;
    zolotarev_data *rdata;
    Real eps;

    switch(params.approximation_type) { 
    case COEFF_TYPE_ZOLOTAREV:
      scale_fac = Real(1) / state.getApproxMax();
      eps = state.getApproxMin() * scale_fac;

      QDPIO::cout << "Initing Linop with Zolotarev Coefficients" << endl;
      rdata = zolotarev(toFloat(eps), params.RatPolyDeg, type);    
      break;

    case COEFF_TYPE_TANH:
      scale_fac = Real(1) / state.getApproxMax();
      eps = state.getApproxMin() * scale_fac;

      QDPIO::cout << "Initing Linop with Higham Rep tanh Coefficients" << endl;
      rdata = higham(toFloat(eps), params.RatPolyDeg);

      break;
    case COEFF_TYPE_TANH_UNSCALED:
      scale_fac = Real(1);
      eps = state.getApproxMin();

      QDPIO::cout << "Initing Linop with Unscaled Higham Rep tanh Coefficients" << endl;
      rdata = higham(toFloat(eps), params.RatPolyDeg);

      break;

    default:
      // The map system should ensure that we never get here but 
      // just for style
      QDPIO::cerr << "Unknown coefficient type: " << params.approximation_type
		  << endl;
      QDP_abort(1);
    }
    
    Real maxerr = (Real)(rdata->Delta);

    // Check N5 is good:
    if( N5 != rdata->db ) { 
      QDPIO::cerr << "Mismatch between N5 and N5 from Coefficient Code" << endl;
      QDPIO::cerr << "N5 = " << N5 << " rdata->db=" << rdata->db << endl;
      QDP_abort(1);
    }

    
    // The coefficients from the continued fraction
    beta.resize(N5);
    for(int i = 0; i < N5; i++) { 
      beta[i] = rdata->beta[i];
    }

    alpha.resize(N5-1);
    for(int i=0; i < N5-1; i++) { 
      alpha[i] = Real(1);
    }

    // The gamma's are the equivalence transforms
    // There are N5-1 of them and they appear in the 
    // diagonal terms as gamma^2 
    // and in the off diagonal terms as gamma_i gamma_i+1
    // except in the last one which is just gamma
    //
    // For the moment choose gamma_i = 1/sqrt(beta_i) */
    
    multi1d<Real> gamma(N5-1);
    for(int i=0; i < N5-1; i++) { 
      gamma[i] = Real(1)/ sqrt( beta[i] );
    }

    // Now perform the equivalence transformation 
    //
    // On the diagonal coefficients
    // Note that beta[N5-1] is NOT changed
    for(int i=0; i < N5-1; i++) {
      beta[i] = beta[i]*gamma[i]*gamma[i];
    }
    
    // and the off diagonal ones
    // from 0..N5-3 we have gamma_i gamma_i+1
    // and on N5-2 we have gamma_i 
    for(int i=0; i < N5-2; i++) {
      alpha[i] *= gamma[i]*gamma[i+1];
    }
    alpha[N5-2] *= gamma[N5-2];
    
    QDPIO::cout << "EvenOddPrecOvlapContfrac5d: " << endl
                << "Degree="<< params.RatPolyDeg
		<< "N5=" << N5 << " scale=" << scale_fac
		<< "Mass=" << params.Mass << endl
		<< "OverMass=" << params.OverMass 
		<< "IsLastZeroP=" << isLastZeroP << endl;

    QDPIO::cout << "Approximation on [-1,eps] U [eps,1] with eps = " << eps <<endl;
    
    QDPIO::cout << "Maximum error | R(x) - sgn(x) | <= Delta = " << maxerr << endl;
    /*
    for(int i=0; i < N5; i++) { 
      QDPIO::cout << "beta["<<i<<"] = " << beta[i] << endl;
    }
    for(int i=0; i < N5; i++) { 
      QDPIO::cout << "alpha["<<i<<"] = " << alpha[i] << endl;
    }
    */

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
      QDP_abort(1);
      break;
    }

    // free the arrays allocated by Tony's Zolo
    zolotarev_free(rdata);
  }

  //! Produce a linear operator for this action
  /*!
   * The operator acts on the entire lattice
   *
   * \param state	    gauge field     	       (Read)
   */
  const EvenOddPrecLinearOperator< multi1d<LatticeFermion>, multi1d<LatticeColorMatrix> >* 
  EvenOddPrecOvlapContFrac5DFermActArray::linOp(Handle<const ConnectState> state_) const
  {
    START_CODE();
    try { 
      
      // This throws a bad cast exception if the cast fails
      // Hence the "try" above
      const OverlapConnectState& state = 
	dynamic_cast<const OverlapConnectState&>(*state_);
      
      multi1d<Real> alpha;
      multi1d<Real> beta;
      Real scale_factor;
      
      init(scale_factor, alpha, beta, state);
      
      return new EvenOddPrecOvlapContFrac5DLinOpArray(state_,
						      params.Mass,
						      params.OverMass,
						      N5,
						      scale_factor,
						      alpha,
						      beta,
						      isLastZeroP);
      
      
      
    }
    catch( bad_cast ) { 
      QDPIO::cerr << "EvenOddPrecOvlapContFrac5DFermActArray::linOp(): ";
      QDPIO::cerr << "Couldnt cast ConnectState to OverlapConnectState " 
		  << endl;
      QDP_abort(1);
    }
    
    // Should never get here... Just to satisfy type system
    return 0;
  }
  
  //! Produce a M^dag.M linear operator for this action
  /*!
   * The operator acts on the entire lattice *
   * \param state	    gauge field     	       (Read)
   */
  const LinearOperator<multi1d<LatticeFermion> >* 
  EvenOddPrecOvlapContFrac5DFermActArray::lMdagM(Handle<const ConnectState> state) const
  {
    return new lmdagm<multi1d<LatticeFermion> >(linOp(state));
  }

  

  //! Propagator of an un-preconditioned Extended-Overlap linear operator
  /*!
   * \param psi      quark propagator ( Modify )
   * \param state    gauge field ( Read )
   * \param chi      source ( Read )
   * \param invType  inverter type ( Read (
   * \param RsdCG    CG (or MR) residual used here ( Read )
   * \param MaxCG    maximum number of CG iterations ( Read )
   * \param ncg_had  number of CG iterations ( Write )
   */
  void 
  EvenOddPrecOvlapContFrac5DFermActArray::qprop(LatticeFermion& psi, 
						Handle<const ConnectState> state, 
						const LatticeFermion& chi, 
						const InvertParam_t& invParam,
						int& ncg_had) const
  {
    
    START_CODE();
    
    const Real Mass = quark_mass();
    int n_count;
    
    int G5 = Ns*Ns - 1;
    
    // Initialize the 5D fields
    multi1d<LatticeFermion> chi5(N5);
    multi1d<LatticeFermion> psi5(N5);
    // Construct the linear operator
    Handle<const EvenOddPrecLinearOperatorBase< multi1d<LatticeFermion> > > A(linOp(state));


  
    switch(invParam.invType) {
    case CG_INVERTER: 
      {
	multi1d<LatticeFermion> tmp5_1(N5);
	
	{
	  multi1d<LatticeFermion> tmp5_2(N5);
	  multi1d<LatticeFermion> tmp5_3(N5);

	  chi5 = zero;
	  psi5 = zero;
	  tmp5_1 = zero;
		
	  // We need to solve D_5 psi = (0,0,0,...,gamma_5 chi)^T
	  // Use both subsets
	  tmp5_1[N5-1] = Gamma(G5)*chi;
	
	
	  // tmp5_3_odd = Qoe Qee^{-1} S_e
	  A->evenEvenInvLinOp(tmp5_2, tmp5_1, PLUS);
	  A->oddEvenLinOp(tmp5_3, tmp5_2, PLUS);


	  // chi5_e = S_e
	  // chi5_o = S_o - QoeQee^{-1} S_e
	  for(int i=0; i < N5; i++) { 
	    chi5[i][rb[0]] = tmp5_1[i];
	    chi5[i][rb[1]] = tmp5_1[i] - tmp5_3[i];
	  }
	}  // tmp5_2 and tmp5_3 go away

	psi5[N5-1][rb[1]] = psi;     
	(*A)(tmp5_1, chi5, MINUS);

	// Solve M^{+}M psi = M^{+} chi
	InvCG2(*A, tmp5_1, psi5, invParam.RsdCG, invParam.MaxCG, n_count);
        
	
	// psi[N5-1]_odd now holds the desired piece.

	// Reconstruct psi[N5-1]_e = Q_ee^{-1} S_e - Q_ee^{-1}Q_eo psi[N5-1]_o
	//         = Q_ee^{-1} ( S_e - Q_eo psi_o )
	{ 

	  // Dont need to initialise as the parts I use 
	  // will be over written the other parts I ignore
	  multi1d<LatticeFermion> tmp5_2(N5);
	  multi1d<LatticeFermion> tmp5_3(N5);

	  // tmp5_2_e = Qeo psi_o
	  A->evenOddLinOp(tmp5_2, psi5, PLUS);
	  for(int i=0; i < N5; i++) {

	    // tmp5_3_e = S_e - Qeo psi_o 
	    //          =     - tmp5_2_e
	    tmp5_3[i][rb[0]] = chi5[i] - tmp5_2[i];

	  }

	  // tmp5_1_e = Qee^{-1} ( S_e - Qeo psi_o )
	  A->evenEvenInvLinOp(tmp5_1, tmp5_3, PLUS);

	  // psi5_e = tmp5_1
	  for(int i=0; i < N5; i++) {
	    psi5[i][rb[0]] = tmp5_1[i];
	  }
	} // tmp5_2, tmp5_3 disappear 
      } // tmp5_1 disappears
      break;
      
    case MR_INVERTER:
      QDP_error_exit("Unsupported inverter type", invParam.invType);
      break;

    case BICG_INVERTER:
      QDP_error_exit("Unsupported inverter type", invParam.invType);
      break;
      
    default:
      QDP_error_exit("Unknown inverter type", invParam.invType);
    }
  
    if ( n_count == invParam.MaxCG )
      QDP_error_exit("no convergence in the inverter", n_count);
    
    ncg_had = n_count;
    
    // Solution now lives in chi5
    
    // Multiply back in factor 2/(1-m) to return to 
    // (1/2)( 1 + m + (1-m) gamma_5 epsilon  )
    // normalisation
    psi5[N5-1] *= Real(2)/(Real(1)-Mass);
    
    // Remove contact term
    psi = psi5[N5-1] - chi;
    
    // Overall normalization
    Real ftmp1 = Real(1) / (Real(1) - Mass);
    psi *= ftmp1;

    END_CODE();
  }
  
  //! Create a ConnectState with just the gauge fields
  const OverlapConnectState*
  EvenOddPrecOvlapContFrac5DFermActArray::createState(const multi1d<LatticeColorMatrix>& u_) const
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
  EvenOddPrecOvlapContFrac5DFermActArray::createState(const multi1d<LatticeColorMatrix>& u_,
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
  EvenOddPrecOvlapContFrac5DFermActArray::createState(const multi1d<LatticeColorMatrix>& u_,
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
  EvenOddPrecOvlapContFrac5DFermActArray::createState(const multi1d<LatticeColorMatrix>& u_,
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
  EvenOddPrecOvlapContFrac5DFermActArray::createState(const multi1d<LatticeColorMatrix>& u_,
				      XMLReader& state_info_xml,
				      const string& state_info_path) const
  {
    multi1d<LatticeColorMatrix> u_tmp = u_;
    
    // HACK UP A LINEAR OPERATOR TO CHECK EIGENVALUES/VECTORS WITH
    getFermBC().modifyU(u_tmp);

    Handle< FermBC<LatticeFermion> > fbc4 = new PeriodicFermBC<LatticeFermion>();
    UnprecWilsonFermAct S_w(fbc4, params.OverMass);

    Handle< const ConnectState > state_aux = new SimpleConnectState(u_tmp);
    Handle< const LinearOperator<LatticeFermion> > Maux = 
      S_w.gamma5HermLinOp(state_aux);
    
    
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
  EvenOddPrecOvlapContFrac5DFermActArray::createState(const multi1d<LatticeColorMatrix>& u_,
				      const OverlapStateInfo& state_info) const
  {
    // HACK UP A LINEAR OPERATOR TO CHECK EIGENVALUES/VECTORS WITH
    multi1d<LatticeColorMatrix> u_tmp = u_;
    getFermBC().modifyU(u_tmp);

    Handle< FermBC<LatticeFermion> > fbc4 = new PeriodicFermBC<LatticeFermion>();
    UnprecWilsonFermAct S_w(fbc4, params.OverMass);

    Handle< const ConnectState > state_aux = new SimpleConnectState(u_tmp);
    Handle< const LinearOperator<LatticeFermion> > Maux = 
      S_w.gamma5HermLinOp(state_aux);
    
    
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

} // End namespace Chroma
