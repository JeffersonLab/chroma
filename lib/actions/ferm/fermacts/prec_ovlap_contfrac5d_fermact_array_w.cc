// $Id: prec_ovlap_contfrac5d_fermact_array_w.cc,v 1.2 2004-10-01 06:29:23 bjoo Exp $
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

#include "actions/ferm/fermacts/fermfactory_w.h"
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
    WilsonTypeFermAct< multi1d<LatticeFermion> >* createFermAct(Handle< FermBC< multi1d<LatticeFermion> > > fbc,
						     XMLReader& xml_in,
						     const std::string& path)
    {
      return new EvenOddPrecOvlapContFrac5DFermActArray(fbc, EvenOddPrecOvlapContFrac5DFermActParams(xml_in, path));
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
  EvenOddPrecOvlapContFrac5DFermActArray::EvenOddPrecOvlapContFrac5DFermActArray(
										 Handle< FermBC< multi1d< LatticeFermion> > > fbc_a_, 
										 const EvenOddPrecOvlapContFrac5DFermActParams& params_) :
    fbc(fbc_a_), params(params_) 
  {

    
    // WHAT IS BELOW ONLY WORKS FOR TYPE=0 approximations
    // which is what we use. Forget TYPE=1
    // the Tanh approximation (Higham) is of type TYPE=0
    // We have two cases.
    bool isEvenRatPolyDeg = ( params.RatPolyDeg % 2 == 0);

    if( isEvenRatPolyDeg ) { 
      QDPIO::cout << "Rat Poly Deg is even " << endl;
      N5 = params.RatPolyDeg+1;
      isLastZeroP = true;     
    }
    else { 
      QDPIO::cout << "Rat Poly Deg is odd " << endl;
      N5 = params.RatPolyDeg;
      isLastZeroP = false;
    }

    QDPIO::cout << "N5 is " << N5 << endl;
    QDPIO::cout << "Last beta coefficient is zero: " << isLastZeroP << endl;
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

    for(int i=0; i < N5; i++) { 
      QDPIO::cout << "beta["<<i<<"] = " << beta[i] << endl;
    }

    alpha.resize(N5);
    for(int i = 0; i < N5; i++) {
      alpha[i] = Real(1);
    } 
    alpha[N5-1] = rdata->beta[N5-1];

    // For the moment choose gamma = 1/sqrt(beta) */
    // except for gamma(N5-1) which always has to be set to 1 */
    // The N5-1 case is special anyway as in certain cases rdata->beta[i]
    // is strictly zero
    multi1d<Real> gamma(N5);
    for(int i=0; i < N5-1; i++) { 
      gamma[i] = Real(1)/ sqrt( rdata->beta[i] );
    }
    gamma[N5-1] = Real(1);

    // Now perform the equivalence transformation on the off diagonal
    // elements 
    for(int i=0; i < N5; i++) {
      beta[i] = beta[i]*gamma[i]*gamma[i];
    }
    
    // and the off diagonal ones
    for(int i=0; i < N5-1; i++) {
      alpha[i] = alpha[i]*gamma[i]*gamma[i+1];
    }
    
    
    QDPIO::cout << "EvenOddPrecOvlapContfrac5d: " 
		<< " N5=" << N5 << " scale=" << scale_fac
		<< " Mass=" << params.Mass 
		<< " OverMass=" << params.OverMass << endl;

    QDPIO::cout << "Approximation on [-1,eps] U [eps,1] with eps = " << eps <<endl;
    
    QDPIO::cout << "Maximum error | R(x) - sgn(x) | <= Delta = " << maxerr << endl;
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
    free( rdata->a );
    free( rdata->ap );
    free( rdata->alpha );
    free( rdata->beta );
    free( rdata );
  }

  //! Produce a linear operator for this action
  /*!
   * The operator acts on the entire lattice
   *
   * \param state	    gauge field     	       (Read)
   */
  const EvenOddPrecLinearOperator<multi1d<LatticeFermion> >* 
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
    multi1d<LatticeFermion> tmp5_1(N5);
    multi1d<LatticeFermion> tmp5_2(N5);
    // Construct the linear operator
    Handle<const EvenOddPrecLinearOperator< multi1d<LatticeFermion> > > A(linOp(state));


    // Construct the source: 
    for(int i=0; i < N5; i++) { 
      chi5[i] = zero;
      psi5[i] = zero;
      tmp5_1[i] = zero;
      tmp5_2[i] = zero;
    }
      
    // Set both even and odd components of psi5[N5-1] = gamma_5 chi 
    // -- use psi as temporary here.
    psi5[N5-1] = Gamma(G5)*chi;
    
    
    // Chi5[N5-1]_e = gamma_5 chi_e = psi5[N5-1]_e
    chi5[N5-1][rb[0]] = psi5[N5-1];

    // Chi5[N5-1]_o = S_o - Q_oe Q_ee^{-1} S_e 
    //              = gamma_5 chi_o - Q_oe Q_ee^{-1} gamma_5 chi_e
    //              = psi5[N5-1]_o - Q_oe Q_ee^{-1} psi5[N5-1]_e
    
    // Get  tmp5_2[N5-1]_o = Q_oe Q_ee^{-1} psi5[N5-1]_e
    A->evenEvenInvLinOp(tmp5_1, psi5, PLUS);
    A->oddEvenLinOp(tmp5_2, tmp5_1, PLUS);

    // Now subtract psi5[N5-1]_o - tmp5_2[N5-1]_o
    for(int i=0; i < N5; i++) { 
      chi5[i][rb[1]] = psi5[i] - tmp5_2[i];
    }
    QDPIO::cout << "|| chi5 || after source prec = " << sqrt(norm2(chi5[N5-1])) << endl;      
    
    // Set psi[N5-1] = psi
    for(int i=0; i < N5-1; i++) { 
      psi5[i]=zero;
    }
    psi5[N5-1] = psi;     
  
    switch(invParam.invType) {
    case CG_INVERTER: 
      {

	// Solve M^{+}M psi = M^{+} chi
        (*A)(tmp5_1, chi5, MINUS);
	QDPIO::cout << "|| chi5 || after source prec2 = " << sqrt(norm2(tmp5_1)) << endl;
	QDPIO::cout << "|| psi5 || after source prec2 = " << sqrt(norm2(psi5)) << endl;

	InvCG2(*A, tmp5_1, psi5, invParam.RsdCG, invParam.MaxCG, n_count);

	// psi[N5-1]_odd now holds the desired piece.

	// Reconstruct psi[N5-1]_e = Q_ee^{-1} S_e - Q_ee^{-1}Q_eo psi[N5-1]_o
	//                        
	A->evenOddLinOp(tmp5_1, psi5, PLUS);
	for(int i=0; i < N5; i++) { 
	  tmp5_2[i][rb[0]] = chi5[i] - tmp5_1[i];
	}
	// This had better leave the odd bits untouched.
	A->evenEvenInvLinOp(psi5, tmp5_2, PLUS);
	break;
      }
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
    // normalisatoin
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
