// $Id: prec_ht_contfrac5d_fermact_array_w.cc,v 2.2 2006-01-17 16:01:46 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned extended-Overlap (5D) (Naryanan&Neuberger) action
 */

#include "chromabase.h"

#include "actions/ferm/fermacts/unprec_wilson_fermact_w.h"
#include "actions/ferm/fermacts/prec_ht_contfrac5d_fermact_array_w.h"
#include "actions/ferm/linop/unprec_wilson_linop_w.h"
#include "actions/ferm/linop/unprec_dwftransf_linop_w.h"
#include "actions/ferm/linop/prec_ht_contfrac5d_linop_array_w.h"
// #include "actions/ferm/linop/prec_ovlap_contfrac5d_pv_linop_array_w.h"
#include "actions/ferm/invert/invcg2_array.h"
#include "zolotarev.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermbcs/fermbcs_w.h"

#include "io/enum_io/enum_io.h"
#include "io/overlap_state_info.h"

namespace Chroma
{

  //! Hooks to register the class with the fermact factory
  namespace EvenOddPrecHtContFrac5DFermActArrayEnv
  {
    //! Callback function
    WilsonTypeFermAct5D< LatticeFermion, multi1d<LatticeColorMatrix> >* createFermAct5D(XMLReader& xml_in,
											const std::string& path)
    {
      return new EvenOddPrecHtContFrac5DFermActArray(WilsonTypeFermBCArrayEnv::reader(xml_in, path), 
							EvenOddPrecHtContFrac5DFermActParams(xml_in, path));
    }

    //! Callback function
    /*! Differs in return type */
    FermionAction<LatticeFermion>* createFermAct(XMLReader& xml_in,
						 const std::string& path)
    {
      return createFermAct5D(xml_in, path);
    }

    //! Name to be used
    const std::string name = "HT_CONTINUED_FRACTION_5D";
    
    //! Register all the factories
    bool registerAll()
    {
      return Chroma::TheFermionActionFactory::Instance().registerObject(name, createFermAct)
	& Chroma::TheWilsonTypeFermAct5DFactory::Instance().registerObject(name, createFermAct5D);
    }

    //! Register the fermact
    const bool registered = registerAll();
  } // End Namespace EvenOddPrecHtContFrac5DFermActArrayEnv

  
  //! Read XML
  EvenOddPrecHtContFrac5DFermActParams::EvenOddPrecHtContFrac5DFermActParams(XMLReader& xml, const std::string& path)
  {
    XMLReader in(xml, path);
    
    try 
    {
      read(in, "Mass", Mass);
      read(in, "RatPolyDeg", RatPolyDeg);
      read(in, "OverMass", OverMass);

      if( in.count("ApproximationType") == 1 ) 
      { 
      	read(in, "ApproximationType", approximation_type);
      }
      else 
      {
	// Default coeffs are unscaled tanh
	approximation_type = COEFF_TYPE_TANH_UNSCALED;
      }

      if (approximation_type == COEFF_TYPE_ZOLOTAREV)
      {
	read(in, "ApproxMin", ApproxMin);
	read(in, "ApproxMax", ApproxMax);
      }
      else
      {
	ApproxMin = ApproxMax = 0.0;
      }

      int count_b5 = in.count("b5");
      int count_c5 = in.count("c5");

      if( count_b5 == 0 && count_c5 == 0 ) 
      {
	QDPIO::cout << "b5 and c5 not specified. Using Shamir values: b5=1 c5=0" << endl;
	b5 = Real(1);
	c5 = Real(0);
      }
      else {
	read(in, "b5", b5);
	read(in, "c5", c5);
      }
    }
    catch( const string &e ) {
      QDPIO::cerr << "Caught Exception reading HT ContFrac Fermact params: " << e << endl;
      QDP_abort(1);
    }
  }

  void read(XMLReader& xml_in, const string& path,
	    EvenOddPrecHtContFrac5DFermActParams& param) 
  {
    
    EvenOddPrecHtContFrac5DFermActParams tmp(xml_in, path);
    param = tmp;
  }
  
  void write(XMLWriter& xml_out, const string& path, const EvenOddPrecHtContFrac5DFermActParams& p)
  {
    if ( path != "." ) { 
      push( xml_out, path);
    }
    
    write(xml_out, "Mass", p.Mass);
    write(xml_out, "OverMass", p.OverMass);
    write(xml_out, "RatPolyDeg", p.RatPolyDeg);
    write(xml_out, "ApproximationType", p.approximation_type);
    write(xml_out, "b5", p.b5);
    write(xml_out, "c5", p.c5);
    write(xml_out, "ApproxMin", p.ApproxMin);
    write(xml_out, "ApproxMax", p.ApproxMax);

    pop(xml_out);
    
    
    if( path != "." ) { 
      pop(xml_out);
    }
  }

  
  
  // Construct the action out of a parameter structure
  EvenOddPrecHtContFrac5DFermActArray::EvenOddPrecHtContFrac5DFermActArray(Handle< FermBC< multi1d< LatticeFermion> > > fbc_a_, 
										 const EvenOddPrecHtContFrac5DFermActParams& params_) :
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
  EvenOddPrecHtContFrac5DFermActArray::init(Real& scale_fac,
					    multi1d<Real>& alpha,
					    multi1d<Real>& beta,
					    const OverlapConnectState& state) const
  {
    int type = 0;
    zolotarev_data *rdata;
    Real epsilon;

    Real approxMin = (state.getEigVal().size() != 0) ? state.getApproxMin() : params.ApproxMin;
    Real approxMax = (state.getEigVal().size() != 0) ? state.getApproxMax() : params.ApproxMax;

    switch(params.approximation_type) 
    {
    case COEFF_TYPE_ZOLOTAREV:
      epsilon = approxMin / approxMax;
      QDPIO::cout << "Initing Linop with Zolotarev Coefficients: epsilon = " << epsilon << endl;
      rdata = zolotarev(toFloat(epsilon), params.RatPolyDeg, type);    
      scale_fac = Real(1) / approxMax;
      break;

    case COEFF_TYPE_TANH_UNSCALED:
      epsilon = approxMin;
      QDPIO::cout << "Initing Linop with Higham Rep tanh Coefficients" << endl;
      rdata = higham(toFloat(epsilon), params.RatPolyDeg);
      scale_fac = Real(1);
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
    
    QDPIO::cout << "EvenOddPrecHtContfrac5d: " << endl
                << "  Degree="<< params.RatPolyDeg
		<< "  N5=" << N5 << " scale=" << scale_fac
		<< "  Mass=" << params.Mass << endl
		<< "  OverMass=" << params.OverMass 
		<< "  IsLastZeroP=" << isLastZeroP << endl;

    QDPIO::cout << "Approximation on [-1,eps] U [eps,1] with eps = " << epsilon <<endl;
    
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
   * \param state	    gauge field     	       (Read)
   */
  const EvenOddPrecConstDetLinearOperator< multi1d<LatticeFermion>, multi1d<LatticeColorMatrix> >* 
  EvenOddPrecHtContFrac5DFermActArray::linOp(Handle<const ConnectState> state_) const
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
      
      return new EvenOddPrecHtContFrac5DLinOpArray(state_,
						   params.Mass,
						   params.OverMass,
						   N5,
						   scale_factor,
						   alpha,
						   beta,
						   params.b5,
						   params.c5,
						   isLastZeroP);
    }
    catch( bad_cast ) { 
      QDPIO::cerr << "EvenOddPrecHtContFrac5DFermActArray::linOp(): ";
      QDPIO::cerr << "Couldnt cast ConnectState to OverlapConnectState " 
		  << endl;
      QDP_abort(1);
    }
    
    // Should never get here... Just to satisfy type system
    return 0;
  }
  

  //! Produce a Pauli-Villars linear operator for this action
  /*!
   * \param state	    gauge field     	       (Read)
   */
  const EvenOddPrecConstDetLinearOperator< multi1d<LatticeFermion>, multi1d<LatticeColorMatrix> >* 
  EvenOddPrecHtContFrac5DFermActArray::linOpPV(Handle<const ConnectState> state_) const
  {
#if 0
    try { 
      
      // This throws a bad cast exception if the cast fails
      // Hence the "try" above
      const OverlapConnectState& state = 
	dynamic_cast<const OverlapConnectState&>(*state_);
      
      multi1d<Real> alpha;
      multi1d<Real> beta;
      Real scale_factor;
      
      init(scale_factor, alpha, beta, state);
      
      // Hmm, not sure about what all the rescaling does to the PV....
      return new EvenOddPrecHtContFrac5DPVLinOpArray(state_,
						     params.Mass,
						     params.OverMass,
						     N5,
						     scale_factor,
						     alpha,
						     beta,
						     isLastZeroP);
    }
    catch( bad_cast ) { 
      QDPIO::cerr << "EvenOddPrecHtContFrac5DFermActArray::linOp(): ";
      QDPIO::cerr << "Couldnt cast ConnectState to OverlapConnectState " 
		  << endl;
      QDP_abort(1);
    }
#endif

    QDPIO::cerr << "Not yet implemented PV stuff for prec_ht_contfrac" << endl;
    QDP_abort(1);

    // Should never get here... Just to satisfy type system
    return 0;
  }
  
  

  //! Propagator of an un-preconditioned Extended-Overlap linear operator
  /*! \ingroup qprop
   *
   * Propagator solver for DWF-like fermions
   */
  template<typename T, typename P>
  class HtContFrac5DQprop : public SystemSolver<T>
  {
  public:
    //! Constructor
    /*!
     * \param A_         Linear operator ( Read )
     * \param Mass_      quark mass ( Read )
     */
    HtContFrac5DQprop(Handle< const EvenOddPrecConstDetLinearOperator<multi1d<T>, P> > A_,
		      Handle< const LinearOperator<T> > D_denum_,
		      const Real& Mass_,
		      const InvertParam_t& invParam_) : 
      A(A_), D_denum(D_denum_), Mass(Mass_),  invParam(invParam_) {}

    //! Destructor is automatic
    ~HtContFrac5DQprop() {}

    //! Return the subset on which the operator acts
    const OrderedSubset& subset() const {return all;}

    //! Solver the linear system
    /*!
     * \param psi      quark propagator ( Modify )
     * \param chi      source ( Read )
     * \return number of CG iterations
     */
    int operator() (T& psi, const T& chi) const
    {
      START_CODE();
    
      const int N5 = A->size();
      int n_count;
    
      int G5 = Ns*Ns - 1;
      
      // Initialize the 5D fields
      multi1d<T> chi5(N5);
      multi1d<T> psi5(N5);

      switch(invParam.invType) 
      {
      case CG_INVERTER: 
      {
	multi1d<T> tmp5_1(N5);
	{
	  multi1d<T> tmp5_2(N5);
	  multi1d<T> tmp5_3(N5);

	  chi5 = zero;
	  psi5 = zero;
	  tmp5_1 = zero;

	  // Need to prepare the source 
	  {
	    LatticeFermion tmp4 = Gamma(G5)*chi;
	    (*D_denum)(tmp5_1[N5-1], tmp4, MINUS);
	  }

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
	  multi1d<T> tmp5_2(N5);
	  multi1d<T> tmp5_3(N5);

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

      return n_count;
    }

  private:
    // Hide default constructor
    HtContFrac5DQprop() {}

    Handle< const EvenOddPrecConstDetLinearOperator<multi1d<T>,P> > A;
    Handle< const LinearOperator<T> > D_denum;
    const Real Mass;
    const InvertParam_t invParam;
  };

 
  //! Propagator of an un-preconditioned Extended-Overlap linear operator
  const SystemSolver<LatticeFermion>* 
  EvenOddPrecHtContFrac5DFermActArray::qprop(Handle<const ConnectState> state,
					     const InvertParam_t& invParam) const
  {
    const multi1d<LatticeColorMatrix>& u = state->getLinks();
    Real a5 = params.b5- params.c5;
    Real WilsonMass = -params.OverMass ;
    Handle<LinearOperator<LatticeFermion> > D_w(new UnprecWilsonLinOp(u, WilsonMass));
    Handle< const LinearOperator<LatticeFermion> > D_denum(new UnprecDWFTransfDenLinOp(u, a5, D_w));

    return new HtContFrac5DQprop<LatticeFermion, multi1d<LatticeColorMatrix> >(Handle< const EvenOddPrecConstDetLinearOperator< multi1d<LatticeFermion> , multi1d<LatticeColorMatrix> > >(linOp(state)),
									       D_denum, 
									       getQuarkMass(),
									       invParam);
  }
  

  //! Create a ConnectState with just the gauge fields
  const OverlapConnectState*
  EvenOddPrecHtContFrac5DFermActArray::createState(const multi1d<LatticeColorMatrix>& u_) const
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
  EvenOddPrecHtContFrac5DFermActArray::createState(const multi1d<LatticeColorMatrix>& u_,
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
  EvenOddPrecHtContFrac5DFermActArray::createState(const multi1d<LatticeColorMatrix>& u_,
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
  EvenOddPrecHtContFrac5DFermActArray::createState(const multi1d<LatticeColorMatrix>& u_,
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
  EvenOddPrecHtContFrac5DFermActArray::createState(const multi1d<LatticeColorMatrix>& u_,
						   XMLReader& state_info_xml,
						   const string& state_info_path) const
  {
    multi1d<LatticeColorMatrix> u_tmp = u_;
    
    // HACK UP A LINEAR OPERATOR TO CHECK EIGENVALUES/VECTORS WITH
    getFermBC().modifyU(u_tmp);

    Handle< FermBC<LatticeFermion> > fbc4 = new PeriodicFermBC<LatticeFermion>();
    UnprecWilsonFermAct S_w(fbc4, Real(-params.OverMass));

    Handle< const ConnectState > state_aux = new SimpleConnectState(u_tmp);
    Handle< const LinearOperator<LatticeFermion> > Maux = 
      S_w.hermitianLinOp(state_aux);
    
    
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
  EvenOddPrecHtContFrac5DFermActArray::createState(const multi1d<LatticeColorMatrix>& u_,
						   const OverlapStateInfo& state_info) const
  {
    // HACK UP A LINEAR OPERATOR TO CHECK EIGENVALUES/VECTORS WITH
    multi1d<LatticeColorMatrix> u_tmp = u_;
    getFermBC().modifyU(u_tmp);

    Handle< FermBC<LatticeFermion> > fbc4 = new PeriodicFermBC<LatticeFermion>();
    UnprecWilsonFermAct S_w(fbc4, Real(-params.OverMass));

    Handle< const ConnectState > state_aux = new SimpleConnectState(u_tmp);
    Handle< const LinearOperator<LatticeFermion> > Maux = 
      S_w.hermitianLinOp(state_aux);
    
    
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
