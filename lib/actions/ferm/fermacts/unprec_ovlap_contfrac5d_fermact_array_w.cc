// $Id: unprec_ovlap_contfrac5d_fermact_array_w.cc,v 3.6 2007-02-22 21:11:45 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned extended-Overlap (5D) (Naryanan&Neuberger) action
 */

#include "chromabase.h"

#include "actions/ferm/fermacts/unprec_wilson_fermact_w.h"
#include "actions/ferm/fermacts/unprec_ovlap_contfrac5d_fermact_array_w.h"
#include "actions/ferm/linop/unprec_ovlap_contfrac5d_linop_array_w.h"
#include "actions/ferm/linop/unprec_ovlap_contfrac5d_nonhermop_array_w.h"
#include "actions/ferm/linop/unprec_ovlap_contfrac5d_pv_linop_array_w.h"

#include "actions/ferm/invert/invcg2_array.h"
#include "zolotarev.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermstates/simple_fermstate.h"
//#include "actions/ferm/fermstates/ferm_createstate_reader_w.h"
#include "actions/ferm/fermbcs/fermbcs_reader_w.h"

#include "io/enum_io/enum_io.h"
#include "io/overlap_state_info.h"

#include "actions/ferm/invert/syssolver_cg_params.h"

namespace Chroma
{

  //! Hooks to register the class with the fermact factory
  namespace UnprecOvlapContFrac5DFermActArrayEnv
  {
    //! Callback function
    WilsonTypeFermAct5D<LatticeFermion,
			multi1d<LatticeColorMatrix>,
			multi1d<LatticeColorMatrix> >* createFermAct5D(XMLReader& xml_in,
								       const std::string& path)
    {
      return new UnprecOvlapContFrac5DFermActArray(WilsonTypeFermBCEnv::reader(xml_in, path), 
						   UnprecOvlapContFrac5DFermActParams(xml_in, path));
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
    const std::string name = "UNPRECONDITIONED_OVERLAP_CONTINUED_FRACTION_5D";
    
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
  } // End Namespace UnprecOvlapContFrac5DFermActArrayEnv

  
  //! Read XML
  UnprecOvlapContFrac5DFermActParams::UnprecOvlapContFrac5DFermActParams(XMLReader& xml, const std::string& path)
  {
    XMLReader in(xml, path);
    
    try { 
      if(in.count("AuxFermAct") == 1 ) {  
	XMLReader xml_tmp(in, "AuxFermAct");
	std::ostringstream os;
	xml_tmp.print(os);
	AuxFermAct = os.str();
      }
      else {
	throw std::string("No auxilliary action");
      }
      
      read(in, "Mass", Mass);
      read(in, "RatPolyDeg", RatPolyDeg);
      
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
    }
    catch( const string &e ) {
      QDPIO::cerr << "Caught Exception reading unprec ContFrac Fermact params: " << e << endl;
      QDP_abort(1);
    }
  }

  void read(XMLReader& xml_in, const string& path,
	    UnprecOvlapContFrac5DFermActParams& param) 
  {
    
    UnprecOvlapContFrac5DFermActParams tmp(xml_in, path);
    param = tmp;
  }
  
  void write(XMLWriter& xml_out, const string& path, const UnprecOvlapContFrac5DFermActParams& p)
  {
    if ( path != "." ) { 
      push( xml_out, path);
    }
    
    xml_out << p.AuxFermAct;
    write(xml_out, "Mass", p.Mass);
    write(xml_out, "RatPolyDeg", p.RatPolyDeg);
    write(xml_out, "ApproximationType", p.approximation_type);
    write(xml_out, "ApproxMin", p.ApproxMin);
    write(xml_out, "ApproxMax", p.ApproxMax);
    
    pop(xml_out);
    
    
    if( path != "." ) { 
      pop(xml_out);
    }
  }

  
  
  // Construct the action out of a parameter structure
  UnprecOvlapContFrac5DFermActArray::UnprecOvlapContFrac5DFermActArray(
    Handle< FermBC<T,P,Q> > fbc_, 
    const UnprecOvlapContFrac5DFermActParams& params_) :
    fbc(fbc_), params(params_) 
  {
    QDPIO::cout << UnprecOvlapContFrac5DFermActArrayEnv::name << ": entering constructor" << endl;

    // Sanity check
    if (fbc.operator->() == 0)
    {
      QDPIO::cerr << UnprecOvlapContFrac5DFermActArrayEnv::name << ": error: fbc is null" << endl;
      QDP_abort(1);
    }

    // Fake a creator. This should be cleaned up
    Handle< CreateFermState<T,P,Q> > cfs_(new CreateSimpleFermState<T,P,Q>(fbc));
    cfs = cfs_;

    
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

    // Construct the fermact 
    std::istringstream  xml_s(params.AuxFermAct);
    XMLReader  fermacttop(xml_s);
    const string fermact_path = "/AuxFermAct";
    
    // In case I fail to upcast to the UnprecWilsonType FermAct
    struct UnprecCastFailure {
      UnprecCastFailure(std::string e) : auxfermact(e) {};
      const string auxfermact;
    };
    
    try {
      string auxfermact;
      read(fermacttop, fermact_path + "/FermAct", auxfermact);
      QDPIO::cout << "AuxFermAct: " << auxfermact << endl;
            
      // Generic Wilson-Type stuff
      FermionAction<T,P,Q>* S_f =
	TheFermionActionFactory::Instance().createObject(auxfermact,
							 fermacttop,
							 fermact_path);
      
      UnprecWilsonTypeFermAct<T,P,Q>* S_aux_ptr = 
	dynamic_cast<UnprecWilsonTypeFermAct<T,P,Q>*>(S_f);
      
      // Dumbass User specifies something that is not UnpreWilsonTypeFermAct
      // dynamic_cast MUST be checked for 0
      if( S_aux_ptr == 0 ) throw UnprecCastFailure(auxfermact);
      
      // Drop AuxFermAct into a Handle immediately.
      // This should free things up at the end
      Handle<UnprecWilsonTypeFermAct<T,P,Q> >  S_w(S_aux_ptr);
      S_aux = S_w;
    }
    catch( const UnprecCastFailure& e) 
    {
      // Breakage Scenario
      QDPIO::cerr << "Unable to upcast auxiliary fermion action to "
		  << "UnprecWilsonTypeFermAct " << endl;
      QDPIO::cerr << UnprecOvlapContFrac5DFermActArrayEnv::name << " does not support even-odd preconditioned "
		  << "auxiliary FermActs" << endl;
      QDPIO::cerr << "You passed : " << endl;
      QDPIO::cerr << e.auxfermact << endl;
      QDP_abort(1);
    }
    catch (const std::exception& e) 
    {
      // General breakage Scenario
      QDPIO::cerr << "Error reading data: " << e.what() << endl;
      throw;
    }
    catch( const std::string& e ) 
    { 
      QDPIO::cerr << "Caught Exception" << e << endl;
      QDP_abort(1);
    }
  }


  void
  UnprecOvlapContFrac5DFermActArray::init(Real& scale_fac,
					  multi1d<Real>& alpha,
					  multi1d<Real>& beta,
					  int& NEig,
					  multi1d<Real>& EigValFunc,
					  const OverlapConnectState& state) const
  {
    int NEigVal = state.getEigVal().size();
    if( NEigVal == 0 ) {
      NEig = 0;
    }
    else {
      NEig = NEigVal;
    }
    
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

    
    QDPIO::cout << "UnprecOvlapContfrac5d: " << endl
                << "  Degree=" << params.RatPolyDeg 
		<< "  N5=" << N5 << " scale=" << scale_fac
		<< "  Nwils = " << NEigVal << " Mass=" << params.Mass << endl ;
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

    // We will also compute te 'function of the eigenvalues
    if( NEig > 0 ) {
      for(int i=0; i <  NEigVal; i++) {
	if( toBool( state.getEigVal()[i] > Real(0) ) ) {
	  EigValFunc[i] = Real(1);
	}
	else if( toBool( state.getEigVal()[i] < Real(0) ) ) {
	  EigValFunc[i] = Real(-1);
	}
	else {
	  EigValFunc[i] = 0;
	}
      }
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
  UnprecLinearOperatorArray<LatticeFermion,
			    multi1d<LatticeColorMatrix>,
			    multi1d<LatticeColorMatrix> >* 
  UnprecOvlapContFrac5DFermActArray::linOp(Handle< FermState<T,P,Q> > state_) const
  {
    START_CODE();
    try { 
      
      // This throws a bad cast exception if the cast fails
      // Hence the "try" above
      const OverlapConnectState& state = 
	dynamic_cast<OverlapConnectState&>(*state_);
      
      if (state.getEigVec().size() != state.getEigVal().size()) {
	QDPIO::cerr << "UnprecOvlapContFrac5DFermActArray::linOp(): ";
	QDPIO::cerr << "UnprecOvlapContFrac5DFermActArray: inconsistent sizes of eigenvectors and values" 
		    << "state.getEigVec.size() = " << state.getEigVec().size() 
		    << " state.getEigVal.size() = " << state.getEigVal().size()
		    << endl;
	QDP_abort(1);
      }
      
      int NEigVal = state.getEigVal().size();
      int NEig;
      
      multi1d<Real> alpha;
      multi1d<Real> beta;
      Real scale_factor;
      multi1d<Real> EigValFunc(NEigVal);
      
      init(scale_factor, alpha, beta, NEig, EigValFunc, state);
      
      return new UnprecOvlapContFrac5DLinOpArray(*S_aux,
						 state_,
						 params.Mass,
						 N5,
						 scale_factor,
						 alpha,
						 beta,
						 NEig,
						 EigValFunc,
						 state.getEigVec(),
						 isLastZeroP);
      
    }
    catch( bad_cast ) 
    {
      QDPIO::cerr << "UnprecOvlapContFrac5DFermActArray::linOp(): ";
      QDPIO::cerr << "Couldnt cast FermState<T,P,Q>  to OverlapFermState<T,P,Q>  " 
		  << endl;
      QDP_abort(1);
    }
    
    // Should never get here... Just to satisfy type system
    return 0;
  }
  

  //! Produce a linear operator for this action
  /*!
   * The operator acts on the entire lattice
   *
   * \param state	    gauge field     	       (Read)
   */
  UnprecLinearOperatorArray<LatticeFermion,
			    multi1d<LatticeColorMatrix>,
			    multi1d<LatticeColorMatrix> >* 
  UnprecOvlapContFrac5DFermActArray::linOpPV(Handle< FermState<T,P,Q> > state_) const
  {
    START_CODE();
    try { 
      
      // This throws a bad cast exception if the cast fails
      // Hence the "try" above
      const OverlapConnectState& state = 
	dynamic_cast<OverlapConnectState&>(*state_);
      
      if (state.getEigVec().size() != state.getEigVal().size()) {
	QDPIO::cerr << "UnprecOvlapContFrac5DFermActArray::linOp(): ";
	QDPIO::cerr << "UnprecOvlapContFrac5DFermActArray: inconsistent sizes of eigenvectors and values" 
		    << "state.getEigVec.size() = " << state.getEigVec().size() 
		    << " state.getEigVal.size() = " << state.getEigVal().size()
		    << endl;
	QDP_abort(1);
      }
      
      int NEigVal = state.getEigVal().size();
      int NEig;
      
      multi1d<Real> alpha;
      multi1d<Real> beta;
      Real scale_factor;
      multi1d<Real> EigValFunc(NEigVal);
      
      init(scale_factor, alpha, beta, NEig, EigValFunc, state);
      
      return new UnprecOvlapContFrac5DPVLinOpArray( *S_aux,
						    state_,
						    params.Mass,
						    N5,
						    scale_factor,
						    alpha,
						    beta,
						    NEig,
						    EigValFunc,
						    state.getEigVec(),
						    isLastZeroP);
      
    }
    catch( bad_cast ) { 
      QDPIO::cerr << "UnprecOvlapContFrac5DFermActArray::linOp(): ";
      QDPIO::cerr << "Couldnt cast FermState<T,P,Q>  to OverlapFermState<T,P,Q>  " 
		  << endl;
      QDP_abort(1);
    }
    
    // Should never get here... Just to satisfy type system
    return 0;
  }
  

  //! Produce the non-hermitian version of the operator
  /*!
   * The operator acts on the entire lattice
   *
   * \param state	    gauge field     	       (Read)
   */
  LinearOperatorArray<LatticeFermion>* 
  UnprecOvlapContFrac5DFermActArray::lnonHermLinOp(Handle< FermState<T,P,Q> > state_) const
  {
    START_CODE();
    
    try { 
      // This cast throws an exception if it fails hence the " try " above
      const OverlapConnectState& state = 
	dynamic_cast<OverlapConnectState&>(*state_);
      
      if (state.getEigVec().size() != state.getEigVal().size()) {
	QDPIO::cerr << "UnprecOvlapContFrac5DFermActArray::lnonHermLinOp(): ";
	QDPIO::cerr << "UnprecOvlapContFrac5DFermActArray: inconsistent sizes of eigenvectors and values" 
		    << "state.getEigVec.size() = " << state.getEigVec().size() 
		    << " state.getEigVal.size() = " << state.getEigVal().size()
		    << endl;
	QDP_abort(1);
      }
      
      
      int NEigVal = state.getEigVal().size();
      int NEig;
      
      multi1d<Real> alpha;
      multi1d<Real> beta;
      Real scale_factor;
      multi1d<Real> EigValFunc(NEigVal);
      
      init(scale_factor, alpha, beta, NEig, EigValFunc, state);
      
      return new UnprecOvlapContFrac5DNonHermOpArray( *S_aux,
						      state_,
						      params.Mass,
						      N5,
						      scale_factor,
						      alpha,
						      beta,
						      NEig,
						      EigValFunc,
						      state.getEigVec() );
    }
    catch( bad_cast ) 
    {
      QDPIO::cerr << "UnprecOvlapContFrac5DFermActArray::lnonHermLinOp(): ";
      QDPIO::cerr << "Couldnt cast FermState<T,P,Q>  to OverlapConnectState  " 
		  << endl;
      QDP_abort(1);
    }

    return 0;
  }


  //! Produce a M^dag.M linear operator for the non hermitian operator
  /*!
   * The operator acts on the entire lattice *
   * \param state	    gauge field     	       (Read)
   */
  LinearOperatorArray<LatticeFermion>* 
  UnprecOvlapContFrac5DFermActArray::lnonHermMdagM(Handle< FermState<T,P,Q> > state) const
  {
    return new MdagMLinOpArray<T>(lnonHermLinOp(state));
  }
  


  //! Propagator of an un-preconditioned Extended-Overlap linear operator
  /*! \ingroup qprop
   *
   * Propagator solver for Extended overlap fermions
   */
  template<typename T>
  class OvUnprecCF5DQprop : public SystemSolver<T>
  {
  public:
    //! Constructor
    /*!
     * \param A_        Linear operator ( Read )
     * \param Mass_      quark mass ( Read )
     */
    OvUnprecCF5DQprop(Handle< LinearOperatorArray<T> > A_,
		      const Real& Mass_,
		      const SysSolverCGParams& invParam_) : 
      A(A_), Mass(Mass_), invParam(invParam_) {}

    //! Destructor is automatic
    ~OvUnprecCF5DQprop() {}

    //! Return the subset on which the operator acts
    const Subset& subset() const {return all;}

    //! Solver the linear system
    /*!
     * \param psi      quark propagator ( Modify )
     * \param chi      source ( Read )
     * \return number of CG iterations
     */
    SystemSolverResults_t operator() (T& psi, const T& chi) const
    {
      START_CODE();
    
      SystemSolverResults_t res;
      const int N5 = A->size();
    
      int G5 = Ns*Ns - 1;
    
      // Initialize the 5D fields
      multi1d<LatticeFermion> chi5(N5);
      multi1d<LatticeFermion> psi5(N5);

      
      // For reasons I do not appreciate doing the solve as
      //  M^{dag} M psi = M^{dag} chi
      //  seems a few iterations faster and more accurate than
      //  Doing M^{dag} M psi = chi
      //  and then applying M^{dag}

      // So first get  M^{dag} gamma_5 chi into the source.
//      if( invType == "CG_INVERTER") 
      {
	multi1d<LatticeFermion> tmp5(N5);
	    
	// Zero out 5D vectors
	for(int i=0; i < N5; i++) {
	  psi5[i] = zero;
	  chi5[i] = zero;
	  tmp5[i] = zero;
	}
	    
	// Set initial guess
	psi5[N5-1] = psi;
	    
	// set up gamma_5 chi into chi-1
	chi5[N5-1] = Gamma(G5)*chi;
	(*A)(tmp5, chi5, MINUS);
	
	// psi5 = (M^{dag} M)^(-1) M^{dag} * gamma_5 * chi5
	// psi5[N5]  = (1 - m)/2 D^{-1}(m) chi [N5]
	res = InvCG2(*A, tmp5, psi5, invParam.RsdCG, invParam.MaxCG);
      }
//      else
//      {
//	QDPIO::cerr << UnprecOvlapContFrac5DFermActArrayEnv::name 
//		    << "Unsupported inverter type =" << invType << endl;
//	QDP_abort(1);
//      }
  
      if ( res.n_count == invParam.MaxCG )
	QDP_error_exit("no convergence in the inverter", res.n_count);
    
      // Compute residual
      {
	multi1d<T>  r(N5);
	(*A)(r, psi5, PLUS);
	r -= chi5;
	res.resid = sqrt(norm2(r));
      }

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

      return res;
    }

  private:
    // Hide default constructor
    OvUnprecCF5DQprop() {}

    Handle< LinearOperatorArray<T> > A;
    Real Mass;
    SysSolverCGParams invParam;
  };

  
  //! Propagator of an un-preconditioned Extended-Overlap linear operator
  SystemSolver<LatticeFermion>* 
  UnprecOvlapContFrac5DFermActArray::qprop(Handle< FermState<T,P,Q> > state,
					   const GroupXML_t& invParam) const
  {
    std::istringstream  is(invParam.xml);
    XMLReader  paramtop(is);

    return new OvUnprecCF5DQprop<T>(Handle< LinearOperatorArray<T> >(linOp(state)),
				    getQuarkMass(),
				    SysSolverCGParams(paramtop,invParam.path));
  }


  
  //! Create a ConnectState with just the gauge fields
  OverlapConnectState*
  UnprecOvlapContFrac5DFermActArray::createState(const multi1d<LatticeColorMatrix>& u_) const
  {
    try { 
      return new OverlapConnectState(fbc, u_);
    } 
    catch(const string& e) { 
      QDPIO::cerr << "Caught Exception: " << e << endl;
      QDP_abort(1);
    }
    
    return 0;
  }
  
  //! Create a ConnectState with just the gauge fields, and a lower
  //  approximation bound
  OverlapConnectState*
  UnprecOvlapContFrac5DFermActArray::createState(const multi1d<LatticeColorMatrix>& u_,
						 const Real& approxMin_) const 
  {
    try { 
      return new OverlapConnectState(fbc, u_, approxMin_);
    } 
    catch(const string& e) { 
      QDPIO::cerr << "Caught Exception: " << e << endl;
      QDP_abort(1);
    }

    return 0;
  }

  //! Create a connect State with just approximation range bounds
  OverlapConnectState*
  UnprecOvlapContFrac5DFermActArray::createState(const multi1d<LatticeColorMatrix>& u_,
						 const Real& approxMin_,
						 const Real& approxMax_) const
  {
    try { 
      return new OverlapConnectState(fbc, u_, approxMin_, approxMax_);
    } 
    catch(const string& e) { 
      QDPIO::cerr << "Caught Exception: " << e << endl;
      QDP_abort(1);
    }

    return 0;
  }
  
  //! Create OverlapConnectState with eigenvalues/vectors
  OverlapConnectState*
  UnprecOvlapContFrac5DFermActArray::createState(const multi1d<LatticeColorMatrix>& u_,
						 const multi1d<Real>& lambda_lo_, 
						 const multi1d<LatticeFermion>& evecs_lo_,
						 const Real& lambda_hi_) const
  {
    try { 
      return new OverlapConnectState(fbc, u_, lambda_lo_, evecs_lo_, lambda_hi_);
    } 
    catch(const string& e) { 
      QDPIO::cerr << "Caught Exception: " << e << endl;
      QDP_abort(1);
    }
    
    return 0;
  }    
  
  //! Create OverlapConnectState from XML
  OverlapConnectState*
  UnprecOvlapContFrac5DFermActArray::createState(const multi1d<LatticeColorMatrix>& u_,
						 XMLReader& state_info_xml,
						 const string& state_info_path) const
  {
    // HACK UP A LINEAR OPERATOR TO CHECK EIGENVALUES/VECTORS WITH
    Handle< FermState<T,P,Q>  > state_aux = new SimpleFermState<T,P,Q> (fbc, u_);
    Handle< LinearOperator<T> > Maux = S_aux->hermitianLinOp(state_aux);
    
    try {
      return new OverlapConnectState(fbc, u_, state_info_xml, state_info_path, *Maux);
    }
    catch(const string& e) { 
      QDPIO::cerr << "Caught Exception: " << e << endl;
      QDP_abort(1);
    }
    
    return 0;
  }

  //! Create OverlapConnectState from XML
  OverlapConnectState*
  UnprecOvlapContFrac5DFermActArray::createState(const multi1d<LatticeColorMatrix>& u_,
						 const OverlapStateInfo& state_info) const
  {
    // HACK UP A LINEAR OPERATOR TO CHECK EIGENVALUES/VECTORS WITH
    Handle< FermState<T,P,Q>  > state_aux = new SimpleFermState<T,P,Q> (fbc, u_);
    Handle< LinearOperator<T> > Maux = S_aux->hermitianLinOp(state_aux);
    
    try {
      return new OverlapConnectState(fbc, u_, state_info, *Maux);
    }
    catch(const string& e) { 
      QDPIO::cerr << "Caught Exception: " << e << endl;
      QDP_abort(1);
    }
    
    return 0;
  }

} // End namespace Chroma
