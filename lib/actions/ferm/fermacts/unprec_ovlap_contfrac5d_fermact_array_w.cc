// $Id: unprec_ovlap_contfrac5d_fermact_array_w.cc,v 1.1 2004-09-27 14:58:43 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned extended-Overlap (5D) (Naryanan&Neuberger) action
 */

#include "chromabase.h"

#include "actions/ferm/fermacts/unprec_wilson_fermact_w.h"
#include "actions/ferm/fermacts/ovlap_contfrac5d_fermact_array_w.h"
#include "actions/ferm/linop/ovlap_contfrac5d_linop_array_w.h"
#include "actions/ferm/linop/ovlap_contfrac5d_nonhermop_array_w.h"

#include "actions/ferm/linop/lmdagm.h"
#include "actions/ferm/invert/invcg2_array.h"
#include "zolotarev.h"

#include "actions/ferm/fermacts/fermfactory_w.h"
#include "io/enum_io/enum_io.h"
#include "io/overlap_state_info.h"

namespace Chroma
{

  //! Hooks to register the class with the fermact factory
  namespace UnprecOvlapContFrac5DFermActEnv
  {
    //! Callback function
    WilsonTypeFermAct<LatticeFermion>* createFermAct(Handle< FermBC< multi1d<LatticeFermion> > > fbc,
						     XMLReader& xml_in,
						     const std::string& path)
    {
      return new UnprecOvlapContFrac5DFermAct(fbc, UnprecOvlapContFrac5DFermActParams(xml_in, path));
    }

    //! Name to be used
    const std::string name = "UNPRECONDITIONED_OVERLAP_CONTINUED_FRACTION_5D";

    //! Register the Wilson fermact
    const bool registered = TheWilsonTypeFermActFactory::Instance().registerObject(name, createFermAct);
  }


  //! Read XML
  UnprecOvlapContFrac5DFermActParams::UnprecOvlapContFrac5DFermActParams(XMLReader& xml, const std::string& path)
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
      else {
	throw "No auxilliary action";
      }


      read(in, "Mass", Mass);
      read(in, "RatPolyDeg", RatPolyDeg);

      if( in.count("ApproximationType") == 1 ) { 
	read(in, "ApproximationType", approximation_type);
      }
      else { 
	// Default coeffs are Zolotarev
	approximation_type = COEFF_TYPE_ZOLOTAREV;
      }

    }
    catch( const string &e ) {
      QDPIO::cerr << "Caught Exception reading Zolo5D Fermact params: " << e << endl;
      QDP_abort(1);
    }
  }

  void read(XMLReader& xml_in, const string& path,
	    UnprecOvlapContFrac5DFermActParams& param) {

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
 
    pop(xml_out);


    if( path != "." ) { 
      pop(xml_out);
    }
  }

  
// Construct the action out of a parameter structure
  UnprecOvlapContFrac5DFermActArray::UnprecOvlapContFrac5DFermActArray(
		        Handle< FermBC< multi1d< LatticeFermion> > > fbc_a_, 
			const UnprecOvlapContFrac5DFermActParams& params_)
			
    fbc(fbc_a_), params(params_) {
  
    // Check RatPolyDeg is even
    if ( params.RatPolyDeg % 2 == 0 ) { 
      QDP_error_exit("For Now (and possibly forever), 5D Operators can only be constructed with ODD approximation order. You gave an even one: =%d\n", params.RatPolyDeg);
    }
    N5 = params.RatPolyDeg;
  
    // Construct the fermact 
      XMLReader  fermacttop(xml_s);
    const string fermact_path = "/AuxFermAct";

    // In case I fail to upcast to the UnprecWilsonType FermAct
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

      // I need a 4D FermAct for this beast so that I can create the
      // 4D Auxiliary action -- WHAT FOLLOWS IS A BEASTLY HACK:
      //
      // I am going to create a 4D trivial (all periodic) BC to wire
      // into the auxiliary operator. This is OK, because the 
      // createState() will use the 5D one and the 4D one is redundant
      // in this case. The auxiliary fermacts createState() with the 
      // trivial BC's will NEVER be called. It would be nice if I could
      // somehow convert the FermBC< multi1d<T> > to FermBC< T >
      // but that is something I need to discuss with Robert.
      Handle< FermBC<LatticeFermion> > fbc4( new PeriodicFermBC<LatticeFermion>() );

     
      // Generic Wilson-Type stuff
      WilsonTypeFermAct<LatticeFermion>* S_f =
	TheWilsonTypeFermActFactory::Instance().createObject(auxfermact,
							     fbc4,
							     fermacttop,
							     fermact_path);

      UnprecWilsonTypeFermAct<LatticeFermion>* S_aux_ptr; 
      S_aux_ptr = dynamic_cast<UnprecWilsonTypeFermAct<LatticeFermion>*>(S_f);

      // Dumbass User specifies something that is not UnpreWilsonTypeFermAct
      // dynamic_cast MUST be checked for 0
      if( S_aux_ptr == 0 ) throw UnprecCastFailure(auxfermact);
     

      // Drop AuxFermAct into a Handle immediately.
      // This should free things up at the end
      Handle<UnprecWilsonTypeFermAct<LatticeFermion> >  S_w(S_aux_ptr);
      S_aux = S_w;
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

    break;
    default:
      QDPIO::cerr << "Auxiliary Fermion Action Unsupported" << endl;
      QDP_abort(1);
    }
  }

  void
  UnprecOvlapContFrac5DFermActArray::init(Real& scale_fac,
				multi1d<Real>& alpha,
				multi1d<Real>& beta,
				int& NEig,
				multi1d<Real>& EigValFunc,
				const OverlapConnectState<LatticeFermion>& state) const
  {
  
    int NEigVal = state.getEigVal().size();
    if( NEigVal == 0 ) {
      NEig = 0;
    }
    else {
      NEig = NEigVal;
    }

    scale_fac = Real(1) / state.getApproxMax();
    Real eps = state.getApproxMin() * scale_fac;
    switch(params.approximation_type) { 
    case COEFF_TYPE_ZOLOTAREV:
      QDPIO::cout << "Initing Linop with Zolotarev Coefficients" << endl;
      QDPIO::cout << "  MaxCGInner =  " << params.invParamInner.MaxCG << endl;
      QDPIO::cout << "  RsdCGInner =  " << params.invParamInner.RsdCG << endl;
      QDPIO::cout << "  NEigVal    =  " << NEigVal << endl;
      
      type = 0;
      rdata = zolotarev(toFloat(eps), params.RatPolyDeg, type);    
      break;

    case COEFF_TYPE_TANH:
      QDPIO::cout << "Initing Linop with Higham Rep tanh Coefficients" << endl;
      QDPIO::cout << "  MaxCGInner =  " << params.invParamInner.MaxCG << endl;
      QDPIO::cout << "  RsdCGInner =  " << params.invParamInner.RsdCG << endl;
      QDPIO::cout << "  NEigVal    =  " << NEigVal << endl;
      
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

    int ZoloN5 = rdata -> db;
    if ( ZoloN5 != N5 ) { 
      QDP_error_exit("ZoloN5 and N5 are mismatched. ZoloN5=%d N5=%d\n", 
		     ZoloN5, N5);
    }

    // The coefficients from the continued fraction
    beta.resize(N5);
    for(int i = 0; i < N5; i++) { 
      beta[i] = rdata->beta[i];
    }


    alpha.resize(N5);
    for(int i = 0; i < N5; i++) {
      alpha[i] = Real(1);
    }

    alpha[N5-1] = rdata->beta[N5-1];

    // For the moment choose gamma = 1/sqrt(beta) */
    // except for gamma(N5-1) which always has to be set to 1 */
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
    
    
    QDPIO::cout << "UnprecOvlapContfrac5d: " 
		<< " N5=" << N5 << " scale=" << scale_fac
		<< " Nwils = " << NEigVal << " Mass=" << Mass << endl ;
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
    default:
      QDPIO::cerr << "Unknown coefficient type " << params.approximation_type 
		  << endl;
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
  const LinearOperator<multi1d<LatticeFermion> >* 
  UnprecOvlapContFrac5DFermActArray::linOp(Handle<const ConnectState> state_) const
  {
    START_CODE();
    const OverlapConnectState<LatticeFermion>& state = dynamic_cast<const OverlapConnectState<LatticeFermion>&>(*state_);

    if (state.getEigVec().size() != state.getEigVal().size())
      QDP_error_exit("UnprecOvlapContFrac5DFermActArray: inconsistent sizes of eigenvectors and values");

    int NEigVal = state.getEigVal().size();
    int NEig;

    multi1d<Real> alpha;
    multi1d<Real> beta;
    Real scale_factor;
    multi1d<Real> EigValFunc(NEigVal);

    init(scale_factor, alpha, beta, NEig, EigValFunc, state);

    return new UnprecOvlapContFrac5DLinOpArray( *S_aux,
				      state_,
				      Mass,
				      N5,
				      scale_factor,
				      alpha,
				      beta,
				      NEig,
				      EigValFunc,
				      state.getEigVec() );
				   
  

  }

//! Produce a linear operator for this action
/*!
 * The operator acts on the entire lattice
 *
 * \param state	    gauge field     	       (Read)
 */
  const LinearOperator<multi1d<LatticeFermion> >* 
  UnprecOvlapContFrac5DFermActArray::lnonHermLinOp(Handle<const ConnectState> state_) const
  {
    START_CODE();
    const OverlapConnectState<LatticeFermion>& state = dynamic_cast<const OverlapConnectState<LatticeFermion>&>(*state_);

    if (state.getEigVec().size() != state.getEigVal().size())
      QDP_error_exit("UnprecOvlapContFrac5DFermActArray: inconsistent sizes of eigenvectors and values");

    int NEigVal = state.getEigVal().size();
    int NEig;

    multi1d<Real> alpha;
    multi1d<Real> beta;
    Real scale_factor;
    multi1d<Real> EigValFunc(NEigVal);

    init(scale_factor, alpha, beta, NEig, EigValFunc, state);

    return new UnprecOvlapContFrac5DNonHermOpArray( *S_aux,
					  state_,
					  Mass,
					  N5,
					  scale_factor,
					  alpha,
					  beta,
					  NEig,
					  EigValFunc,
					  state.getEigVec() );
				   
  

  }


//! Produce a M^dag.M linear operator for this action
/*!
 * The operator acts on the entire lattice *
 * \param state	    gauge field     	       (Read)
 */
  const LinearOperator<multi1d<LatticeFermion> >* 
  UnprecOvlapContFrac5DFermActArray::lMdagM(Handle<const ConnectState> state) const
  {
    return new lmdagm<multi1d<LatticeFermion> >(linOp(state));
  }

//! Produce a M^dag.M linear operator for this action
/*!
 * The operator acts on the entire lattice *
 * \param state	    gauge field     	       (Read)
 */
  const LinearOperator<multi1d<LatticeFermion> >* 
  UnprecOvlapContFrac5DFermActArray::lnonHermMdagM(Handle<const ConnectState> state) const
  {
    return new lmdagm<multi1d<LatticeFermion> >(lnonHermLinOp(state));
  }


  const OverlapConnectState<LatticeFermion>*
  UnprecOvlapContFrac5DFermActArray::createState(const multi1d<LatticeColorMatrix>& u_) const
  {
    multi1d<LatticeColorMatrix> u_tmp = u_;
    getFermBC().modifyU(u_tmp);

    Real approxMin = 0.0;
    Real approxMax = Real(2)*Real(Nd);
    return new OverlapConnectState<LatticeFermion>(u_tmp, approxMin, approxMax);
  }


  const OverlapConnectState<LatticeFermion>*
  UnprecOvlapContFrac5DFermActArray::createState(const multi1d<LatticeColorMatrix>& u_,
				       const Real& approxMin_) const 
  {
    if ( toBool( approxMin_ < Real(0) )) { 
      ostringstream error_str;
      error_str << "UnprecOvlapContFrac5DFermActArray: approxMin_ has to be positive" << endl;
      throw error_str.str();
    }
 

    // First put in the BC
    multi1d<LatticeColorMatrix> u_tmp = u_;
    getFermBC().modifyU(u_tmp);

    Real approxMax = Real(2)*Real(Nd);
    return new OverlapConnectState<LatticeFermion>(u_tmp, approxMin_, approxMax);
  }


  const OverlapConnectState<LatticeFermion>*
  UnprecOvlapContFrac5DFermActArray::createState(const multi1d<LatticeColorMatrix>& u_,
				       const Real& approxMin_,
				       const Real& approxMax_) const
  {
    ostringstream error_str;
  
 
    if ( toBool(approxMin_ < 0 )) { 
      error_str << "UnprecOvlapContFrac5DFermActArray: approxMin_ has to be positive" << endl;
      throw error_str.str();
    }

    if ( toBool(approxMax_ < approxMin_) ) { 
      error_str << "UnprecOvlapContFrac5DFermActArray: approxMax_ has to be larger than approxMin_" << endl;
      throw error_str.str();
    }
 

    // First put in the BC
    multi1d<LatticeColorMatrix> u_tmp = u_;
    getFermBC().modifyU(u_tmp);

    return new OverlapConnectState<LatticeFermion>(u_tmp, approxMin_, approxMax_);
  }



  const OverlapConnectState<LatticeFermion>*
  UnprecOvlapContFrac5DFermActArray::createState(const multi1d<LatticeColorMatrix>& u_,
				       const multi1d<Real>& lambda_lo_, 
				       const multi1d<LatticeFermion>& evecs_lo_,
				       const Real& lambda_hi_) const
  {
    ostringstream error_str;

    if ( lambda_lo_.size() == 0 ) {
      error_str << "Attempt to createState with 0 e-values and no approxMin" << endl;
      throw error_str.str();
    }

    if ( lambda_lo_.size() != evecs_lo_.size() ) {
      error_str << "Attempt to createState with no of low eigenvalues != no of low eigenvectors" << endl;
      throw error_str.str();
    }

    Real approxMax = lambda_hi_;
    Real approxMin = fabs(lambda_lo_[ lambda_lo_.size() - 1 ]);

    // First put in the BC
    multi1d<LatticeColorMatrix> u_tmp = u_;
    getFermBC().modifyU(u_tmp);

    return new OverlapConnectState<LatticeFermion>(u_tmp, lambda_lo_, evecs_lo_, lambda_hi_, approxMin, approxMax);
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
  UnprecOvlapContFrac5DFermActArray::qprop(LatticeFermion& psi, 
				 Handle<const ConnectState> state, 
				 const LatticeFermion& chi, 
				 enum InvType invType,
				 const Real& RsdCG, 
				 int MaxCG, int& ncg_had) const
  {

    START_CODE();

    const int  N5 = size();   // array size better match
    const Real Mass = quark_mass();
    int n_count;
  
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

    // Construct the linear operator
    Handle<const LinearOperator< multi1d<LatticeFermion> > > A(linOp(state));
    LinearOperator< multi1d<LatticeFermion> >* AdagA;

    switch(invType) {
    case CG_INVERTER: 

      // Zero out 5D vectors
      for(int i=0; i < N5; i++) {
	psi5[i] = zero;
	chi5[i] = zero;
      }


   
      // Use psi5 as temporary
      psi5[N5-1] = Gamma(G5) * chi;

      // chi5 now holds  D^{dag} gamma_5 chi
      (*A)(chi5, psi5, MINUS);

      // Now use psi5 as it was meant to be
      psi5[N5-1] = psi;

      // psi5 = (M^{dag} M)^(-1) M^{dag} * gamma_5 * chi5
      // psi5[N5]  = (1 - m)/2 D^{-1}(m) chi [N5]
      InvCG2(*A, chi5, psi5, RsdCG, MaxCG, n_count);


   

      break;
  
    case MR_INVERTER:
      QDP_error_exit("Unsupported inverter type", invType);
      break;

    case BICG_INVERTER:
      QDP_error_exit("Unsupported inverter type", invType);
      break;
  
    default:
      QDP_error_exit("Unknown inverter type", invType);
    }
  
    if ( n_count == MaxCG )
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


  const OverlapConnectState<LatticeFermion>*
  UnprecOvlapContFrac5DFermActArray::createState(const multi1d<LatticeColorMatrix>& u_,
				       const OverlapStateInfo& state_info,
				       XMLWriter& xml_out,
				       Real wilsonMass) const
  {
    push(xml_out, "Zolo5DCreateState");


    // If No eigen values specified use min and max
    if ( state_info.getNWilsVec() == 0 ) 
    { 
      write(xml_out, "ApproxMin", state_info.getApproxMin());
      write(xml_out, "ApproxMax", state_info.getApproxMax());
      pop(xml_out);

      return createState(u_,
			 state_info.getApproxMin(),
			 state_info.getApproxMax());
    }
    else
    {
      QDPIO::cout << "Warning!!!! Using 5D op with projected e-values is dubious "
		  << " at this time " << endl;

      // If there are eigen values, either load them, 
      if( state_info.loadEigVec() ) 
      {
	ChromaWilsonRitz_t ritz_header;
	multi1d<Real> lambda_lo;
	multi1d<LatticeFermion> eigv_lo;
	Real lambda_hi;
	const EigenIO_t& eigen_io = state_info.getEigenIO();

	push(xml_out, "EigenSystem");
	if( eigen_io.eigen_filefmt == EVEC_TYPE_SCIDAC ) { 
	  readEigen(ritz_header, lambda_lo, eigv_lo, lambda_hi, 
		    eigen_io.eigen_file,
		    state_info.getNWilsVec(),
		    QDPIO_SERIAL);
	  write(xml_out, "OriginalRitzHeader", ritz_header);
	}
	else if ( eigen_io.eigen_filefmt == EVEC_TYPE_SZIN ) { 

	  if( toBool( fabs(wilsonMass) > 8 ) ){
	    QDPIO::cerr << "OverMass unspecified, or | OverMass | > 8" << endl;
	    QDPIO::cerr << "The wilson mass is needed to rescale the eigenvalues" << endl;
	    QDP_abort(1);
	  }

	  readEigenSzin(lambda_lo, eigv_lo, lambda_hi, state_info.getNWilsVec(), eigen_io.eigen_file);
	
	  // Now I need to scale things by the wilson mass (Nd + m)
	  for(int i=0; i < lambda_lo.size(); i++) { 
	    lambda_lo[i] *= (Real(Nd) + wilsonMass);
	  }

	  lambda_hi *= (Real(Nd) + wilsonMass);

	}
	else {
	  QDPIO::cerr << "Unsupported Eigenvector format for reading " << endl;
	  QDP_abort(1);
	}

	write(xml_out, "lambda_lo", lambda_lo);
	write(xml_out, "lambda_high", lambda_hi);
         
	Handle< const ConnectState > wils_connect_state = S_aux->createState(u_);
	Handle< const LinearOperator<LatticeFermion> > H = S_aux->gamma5HermLinOp(wils_connect_state);

      	      
	multi1d<Double> check_norm(state_info.getNWilsVec());
	multi1d<Double> check_norm_rel(state_info.getNWilsVec());
	for(int i=0; i < state_info.getNWilsVec() ; i++) { 
	  LatticeFermion Me;
	  (*H)(Me, eigv_lo[i], PLUS);
	
	  LatticeFermion lambda_e;
	
	  lambda_e = lambda_lo[i]*eigv_lo[i];
	  LatticeFermion r_norm = Me - lambda_e;
	  check_norm[i] = sqrt(norm2(r_norm));
	  check_norm_rel[i] = check_norm[i]/fabs(Double(lambda_lo[i]));
	
	  QDPIO::cout << "Eigenpair " << i << " Resid Norm = " 
		      << check_norm[i] << " Resid Rel Norm = " << check_norm_rel[i] << endl;
	}
	write(xml_out, "eigen_norm", check_norm);
	write(xml_out, "eigen_rel_norm", check_norm_rel);
      
	pop(xml_out); // Eigensystem

	pop(xml_out); // Zolo5DCreateState
	return createState(u_, lambda_lo, eigv_lo, lambda_hi);
      }
      else if( state_info.computeEigVec() ) {
	QDPIO::cerr << "Recomputation of eigensystem not yet implemented" << endl;
	QDP_abort(1);
      }
      else {
	QDPIO::cerr << "I have to create a state without min/max, loading/computing eigenvectors/values. How do I do that? "<< endl;
	QDP_abort(1);
      }
    }

    return 0;  // Should never get here. Make compilers happy.
  }

}
