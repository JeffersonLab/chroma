// $Id: unprec_ht_contfrac5d_fermact_array_w.cc,v 1.1 2005-01-05 05:39:16 edwards Exp $
/*! \file
 *  \brief Unpreconditioned H_T kernel continued fraction (5D) action
 */

#include "chromabase.h"

#include "actions/ferm/fermacts/unprec_wilson_fermact_w.h"
#include "actions/ferm/fermacts/unprec_ht_contfrac5d_fermact_array_w.h"
#include "actions/ferm/linop/unprec_ht_contfrac5d_linop_array_w.h"

#include "actions/ferm/linop/lmdagm.h"
#include "actions/ferm/invert/invcg2_array.h"
#include "zolotarev.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermbcs/fermbcs_w.h"

#include "io/enum_io/enum_io.h"

namespace Chroma
{

  //! Hooks to register the class with the fermact factory
  namespace UnprecHTContFrac5DFermActArrayEnv
  {
    //! Callback function
    WilsonTypeFermAct5D<LatticeFermion,multi1d<LatticeColorMatrix> >* createFermAct5D(XMLReader& xml_in,
										      const std::string& path)
    {
      return new UnprecHTContFrac5DFermActArray(WilsonTypeFermBCArrayEnv::reader(xml_in, path), 
						UnprecHTContFrac5DFermActParams(xml_in, path));
    }

    //! Callback function
    /*! Differs in return type */
    FermionAction<LatticeFermion>* createFermAct(XMLReader& xml_in,
						 const std::string& path)
    {
      return createFermAct5D(xml_in, path);
    }

    //! Name to be used
    const std::string name = "UNPRECONDITIONED_HT_CONTINUED_FRACTION_5D";
    
    //! Register all the factories
    bool registerAll()
    {
      return Chroma::TheFermionActionFactory::Instance().registerObject(name, createFermAct)
	& Chroma::TheWilsonTypeFermAct5DFactory::Instance().registerObject(name, createFermAct5D);
    }

    //! Register the fermact
    const bool registered = registerAll();
  } // End Namespace UnprecHTContFrac5DFermActArrayEnv

  
  //! Read XML
  UnprecHTContFrac5DFermActParams::UnprecHTContFrac5DFermActParams(XMLReader& xml, const std::string& path)
  {
    XMLReader in(xml, path);
    
    try 
    {
      read(in, "OverMass", OverMass);
      read(in, "Mass", Mass);
      read(in, "b5", b5);
      read(in, "c5", c5);
      read(in, "RatPolyDeg", RatPolyDeg);

      read(in, "ApproxMin", ApproxMin);
      read(in, "ApproxMax", ApproxMax);
      
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
	    UnprecHTContFrac5DFermActParams& param) 
  {
    UnprecHTContFrac5DFermActParams tmp(xml_in, path);
    param = tmp;
  }

  
  void write(XMLWriter& xml_out, const string& path, const UnprecHTContFrac5DFermActParams& p)
  {
    if ( path != "." ) { 
      push( xml_out, path);
    }
    
    xml_out << p.AuxFermAct;
    write(xml_out, "OverMass", p.OverMass);
    write(xml_out, "Mass", p.Mass);
    write(xml_out, "RatPolyDeg", p.RatPolyDeg);
    write(xml_out, "ApproximationType", p.approximation_type);
    write(xml_out, "b5", p.b5);
    write(xml_out, "c5", p.c5);
    
    pop(xml_out);
    
    
    if( path != "." ) { 
      pop(xml_out);
    }
  }

  
  
  // Construct the action out of a parameter structure
  UnprecHTContFrac5DFermActArray::UnprecHTContFrac5DFermActArray(
    Handle< FermBC< multi1d< LatticeFermion> > > fbc_a_, 
    const UnprecHTContFrac5DFermActParams& params_) :
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
  UnprecHTContFrac5DFermActArray::init(Real& scale_fac,
				       Real& a5,
				       multi1d<Real>& alpha,
				       multi1d<Real>& beta) const
  {
    int type = 0;
    zolotarev_data *rdata;
    Real eps;

    switch(params.approximation_type) { 
    case COEFF_TYPE_ZOLOTAREV:
      scale_fac = Real(1) / params.ApproxMax;
      eps = params.ApproxMin * scale_fac;

      QDPIO::cout << "Initing Linop with Zolotarev Coefficients" << endl;
      rdata = zolotarev(toFloat(eps), params.RatPolyDeg, type);    
      break;

    case COEFF_TYPE_TANH:
      scale_fac = Real(1) / params.ApproxMax;
      eps = params.ApproxMin * scale_fac;

      QDPIO::cout << "Initing Linop with Higham Rep tanh Coefficients" << endl;
      rdata = higham(toFloat(eps), params.RatPolyDeg);

      break;
    case COEFF_TYPE_TANH_UNSCALED:
      scale_fac = Real(1);
      eps = params.ApproxMin;

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

    
    QDPIO::cout << "UnprecHTContfrac5d: " << endl
                << "Degree=" << params.RatPolyDeg 
		<< "N5=" << N5 << " scale=" << scale_fac
		<< " Mass=" << params.Mass << endl ;
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
	QDPIO::cout << "Approximation type " << type << " with R(0) =  infinity"
                    << endl;
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
  const UnprecLinearOperator< multi1d<LatticeFermion>, multi1d<LatticeColorMatrix> >* 
  UnprecHTContFrac5DFermActArray::linOp(Handle<const ConnectState> state) const
  {
    START_CODE();

    multi1d<Real> alpha;
    multi1d<Real> beta;
    Real scale_factor;
    Real a5;
      
    init(scale_factor, a5, alpha, beta);
      
    return new UnprecHTContFrac5DLinOpArray( state->getLinks(),
					     params.OverMass,
					     params.Mass,
					     N5,
					     b5,
					     c5,
					     scale_factor,
					     alpha,
					     beta,
					     isLastZeroP);
  }
  



  //! Produce a M^dag.M linear operator for this action
  /*!
   * The operator acts on the entire lattice *
   * \param state	    gauge field     	       (Read)
   */
  const LinearOperator< multi1d<LatticeFermion> >* 
  UnprecHTContFrac5DFermActArray::lMdagM(Handle<const ConnectState> state) const
  {
    return new lmdagm<multi1d<LatticeFermion> >(linOp(state));
  }


  //! Propagator of unpreconditioned H_T kernel continued fraction (5D) operator
  /*! \ingroup qprop
   *
   * Propagator solver for Extended overlap fermions
   */
  template<typename T>
  class OvHTCFZ5DQprop : public SystemSolver<T>
  {
  public:
    //! Constructor
    /*!
     * \param A_        Linear operator ( Read )
     * \param Mass_      quark mass ( Read )
     */
    OvHTCFZ5DQprop(Handle< const LinearOperator< multi1d<T> > > A_,
		   Handle< LinearOperator<LatticeFermion> > D_denum_,
		   const Real& Mass_,
		   const InvertParam_t& invParam_) : 
      A(A_), D_denum(D_denum_), Mass(Mass_), invParam(invParam_) {}

    //! Destructor is automatic
    ~OvHTCFZ5DQprop() {}

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
      multi1d<LatticeFermion> chi5(N5);
      multi1d<LatticeFermion> psi5(N5);

      
      // For reasons I do not appreciate doing the solve as
      //  M^{dag} M psi = M^{dag} chi
      //  seems a few iterations faster and more accurate than
      //  Doing M^{dag} M psi = chi
      //  and then applying M^{dag}

      // So first get  M^{dag} gamma_5 chi into the source.
    
      switch(invParam.invType) 
      {
      case CG_INVERTER: 
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
	
	// set up D_denum^dag gamma_5 chi into chi-1
	{
	  LatticeFermion tmp = Gamma(G5)*chi;
	  (*D_denum)(chi5[N5-1], tmp, MINUS);
	}
	(*A)(tmp5, chi5, MINUS);
	
	// psi5 = (M^{dag} M)^(-1) M^{dag} * gamma_5 * chi5
	// psi5[N5]  = (1 - m)/2 D^{-1}(m) chi [N5]
	InvCG2(*A, tmp5, psi5, invParam.RsdCG, invParam.MaxCG, n_count);
      }
      break;
      
      default:
	QDP_error_exit("Unknown inverter type", invParam.invType);
      }
  
      if ( n_count == invParam.MaxCG )
	QDP_error_exit("no convergence in the inverter", n_count);
    
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

      return n_count;
    }

  private:
    // Hide default constructor
    OvHTCFZ5DQprop() {}

    Handle< const LinearOperator< multi1d<T> > > A;
    Handle< LinearOperator<LatticeFermion> > D_denum;  
    const Real Mass;
    const InvertParam_t invParam;
  };

  
  //! Propagator of unpreconditioned H_T kernel continued fraction (5D) operator
  const SystemSolver<LatticeFermion>* 
  UnprecHTContFrac5DFermActArray::qprop(Handle<const ConnectState> state,
					const InvertParam_t& invParam) const
  {
    const multi1d<LatticeColorMatrix>& u = state->getLinks();
    Handle< const LinearOperator<LatticeFermion> > D_w(new UnprecWilsonLinOp(u, Real(-OverMass)));
    Handle< const LinearOperator<LatticeFermion> > D_denum(new UnprecDWFTransfDenLinOp(u, a5, D_w));

    return new OvHTCFZ5DQprop<LatticeFermion>(Handle< const LinearOperator< multi1d<LatticeFermion> > >(linOp(state)),
					      D_denum,
					      getQuarkMass(),
					      invParam);
  }


} // End namespace Chroma
