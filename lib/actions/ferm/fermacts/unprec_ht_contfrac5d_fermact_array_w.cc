/*! \file
 *  \brief Unpreconditioned H_T kernel continued fraction (5D) action
 */

#include "chromabase.h"

#include "actions/ferm/fermacts/unprec_wilson_fermact_w.h"
#include "actions/ferm/fermacts/unprec_ht_contfrac5d_fermact_array_w.h"

#include "actions/ferm/linop/unprec_ht_contfrac5d_linop_array_w.h"
#include "actions/ferm/linop/unprec_wilson_linop_w.h"
#include "actions/ferm/linop/unprec_dwftransf_linop_w.h"

#include "actions/ferm/invert/invcg2_array.h"
#include "zolotarev.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermstates/ferm_createstate_reader_w.h"

#include "io/enum_io/enum_io.h"

#include "actions/ferm/invert/syssolver_cg_params.h"

namespace Chroma
{

  //! Hooks to register the class with the fermact factory
  namespace UnprecHTContFrac5DFermActArrayEnv
  {
    //! Callback function
    WilsonTypeFermAct5D<LatticeFermion,
			multi1d<LatticeColorMatrix>,
			multi1d<LatticeColorMatrix> >* createFermAct5D(XMLReader& xml_in,
								       const std::string& path)
    {
      return new UnprecHTContFrac5DFermActArray(CreateFermStateEnv::reader(xml_in, path), 
						UnprecHTContFrac5DFermActParams(xml_in, path));
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
    const std::string name = "UNPRECONDITIONED_HT_CONTINUED_FRACTION_5D";
    
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
  } // End Namespace UnprecHTContFrac5DFermActArrayEnv

  
  //! Read XML
  UnprecHTContFrac5DFermActParams::UnprecHTContFrac5DFermActParams(XMLReader& xml, const std::string& path)
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
	QDPIO::cout << "b5 and c5 not specified. Using Shamir values: b5=1 c5=0" << std::endl;
	b5 = Real(1);
	c5 = Real(0);
      }
      else {
	read(in, "b5", b5);
	read(in, "c5", c5);
      }

    }
    catch( const std::string &e ) {
      QDPIO::cerr << "Caught Exception reading Unprec HT ContFrac Fermact params: " << e << std::endl;
      QDP_abort(1);
    }
  }


  void read(XMLReader& xml_in, const std::string& path,
	    UnprecHTContFrac5DFermActParams& param) 
  {
    UnprecHTContFrac5DFermActParams tmp(xml_in, path);
    param = tmp;
  }

  
  void write(XMLWriter& xml_out, const std::string& path, const UnprecHTContFrac5DFermActParams& p)
  {
    if ( path != "." ) { 
      push( xml_out, path);
    }
    
    write(xml_out, "OverMass", p.OverMass);
    write(xml_out, "Mass", p.Mass);
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
  UnprecHTContFrac5DFermActArray::UnprecHTContFrac5DFermActArray(
    Handle< CreateFermState<T,P,Q> > cfs_a_, 
    const UnprecHTContFrac5DFermActParams& params_) :
    cfs(cfs_a_), params(params_) 
  {
    // WHAT IS BELOW ONLY WORKS FOR TYPE=0 approximations
    // which is what we use. Forget TYPE=1
    // the Tanh approximation (Higham) is of type TYPE=0
    // We have two cases.
    bool isEvenRatPolyDeg = ( params.RatPolyDeg % 2 == 0);

    if( isEvenRatPolyDeg ) 
    {
      N5 = params.RatPolyDeg+1;
      isLastZeroP = true;     
    }
    else 
    {
      N5 = params.RatPolyDeg;
      isLastZeroP = false;
    }
  }


  void
  UnprecHTContFrac5DFermActArray::init(multi1d<Real>& alpha,
				       multi1d<Real>& beta) const
  {
    int type = 0;

    zolotarev_data *rdata;
    Real epsilon;
    Real scale_fac;
				       
    switch(params.approximation_type) 
    {
    case COEFF_TYPE_ZOLOTAREV:
      epsilon = params.ApproxMin / params.ApproxMax;
      QDPIO::cout << "Initing Linop with Zolotarev Coefficients: epsilon = " << epsilon << std::endl;
      rdata = zolotarev(toFloat(epsilon), params.RatPolyDeg, type);    
      scale_fac = Real(1) / params.ApproxMax;
      break;

    case COEFF_TYPE_TANH_UNSCALED:
      epsilon = params.ApproxMin;
      QDPIO::cout << "Initing Linop with Higham Rep tanh Coefficients" << std::endl;
      rdata = higham(toFloat(epsilon), params.RatPolyDeg);
      scale_fac = Real(1);
      break;

    default:
      // The std::map system should ensure that we never get here but 
      // just for style
      QDPIO::cerr << "Unknown coefficient type: " << params.approximation_type
		  << std::endl;
      QDP_abort(1);
    }
    
    Real maxerr = (Real)(rdata->Delta);

    // Check N5 is good:
    if( N5 != rdata->db ) { 
      QDPIO::cerr << "Mismatch between N5 and N5 from Coefficient Code" << std::endl;
      QDPIO::cerr << "N5 = " << N5 << " rdata->db=" << rdata->db << std::endl;
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

    
    // Rescale the diagonal terms
    for(int i=0; i < beta.size(); i++) {
      beta[i] *= scale_fac;
    }
    
    QDPIO::cout << "UnprecHTContfrac5d:" << std::endl
                << " Degree=" << params.RatPolyDeg
		<< " N5=" << N5 << " scale=" << scale_fac
		<< " Mass=" << params.Mass << std::endl ;
    QDPIO::cout << "Approximation on [-1,eps] U [eps,1] with eps = " << epsilon <<std::endl;
    
    QDPIO::cout << "Maximum error | R(x) - sgn(x) | <= Delta = " << maxerr << std::endl;
    /*
      for(int i=0; i < N5; i++) { 
      QDPIO::cout << "beta["<<i<<"] = " << beta[i] << std::endl;
      }
      for(int i=0; i < N5; i++) { 
      QDPIO::cout << "alpha["<<i<<"] = " << alpha[i] << std::endl;
      }
    */

    switch( params.approximation_type) {
    case COEFF_TYPE_ZOLOTAREV:
      QDPIO::cout << "Coefficients from Zolotarev" << std::endl;
      
      if(type == 0) {
	QDPIO::cout << "Approximation type " << type << " with R(0) = 0"
		    << std::endl;
      }
      else {
	QDPIO::cout << "Approximation type " << type << " with R(0) =  infinity"
                    << std::endl;
      }
      
      break;
    case COEFF_TYPE_TANH:
      QDPIO::cout << "Coefficients from Higham Tanh representation" << std::endl;
      break;
    case COEFF_TYPE_TANH_UNSCALED:
      QDPIO::cout << "Coefficients from Unscaled Higham Tanh representation" << std::endl;
      break;
    default:
      QDPIO::cerr << "Unknown coefficient type " << params.approximation_type 
		  << std::endl;
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
  UnprecLinearOperatorArray<LatticeFermion,
			    multi1d<LatticeColorMatrix>,
			    multi1d<LatticeColorMatrix> >* 
  UnprecHTContFrac5DFermActArray::linOp(Handle< FermState<T,P,Q> > state) const
  {
    START_CODE();

    multi1d<Real> alpha;
    multi1d<Real> beta;
      
    init(alpha, beta);
      
    return new UnprecHTContFrac5DLinOpArray(state,
					    params.OverMass,
					    params.Mass,
					    N5,
					    params.b5,
					    params.c5,
					    alpha,
					    beta,
					    isLastZeroP);
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
    OvHTCFZ5DQprop(Handle< LinearOperatorArray<T> > A_,
		   Handle< LinearOperator<T> > D_denum_,
		   const Real& Mass_,
		   const SysSolverCGParams& invParam_) : 
      A(A_), D_denum(D_denum_), Mass(Mass_), invParam(invParam_) {}

    //! Destructor is automatic
    ~OvHTCFZ5DQprop() {}

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
	
	// set up D_denum^dag gamma_5 chi into chi-1
	{
	  LatticeFermion tmp = Gamma(G5)*chi;
	  (*D_denum)(chi5[N5-1], tmp, MINUS);
	}
	(*A)(tmp5, chi5, MINUS);
	
	// psi5 = (M^{dag} M)^(-1) M^{dag} * gamma_5 * chi5
	// psi5[N5]  = (1 - m)/2 D^{-1}(m) chi [N5]
	res = InvCG2(*A, tmp5, psi5, invParam.RsdCG, invParam.MaxCG);
      }
//      else
//      {
//	QDPIO::cerr << UnprecHTContFrac5DFermActArrayEnv::name 
//		    << "Unsupported inverter type =" << invType << std::endl;
//	QDP_abort(1);
//     }
  
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
    OvHTCFZ5DQprop() {}

    Handle< LinearOperatorArray<T> > A;
    Handle< LinearOperator<T> > D_denum;  
    Real Mass;
    SysSolverCGParams invParam;
  };

  
  //! Propagator of unpreconditioned H_T kernel continued fraction (5D) operator
  SystemSolver<LatticeFermion>* 
  UnprecHTContFrac5DFermActArray::qprop(Handle< FermState<T,P,Q> > state,
					const GroupXML_t& invParam) const
  {
    Real a5 = params.b5 - params.c5;
    Real WilsonMass = -params.OverMass;

    Handle< LinearOperator<T> > D_w(new UnprecWilsonLinOp(state, WilsonMass));
    Handle< LinearOperator<T> > D_denum(new UnprecDWFTransfDenLinOp(a5, D_w));

    std::istringstream  is(invParam.xml);
    XMLReader  paramtop(is);
	
    return new OvHTCFZ5DQprop<T>(Handle< LinearOperatorArray<T> >(linOp(state)),
				 D_denum,
				 getQuarkMass(),
				 SysSolverCGParams(paramtop,invParam.path));
  }


} // End namespace Chroma
