// $Id: eoprec_ovext_fermact_array_w.cc,v 3.2 2007-02-22 21:11:45 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned extended-Overlap (5D) (Naryanan&Neuberger) action
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/eoprec_ovext_fermact_array_w.h"
#include "actions/ferm/linop/eoprec_ovext_linop_array_w.h"
#include "actions/ferm/invert/invcg2_array.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermstates/ferm_createstate_reader_w.h"

#include "io/enum_io/enum_io.h"
#include "io/overlap_state_info.h"
#include "actions/ferm/fermacts/zolotarev.h"

#include "actions/ferm/linop/unprec_wilson_linop_w.h"
#include "actions/ferm/linop/unprec_dwftransf_linop_w.h"

#include "actions/ferm/invert/syssolver_cg_params.h"

namespace Chroma
{
  //! Hooks to register the class with the fermact factory
  namespace EvenOddPrecOvExtFermActArrayEnv
  {
    //! Callback function
    WilsonTypeFermAct5D<LatticeFermion,
			multi1d<LatticeColorMatrix>,
			multi1d<LatticeColorMatrix> >* createFermAct5D(XMLReader& xml_in,
								       const std::string& path)
    {
      return new EvenOddPrecOvExtFermActArray(CreateFermStateEnv::reader(xml_in, path), 
					      EvenOddPrecOvExtFermActArrayParams(xml_in, path));
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
    const std::string name = "OVEXT";

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
  }


  //! Read parameters
  EvenOddPrecOvExtFermActArrayParams::EvenOddPrecOvExtFermActArrayParams(XMLReader& xml, 
									 const std::string& path)
  {
    try 
    {
      XMLReader paramtop(xml, path);
      
      // Read the stuff for the action
      read(paramtop, "OverMass", OverMass);
      read(paramtop, "b5", b5);
      read(paramtop, "c5", c5);
      read(paramtop, "Mass", Mass);
      read(paramtop, "RatPolyDeg", RatPolyDeg);
      read(paramtop, "ApproximationType", approximation_type);
      if (approximation_type == COEFF_TYPE_ZOLOTAREV) {
	read(paramtop, "ApproxMin", ApproxMin);
	read(paramtop, "ApproxMax", ApproxMax);
      }
      else {
	ApproxMin = ApproxMax = 0.0;
      }

      XMLReader tuning_strategy_reader(paramtop, "TuningStrategy");
      std::ostringstream os;
      tuning_strategy_reader.print(os);
      tuning_strategy_xml = os.str();    
    }
    catch( const std::string& e) { 
      QDPIO::cout << "Caught Exception while reading XML: " << e << endl;
      QDP_abort(1);
    }
  }


  //! Read parameters
  void read(XMLReader& xml, const string& path, EvenOddPrecOvExtFermActArrayParams& param)
  {
    EvenOddPrecOvExtFermActArrayParams tmp(xml, path);
    param = tmp;
  }

  void write(XMLWriter& xml, const string& path, const EvenOddPrecOvExtFermActArrayParams& p) 
  {
    push(xml, path);
    write(xml, "OverMass", p.OverMass);
    write(xml, "b5" , p.b5);
    write(xml, "c5" , p.c5);
    write(xml, "Mass", p.Mass);
    write(xml, "RatPolyDeg", p.RatPolyDeg);
    write(xml, "ApproximationType", p.approximation_type);
    if (p.approximation_type == COEFF_TYPE_ZOLOTAREV) {
      write(xml, "ApproxMin", p.ApproxMin);
      write(xml, "ApproxMax", p.ApproxMax);
    }

    //    This may be broken here...
    QDP::write(xml, "TuningStrategy", p.tuning_strategy_xml);

    pop(xml);
  }


  // General constructor
  EvenOddPrecOvExtFermActArray::EvenOddPrecOvExtFermActArray(Handle< CreateFermState<T,P,Q> > cfs_, 
							     const EvenOddPrecOvExtFermActArrayParams& param_) :
    cfs(cfs_), param(param_) 
  {
    // Set up the strategy for tuning the betas
    std::istringstream ts_is(param.tuning_strategy_xml);
    XMLReader tuning_xml(ts_is);
    std::string strategy_name;
    try { 
      read(tuning_xml, "/TuningStrategy/Name", strategy_name);
    }
    catch(const std::string& e) { 
      QDPIO::cerr << "Caught exception processing TuningStrategy: " << e << endl;
    }
      
      
    theTuningStrategy = TheAbsOvExtTuningStrategyFactory::Instance().createObject(strategy_name, tuning_xml, "/TuningStrategy");
  }


  // Initializer
  int EvenOddPrecOvExtFermActArray::getN5FromRatPolyDeg(const int& RatPolyDeg) const 
  {
    // Type 0 and Tanh approximations: 

    // If RatPolyDeg is even: => 2*(RatPolyDeg/2) + 1 = RatPolyDeg+1
    // If RatPolyDeg is odd: =>  2*((RatPolyDeg-1)/2 + 1 = RatPolyDeg
    if( RatPolyDeg % 2 == 0 ) { 
      return RatPolyDeg+1; 
    }
    else { 
      return RatPolyDeg;
    }
  }


  //! Get the rational approximation coefficients
  void EvenOddPrecOvExtFermActArray::init(int& Npoles, 
					  Real& coeffP, 
					  multi1d<Real>& resP,
					  multi1d<Real>& rootQ) const
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


    /* Hermitian 4D overlap operator 1/2 ( 1 + Mass + (1 - Mass) gamma5 * sgn(H)) 
       using a partial fraction expansion of the optimal rational function
       approximation to sgn. Here, H = 1/2 * gamma_5 * (1/kappa - D'). 
       The coefficients are computed by Zolotarev's formula. */

    switch(param.approximation_type) { 
    case COEFF_TYPE_ZOLOTAREV:
    {
      // Rescale approximation: (approxMin, approxMax)
      //                  ->(alpha*approxMin, alpha*approxMax)
      scale_fac = Real(1) / (param.ApproxMax);
      eps = (param.ApproxMin) / (param.ApproxMax);
	
      QDPIO::cout << "Initing Linop with Zolotarev Coefficients" << endl;
	
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
      rdata = zolotarev(toFloat(eps), param.RatPolyDeg, type);
      if( rdata == 0x0 ) { 
	QDPIO::cerr << "Failed to get Zolo Coeffs" << endl;
	QDP_abort(1);
      } 
    }
    break;

    case COEFF_TYPE_TANH_UNSCALED:
    {
      scale_fac = Real(1) ;
      eps = param.ApproxMin;
	
      QDPIO::cout << "Initing Linop with Unscaled Higham Rep tanh Coefficients" << endl;
	
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
      rdata = higham(toFloat(eps), param.RatPolyDeg);
    }
    break;

    default:
      // The map system should ensure that we never get here but 
      // just for style
      QDPIO::cerr << "Unknown coefficient type: " << param.approximation_type
		  << endl;
      QDP_abort(1);
    }
    
    maxerr = (Real)(rdata -> Delta);
    QDPIO::cout << "Maxerr " << maxerr << flush << endl; 

    /* The number of residuals and poles */
    /* Allocate the roots and residua */
    Npoles = rdata -> dd;

    if ( (2*Npoles+1) != getN5FromRatPolyDeg(param.RatPolyDeg)) { 
      QDPIO::cout << "Oops. 2Npoles+1 = " << (2*Npoles+1)
		  << " but N5=" << getN5FromRatPolyDeg(param.RatPolyDeg)
		  << " this is inconsitent" << endl;
      QDP_abort(1);
    }


    /* The roots, i.e., the shifts in the partial fraction expansion */
    rootQ.resize(Npoles);
    /* The residuals in the partial fraction expansion */
    resP.resize(Npoles);
 
    /* The coefficients from the partial fraction.
       -- reverse order so biggest is near the physical field */

    /* Fill in alpha[0] = alpha[da] if it is not zero*/
    coeffP = rdata -> alpha[rdata -> da - 1] * scale_fac;
    /* Fill in the coefficients for the roots and the residua */
    /* Make sure that the smallest shift is in the last value rootQ(Npoles-1)*/
    Real t = Real(1) / (scale_fac * scale_fac);
    for(int n=0; n < Npoles; ++n) {
    
      resP[Npoles-1-n] = rdata -> alpha[n] / scale_fac;
      rootQ[Npoles-1-n] = rdata -> ap[n];
      rootQ[Npoles-1-n] = -(rootQ[Npoles-1-n] * t);
    }
    
    QDPIO::cout << "PartFracApprox n=" << param.RatPolyDeg 
		<<" scale=" << scale_fac
		<<" Mass=" << param.Mass
		<< endl;
  
    QDPIO::cout << "Approximation on [-1,-eps] U [eps,1] with eps = " << eps <<
      endl;
    QDPIO::cout << "Maximum error |R(x) - sqn(x)| <= " << maxerr << endl;
  
    switch( param.approximation_type) {
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

    case COEFF_TYPE_TANH_UNSCALED:
      QDPIO::cout << "Coefficients from Unscaled Higham Tanh representation" << endl;
      break;

    default:
      QDPIO::cerr << "Unknown coefficient type " << param.approximation_type 
		  << endl;
      break;
    }

    QDPIO::cout << "Number of poles= " << Npoles << endl;
    QDPIO::cout << "Overall Factor=  " << coeffP << endl;
    QDPIO::cout << "Numerator coefficients:" << endl;
    for(int n=0; n < Npoles; n++) { 
      QDPIO::cout <<"  resP[" << n << "]= " << resP[n] << endl;
    }
    QDPIO::cout << "Denominator roots: " << endl;
    for(int n=0; n < Npoles; n++) { 
      QDPIO::cout <<"  rootQ[" << n<< "]= " << rootQ[n] << endl;
    }
 
    // Free the arrays allocate by Tony's zolo
    zolotarev_free(rdata);
  }


  //! Produce a linear operator for this action
  /*! 
   * The operator acts on the entire lattice
   *
   * \param state	    gauge field     	       (Read)
   */
  EvenOddPrecConstDetLinearOperatorArray<LatticeFermion, 
					 multi1d<LatticeColorMatrix>,
					 multi1d<LatticeColorMatrix> >* 
  EvenOddPrecOvExtFermActArray::linOp(Handle< FermState<T,P,Q> > state) const
  {
    int Npoles;
    Real coeffP;
    multi1d<Real> resP;
    multi1d<Real> rootQ;
    
    // Get the coefficients
    
    init(Npoles, coeffP, resP, rootQ);

    multi1d<Real> beta(Npoles);
    (*theTuningStrategy)(beta, coeffP, resP, rootQ, param.Mass);

    return new EvenOddPrecOvExtLinOpArray(state,
					  Npoles,
					  coeffP, 
					  resP, 
					  rootQ,
					  beta,
					  param.OverMass,
					  param.Mass,
					  param.b5,
					  param.c5);
    
  }


  //! Propagator of an un-preconditioned Extended-Overlap linear operator
  /*! \ingroup qprop
   *
   * Propagator solver for Extended overlap fermions
   */
  template<typename T, typename P, typename Q>
  class PrecOvExt5DQprop : public SystemSolver<T>
  {
  public:
    //! Constructor
    /*!
     * \param A_         Linear operator ( Read )
     * \param Mass_      quark mass ( Read )
     */
    PrecOvExt5DQprop(Handle< EvenOddPrecConstDetLinearOperatorArray<T,P,Q> > A_,
		     Handle< LinearOperator<T> > D_denum_,
		     const Real& Mass_,
		     const SysSolverCGParams& invParam_) : 
      A(A_), D_denum(D_denum_), Mass(Mass_),  invParam(invParam_) {}

    //! Destructor is automatic
    ~PrecOvExt5DQprop() {}
    
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
	multi1d<T> chi5(N5);
	multi1d<T> psi5(N5);

//	if( invType == "CG_INVERTER") 
	{
	  multi1d<T> tmp5_1(N5);
	  {
	    multi1d<T> tmp5_2(N5);
	    multi1d<T> tmp5_3(N5);

	    chi5 = zero;
	    psi5 = zero;
	    tmp5_1 = zero;

	    // Need to prepare the source 
	    psi5[N5-1] = Gamma(G5)*chi;

	    A->evenEvenInvLinOp(tmp5_2, psi5, PLUS);
	    A->oddEvenLinOp(tmp5_3, tmp5_2, PLUS);


	    // chi5_e = S_e
	    // chi5_o = S_o - QoeQee^{-1} S_e
	    for(int i=0; i < N5; i++) { 
	      chi5[i][rb[0]] = psi5[i];
	      chi5[i][rb[1]] = psi5[i] - tmp5_3[i];
	    }
	  }  // tmp5_2 and tmp5_3 go away

	
	  // CGNE tmp5_1 holds source
	  (*A)(tmp5_1, chi5, MINUS);

	  // Initial guess
	  psi5[N5-1][rb[1]] = psi;     

	  // Solve M^{+}M psi = M^{+} chi
	  res = InvCG2(*A, tmp5_1, psi5, invParam.RsdCG, invParam.MaxCG);
        
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
//	else
//	{
//	  QDPIO::cerr << EvenOddPrecOvExtFermActArrayEnv::name 
//		      << "Unsupported inverter type =" << invParam.invType << endl;
//	  QDP_abort(1);
//	}
  
	if ( res.n_count == invParam.MaxCG )
	  QDP_error_exit("no convergence in the inverter", res.n_count);
    
	// Need to compute residual someday...
	res.resid = zero;

	// Solution now lives in psi5
	{
	  LatticeFermion tmp4;
	
	  // Take care of H_t scaling
	  (*D_denum)(tmp4, psi5[N5-1], PLUS);
	  psi5[N5-1] = tmp4;
	}

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

	return res;
      }

  private:
    // Hide default constructor
    PrecOvExt5DQprop() {}

    Handle< EvenOddPrecConstDetLinearOperatorArray<T,P,Q> > A;
    Handle< LinearOperator<T> > D_denum;
    Real Mass;
    SysSolverCGParams invParam;
  };

 
  //! Propagator of an un-preconditioned Extended-Overlap linear operator
  SystemSolver<LatticeFermion>* 
  EvenOddPrecOvExtFermActArray::qprop(Handle< FermState<T,P,Q> > state,
				      const GroupXML_t& invParam) const
  {
    Real a5 = param.b5- param.c5;

    Handle< LinearOperator<LatticeFermion> > D_w(new UnprecWilsonLinOp(state, -param.OverMass));
    Handle< LinearOperator<LatticeFermion> > D_denum(new UnprecDWFTransfDenLinOp(a5, D_w));
    Handle< EvenOddPrecConstDetLinearOperatorArray<T,P,Q> > A(linOp(state));

    std::istringstream  is(invParam.xml);
    XMLReader  paramtop(is);
	
    return new PrecOvExt5DQprop<T,P,Q>(A, D_denum, param.Mass, SysSolverCGParams(paramtop,invParam.path));
  }

}
