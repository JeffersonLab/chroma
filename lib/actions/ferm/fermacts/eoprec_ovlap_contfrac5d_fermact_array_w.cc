// $Id: eoprec_ovlap_contfrac5d_fermact_array_w.cc,v 3.2 2007-02-22 21:11:45 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned extended-Overlap (5D) (Naryanan&Neuberger) action
 */

#include "chromabase.h"

#include "actions/ferm/fermacts/unprec_wilson_fermact_w.h"
#include "actions/ferm/fermacts/eoprec_ovlap_contfrac5d_fermact_array_w.h"
//#include "actions/ferm/linop/eoprec_ovlap_contfrac5d_linop_array_w.h"
#include "actions/ferm/linop/ovlap_contfrac5d_w.h"
#include "actions/ferm/linop/eoprec_ovlap_contfrac5d_pv_linop_array_w.h"
#include "actions/ferm/invert/invcg2_array.h"
#include "actions/ferm/fermacts/zolotarev.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermstates/ferm_createstate_reader_w.h"

#include "io/enum_io/enum_io.h"
#include "io/overlap_state_info.h"

#include "actions/ferm/invert/syssolver_cg_params.h"

namespace Chroma
{

  //! Hooks to register the class with the fermact factory
  namespace EvenOddPrecOvlapContFrac5DFermActArrayEnv
  {
    //! Callback function
    WilsonTypeFermAct5D<LatticeFermion, 
			multi1d<LatticeColorMatrix>,
			multi1d<LatticeColorMatrix> >* createFermAct5D(XMLReader& xml_in,
								       const std::string& path)
    {
      return new EvenOddPrecOvlapContFrac5DFermActArray(CreateFermStateEnv::reader(xml_in, path), 
							EvenOddPrecOvlapContFrac5DFermActParams(xml_in, path));
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
    const std::string name = "OVERLAP_CONTINUED_FRACTION_5D";
    
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
  } // End Namespace EvenOddPrecOvlapContFrac5DFermActArrayEnv

  
  //! Read XML
  EvenOddPrecOvlapContFrac5DFermActParams::EvenOddPrecOvlapContFrac5DFermActParams(XMLReader& xml, 
										   const std::string& path)
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
    }
    catch( const string &e ) {
      QDPIO::cerr << "Caught Exception reading prec ContFrac Fermact params: " << e << endl;
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
    write(xml_out, "ApproxMin", p.ApproxMin);
    write(xml_out, "ApproxMax", p.ApproxMax);
    
    pop(xml_out);
    
    if( path != "." ) { 
      pop(xml_out);
    }
  }

  
  
  // Construct the action out of a parameter structure
  EvenOddPrecOvlapContFrac5DFermActArray::EvenOddPrecOvlapContFrac5DFermActArray(
    Handle< CreateFermState<T,P,Q> > cfs_a_, 
    const EvenOddPrecOvlapContFrac5DFermActParams& params_) :
    cfs(cfs_a_), params(params_) 
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
					       multi1d<Real>& beta) const
  {
    int type = 0;
    zolotarev_data *rdata;
    Real epsilon;

    Real approxMin = params.ApproxMin;
    Real approxMax = params.ApproxMax;

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
    
    QDPIO::cout << "EvenOddPrecOvlapContfrac5d: " << endl
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
  EvenOddPrecConstDetLinearOperatorArray<LatticeFermion, 
					 multi1d<LatticeColorMatrix>,
					 multi1d<LatticeColorMatrix> >* 
  EvenOddPrecOvlapContFrac5DFermActArray::linOp(Handle< FermState<T,P,Q> > state_) const
  {
    START_CODE();
      
    multi1d<Real> alpha;
    multi1d<Real> beta;
    Real scale_factor;
      
    init(scale_factor, alpha, beta);
      
    return new EvenOddPrecOvlapContFrac5DLinOpArray(state_,
						    params.Mass,
						    params.OverMass,
						    N5,
						    scale_factor,
						    alpha,
						    beta,
						    isLastZeroP);
  }
  

  //! Produce a Pauli-Villars linear operator for this action
  /*!
   * \param state	    gauge field     	       (Read)
   */
  EvenOddPrecConstDetLinearOperatorArray<LatticeFermion, 
					 multi1d<LatticeColorMatrix>,
					 multi1d<LatticeColorMatrix> >* 
  EvenOddPrecOvlapContFrac5DFermActArray::linOpPV(Handle< FermState<T,P,Q> > state_) const
  {
    multi1d<Real> alpha;
    multi1d<Real> beta;
    Real scale_factor;
      
    init(scale_factor, alpha, beta);
      
    // Hmm, not sure about what all the rescaling does to the PV....
    return new EvenOddPrecOvlapContFrac5DPVLinOpArray(state_,
						      params.Mass,
						      params.OverMass,
						      N5,
						      scale_factor,
						      alpha,
						      beta,
						      isLastZeroP);
  }
  
  

  //! Propagator of an un-preconditioned Extended-Overlap linear operator
  /*! \ingroup qprop
   *
   * Propagator solver for DWF-like fermions
   */
  template<typename T, typename P, typename Q>
  class ContFrac5DQprop : public SystemSolver<T>
  {
  public:
    //! Constructor
    /*!
     * \param A_         Linear operator ( Read )
     * \param Mass_      quark mass ( Read )
     */
    ContFrac5DQprop(Handle< EvenOddPrecConstDetLinearOperatorArray<T,P,Q> > A_,
		    const Real& Mass_,
		    const SysSolverCGParams& invParam_) : 
      A(A_), Mass(Mass_), invParam(invParam_) {}

    //! Destructor is automatic
    ~ContFrac5DQprop() {}

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

//      if( invType == "CG_INVERTER") 
      {
	multi1d<T> tmp5_1(N5);
	{
	  multi1d<T> tmp5_2(N5);
	  multi1d<T> tmp5_3(N5);

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
//      else
//      {
//	QDPIO::cerr << EvenOddPrecOvlapContFrac5DFermActArrayEnv::name 
//		    << "Unsupported inverter type =" << invType << endl;
//	QDP_abort(1);
//     }
  
      if ( res.n_count == invParam.MaxCG )
	QDP_error_exit("no convergence in the inverter", res.n_count);
    
      // Compute residual
      {
	multi1d<T>  r(N5);
	A->unprecLinOp(r, psi5, PLUS);
	chi5 = zero;
	chi5[N5-1] = Gamma(G5)*chi;
	r -= chi5;
	res.resid = sqrt(norm2(r));
      }

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

      return res;
    }

  private:
    // Hide default constructor
    ContFrac5DQprop() {}

    Handle< EvenOddPrecConstDetLinearOperatorArray<T,P,Q> > A;
    Real Mass;
    SysSolverCGParams invParam;
  };

 
  //! Propagator of an un-preconditioned Extended-Overlap linear operator
  SystemSolver<LatticeFermion>* 
  EvenOddPrecOvlapContFrac5DFermActArray::qprop(Handle< FermState<T,P,Q> > state,
						const GroupXML_t& invParam) const
  {
    std::istringstream  is(invParam.xml);
    XMLReader  paramtop(is);
	
    return new ContFrac5DQprop<T,P,Q>(
      Handle< EvenOddPrecConstDetLinearOperatorArray<T,P,Q> >(linOp(state)),
      getQuarkMass(),
      SysSolverCGParams(paramtop,invParam.path));
  }

} // End namespace Chroma
