// $Id: zolotarev5d_fermact_array_w.cc,v 1.8 2004-04-26 11:19:12 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned extended-Overlap (5D) (Naryanan&Neuberger) action
 */

#include "chromabase.h"

#include "actions/ferm/fermacts/unprec_wilson_fermact_w.h"
#include "actions/ferm/fermacts/zolotarev5d_fermact_array_w.h"
#include "actions/ferm/linop/zolotarev5d_linop_array_w.h"
#include "actions/ferm/linop/lmdagm.h"
#include "actions/ferm/invert/invcg2_array.h"
#include "zolotarev.h"

  // Construct the action out of a parameter structure
Zolotarev5DFermActArray::Zolotarev5DFermActArray(Handle< FermBC< multi1d< LatticeFermion> > > fbc_a_, 
						  Handle< FermBC< LatticeFermion > > fbc_,
						 const Zolotarev5DFermActParams& params,
						 XMLWriter& writer_) :
  fbc(fbc_a_), m_q(params.Mass), RatPolyDeg(params.RatPolyDeg), writer(writer_)  {
  
  // Check RatPolyDeg is even
  if ( params.RatPolyDeg % 2 == 0 ) { 
    QDP_error_exit("For Now (and possibly forever), 5D Operators can only be constructed with ODD approximation order. You gave an even one: =%d\n", params.RatPolyDeg);
  }
  N5 = params.RatPolyDeg;
  
  // Get the auxiliary fermion action
  UnprecWilsonTypeFermAct<LatticeFermion>* S_w;
  switch( params.AuxFermActHandle->getFermActType() ) {
  case FERM_ACT_WILSON:
    {
      // Upcast
      const WilsonFermActParams& wils = dynamic_cast<const WilsonFermActParams &>( *(params.AuxFermActHandle));
      
      //Get the FermAct
      S_w = new UnprecWilsonFermAct(fbc_, wils.Mass);
      if( S_w == 0x0 ) { 
	QDPIO::cerr << "Unable to instantiate S_aux " << endl;
	QDP_abort(1);
      }
      
    }
    break;
  default:
    QDPIO::cerr << "Auxiliary Fermion Action Unsupported" << endl;
    QDP_abort(1);
  }
  
  // Drop AuxFermAct into a Handle immediately.
  // This should free things up at the end
  Handle<UnprecWilsonTypeFermAct<LatticeFermion> >  Handle_S_w(S_w);
  S_aux = Handle_S_w;
}

void
Zolotarev5DFermActArray::init(Real& scale_fac,
			      multi1d<Real>& alpha,
			      multi1d<Real>& beta,
			      int& NEig,
			      multi1d<Real>& EigValFunc,
			      const OverlapConnectState<LatticeFermion>& state) const
{
  XMLBufferWriter my_writer;
  push( my_writer, "Zolo5DInit" );

  int NEigVal = state.getEigVal().size();
  if( NEigVal == 0 ) {
    NEig = 0;
  }
  else {
    NEig = NEigVal - 1;
  }

  scale_fac = Real(1) / state.getApproxMax();
  Real eps = state.getApproxMin() * scale_fac;
  

  int type = 0;
  zolotarev_data *rdata = zolotarev(toFloat(eps), RatPolyDeg, type);

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


  push(my_writer, "ZoloContFracCoeff");
  write(my_writer, "beta", beta);
  pop(my_writer);

  QDPIO::cout << "Did Beta" << endl;

  alpha.resize(N5);
  for(int i = 0; i < N5; i++) {
    alpha[i] = Real(1);
  }

  alpha[N5-1] = rdata->beta[N5-1];

  QDPIO::cout << "Did alpha" << endl;

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

  QDPIO::cout << "Did alpha2 " << endl;

  push(my_writer, "ZoloEquivTransCoeff");
  write(my_writer, "beta", beta);
  write(my_writer, "alpha", alpha);
  pop(my_writer);

  QDPIO::cout << "Zolotarev 5d: " 
	      << " N5=" << N5 << " scale=" << scale_fac
	      << " Nwils = " << NEigVal << " m_q=" << m_q << endl ;
  QDPIO::cout << "Approximation on [-1,eps] U [eps,1] with eps = " << eps <<endl;
 
  QDPIO::cout << "Maximum error | R(x) - sgn(x) | <= Delta = " << maxerr << endl;

  if( type == 0 ) {
    QDPIO::cout << "Approximation type " << type << " with R(0) = 0" << endl;
  }
  else {
    QDPIO::cout << "Approximation type " << type << " with R(0) = infinity" << endl;
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

  pop(my_writer);

  writer << my_writer;

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
Zolotarev5DFermActArray::linOp(Handle<const ConnectState> state_) const
{
  START_CODE("Zolotarev5DFermActArray::linOp");
   const OverlapConnectState<LatticeFermion>& state = dynamic_cast<const OverlapConnectState<LatticeFermion>&>(*state_);

  if (state.getEigVec().size() != state.getEigVal().size())
    QDP_error_exit("Zolotarev5DFermActArray: inconsistent sizes of eigenvectors and values");

  int NEigVal = state.getEigVal().size();
  int NEig;

  multi1d<Real> alpha;
  multi1d<Real> beta;
  Real scale_factor;
  multi1d<Real> EigValFunc(NEigVal);

  init(scale_factor, alpha, beta, NEig, EigValFunc, state);

  return new Zolotarev5DLinOpArray( *S_aux,
				    state_,
				    m_q,
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
Zolotarev5DFermActArray::lMdagM(Handle<const ConnectState> state) const
{
  return new lmdagm<multi1d<LatticeFermion> >(linOp(state));
}


const OverlapConnectState<LatticeFermion>*
Zolotarev5DFermActArray::createState(const multi1d<LatticeColorMatrix>& u_,
				  const Real& approxMin_) const 
{
  if ( toBool( approxMin_ < Real(0) )) { 
    ostringstream error_str;
    error_str << "Zolotarev5DFermActArray: approxMin_ has to be positive" << endl;
    throw error_str.str();
  }
 

  // First put in the BC
  multi1d<LatticeColorMatrix> u_tmp = u_;
  getFermBC().modifyU(u_tmp);

  Real approxMax = Real(2)*Real(Nd);
  return new OverlapConnectState<LatticeFermion>(u_tmp, approxMin_, approxMax);
}


const OverlapConnectState<LatticeFermion>*
Zolotarev5DFermActArray::createState(const multi1d<LatticeColorMatrix>& u_,
				  const Real& approxMin_,
				  const Real& approxMax_) const
{
  ostringstream error_str;
  
 
  if ( toBool(approxMin_ < 0 )) { 
    error_str << "Zolotarev5DFermActArray: approxMin_ has to be positive" << endl;
    throw error_str.str();
  }

  if ( toBool(approxMax_ < approxMin_) ) { 
    error_str << "Zolotarev5DFermActArray: approxMax_ has to be larger than approxMin_" << endl;
    throw error_str.str();
  }
 

  // First put in the BC
  multi1d<LatticeColorMatrix> u_tmp = u_;
  getFermBC().modifyU(u_tmp);

  return new OverlapConnectState<LatticeFermion>(u_tmp, approxMin_, approxMax_);
}



const OverlapConnectState<LatticeFermion>*
Zolotarev5DFermActArray::createState(const multi1d<LatticeColorMatrix>& u_,
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
Zolotarev5DFermActArray::qprop(LatticeFermion& psi, 
			       Handle<const ConnectState> state, 
			       const LatticeFermion& chi, 
			       enum InvType invType,
			       const Real& RsdCG, 
			       int MaxCG, int& ncg_had) const
{

  START_CODE("Zolotarev5DFermActArray::qprop");

  const int  N5 = size();   // array size better match
  const Real m_q = quark_mass();
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
  psi5[N5-1] *= Real(2)/(Real(1)-m_q);
  
  // Remove contact term
  psi = psi5[N5-1] - chi;

  // Overall normalization
  Real ftmp1 = Real(1) / (Real(1) - m_q);
  psi *= ftmp1;

  END_CODE("Zolotarev5DFermActArray::qprop");
}


const OverlapConnectState<LatticeFermion>*
Zolotarev5DFermActArray::createState(const multi1d<LatticeColorMatrix>& u_,
				const OverlapStateInfo& state_info,
				XMLWriter& xml_out,
				Real wilsonMass) const
{
  push(xml_out, "Zolo4DCreateState");


  // If No eigen values specified use min and max
 if ( state_info.getNWilsVec() == 0 ) { 
    write(xml_out, "ApproxMin", state_info.getApproxMin());
    write(xml_out, "ApproxMax", state_info.getApproxMax());
    pop(xml_out);

    return createState(u_,
		       state_info.getApproxMin(),
		       state_info.getApproxMax());
  }
  else {

    QDPIO::cout << "Warning!!!! Using 5D op with projected e-values is dubious "
		<< " at this time " << endl;

    // If there are eigen values, either load them, 
    if( state_info.loadEigVec() ) {
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
	  QDPIO::cerr << "WilsonMass unspecified, or | WilsonMass | > 8" << endl;
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

      pop(xml_out); // Zolo4DCreateState
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
}
