// $Id: zolotarev5d_fermact_array_w.cc,v 1.4 2004-01-13 17:52:15 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned extended-Overlap (5D) (Naryanan&Neuberger) action
 */

#include "chromabase.h"

#include "actions/ferm/fermacts/zolotarev5d_fermact_array_w.h"
#include "actions/ferm/linop/zolotarev5d_linop_array_w.h"
#include "actions/ferm/linop/lmdagm.h"
#include "actions/ferm/invert/invcg2_array.h"
#include "zolotarev.h"


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
  Write(my_writer, beta);
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
  Write(my_writer, beta);
  Write(my_writer, alpha);
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

