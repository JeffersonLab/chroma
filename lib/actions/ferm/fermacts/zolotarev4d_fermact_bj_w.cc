// $Id: zolotarev4d_fermact_bj_w.cc,v 1.11 2004-01-08 03:11:29 edwards Exp $
/*! \file
 *  \brief 4D Zolotarev variant of Overlap-Dirac operator
 */
#include <iostream>
#include <sstream>
#include <chromabase.h>
#include <linearop.h>

#include "fermacts.h"
#include "actions/ferm/fermacts/zolotarev.h"
#include "actions/ferm/linop/linop.h"
#include "meas/eig/eig.h"

using namespace std;

//! Creation routine
/*! */
void 
Zolotarev4DFermActBj::init(int& numroot, 
			   Real& coeffP, 
			   multi1d<Real>& resP,
			   multi1d<Real>& rootQ, 
			   int& NEig, 
			   multi1d<Real>& EigValFunc,
			   const ZolotarevConnectState<LatticeFermion>& state ) const
{
  /* A scale factor which should bring the spectrum of the hermitian
     Wilson Dirac operator H into |H| < 1. */
  Real scale_fac;
  XMLBufferWriter my_writer;
  
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

  // The residual for the solutions of the multi-shift linear system
  // I put this in the class constructor 

  //  RsdCGinner = 1.0e-7;  // Hardwired the accuracy


  /* Hermitian 4D overlap operator 1/2 ( 1 + m_q + (1 - m_q) gamma5 * sgn(H)) 
     using a partial fraction expansion of the optimal rational function
     approximation to sgn. Here, H = 1/2 * gamma_5 * (1/kappa - D'). 
     The coefficients are computed by Zolotarev's formula. */

  int NEigVal = state.getEigVal().size();

  /* The operator gamma_5 * M with the M constructed here has its eigenvalues
     in the range m/(m + Nd) <= |gamma_5 * M| <= (m + 2*Nd)/(m + Nd) (in the 
     free case) where here m is arbitrary.
     So if we multiply M by a factor scale_fac = (m + Nd)/(m + 2*Nd) we have
     |gamma_5 * M| <= 1. */
  if(NEigVal == 0 ) 
  {
    NEig = 0;
  }
  else
  {
      NEig = NEigVal - 1;
  }

  scale_fac = 1 / state.getApproxMax();
  eps = state.getApproxMin() * scale_fac;


  push(my_writer, "Zolotarev4D");
  Write(my_writer, MaxCGinner);
  Write(my_writer, RsdCGinner);
  Write(my_writer, NEigVal);
  Write(my_writer, NEig);

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
    rdata = zolotarev(toFloat(eps), RatPolyDeg, type);
    maxerr = (Real)(rdata -> Delta);

    push(my_writer, "ZolotarevApprox");
    Write(my_writer, eps);
    Write(my_writer, scale_fac);
    Write(my_writer, RatPolyDeg);
    Write(my_writer, type);
    Write(my_writer, maxerr);
    pop(my_writer);

    /* The number of residuals and poles */
    /* Allocate the roots and residua */
    numroot = rdata -> dd;
    /* The roots, i.e., the shifts in the partial fraction expansion */
    rootQ.resize(numroot);

    /* The residuals in the partial fraction expansion */
    resP.resize(numroot);

    /* Fill in alpha[0] = alpha[da] if it is not zero*/
    coeffP = 0;
    coeffP = rdata -> alpha[rdata -> da - 1];
    /* The coefficients from the partial fraction.
       Here, we write them out for the sake of bookkeeping. */
    resP = 0;
    rootQ = 0;
    for(int n=0; n < numroot; ++n) {
      resP[n] = rdata -> alpha[n];
      rootQ[n] = rdata -> ap[n];
      rootQ[n] = -rootQ[n];
    }


    push(my_writer,"ZolotarevPartFrac");
    Write(my_writer, scale_fac);
    Write(my_writer, coeffP);
    Write(my_writer, resP);
    Write(my_writer, rootQ);
    pop(my_writer);

    /* Now fill in the coefficients for real, i.e., taking the rescaling
       into account */
    /* Fill in alpha[0] = alpha[da] if it is not zero*/
    coeffP = rdata -> alpha[rdata -> da - 1] * scale_fac;
    /* Fill in the coefficients for the roots and the residua */
    /* Make sure that the smallest shift is in the last value rootQ(numroot-1)*/
    resP = 0;
    rootQ = 0;
    Real t = Real(1) / (scale_fac * scale_fac);
    for(int n=0; n < numroot; ++n) {
      
      resP[n] = rdata -> alpha[n] / scale_fac;
      rootQ[n] = rdata -> ap[n];
      rootQ[n] = -(rootQ[n] * t);
    }
    

  /* Write them out into the namelist */
    push(my_writer,"ZolotarevPartFracResc");
    Write(my_writer, scale_fac);
    Write(my_writer, coeffP);
    Write(my_writer, resP);
    Write(my_writer, rootQ);
    pop(my_writer);

    pop(my_writer);


    QDP_info("ZOLOTAREV_4d: n= %d scale= %g coeff= %g  Nwils= %d  m_q= %g  Rsd= %g",
	     RatPolyDeg,toFloat(scale_fac),toFloat(coeffP),NEigVal,
	     toFloat(m_q),toFloat(RsdCGinner));

//  QDP_info("Auxiliary fermion action: OverAuxAct = %d",OverAuxAct);
    QDP_info("Approximation on [-1,-eps] U [eps,1] with eps = %g",toFloat(eps));
    /* maxerr = rdata -> Delta; */
    QDP_info("Maximum error |R(x) - sgn(x)| <= Delta = %g",toFloat(maxerr));
    if(type == 0) {QDP_info("Approximation type %d with R(0) = 0",type);}
    else {QDP_info("Approximation type %d with R(0) =  infinity",type);}
  
  /* We will also compute the 'function' of the eigenvalues */
  /* for the Wilson vectors to be projected out. */
  if (NEig > 0)
  {
    for(int i = 0; i < NEigVal; i++)
    {
      if (toBool(state.getEigVal()[i] > 0.0))
	EigValFunc[i] = 1.0;
      else if (toBool(state.getEigVal()[i] < 0.0))
	EigValFunc[i] = -1.0;
      else
	EigValFunc[i] = 0.0;
    }
  }

  writer << my_writer;
}

//! Produce a linear operator for this action
/*!
 * The operator acts on the entire lattice
 *
 * \param state_	 gauge field state  	 (Read)
 */
const LinearOperator<LatticeFermion>* 
Zolotarev4DFermActBj::linOp(Handle<const ConnectState> state_) const
{
  START_CODE("Zolotarev4DLinOp::create");

  const ZolotarevConnectState<LatticeFermion>& state = dynamic_cast<const ZolotarevConnectState<LatticeFermion>&>(*state_);

  if (state.getEigVec().size() != state.getEigVal().size())
    QDP_error_exit("Zolotarev4DLinOp: inconsistent sizes of eigenvectors and values");

  int NEigVal = state.getEigVal().size();

  /* The actual number of eigenvectors to project out.
     The highest of the valid low eigenmodes is not
     projected out. So we will put NEig = NEigVal - 1 */  
  int NEig;

  /* The number of residuals and poles */
  int numroot;

  /* The roots, i.e., the shifts in the partial fraction expansion */
  multi1d<Real> rootQ;

  /* The residuals in the partial fraction expansion */
  multi1d<Real> resP;

  /* This will be our alpha(0) which can be 0 depending on type */
  /* an even- or oddness of RatPolyDeg*/
  Real coeffP; 

  /* Array of values of the sign function evaluated on the eigenvectors of H */
  multi1d<Real> EigValFunc(NEigVal);

  // Common initialization
  init(numroot, coeffP, resP, rootQ, NEig, EigValFunc, state);

  
  /* Finally construct and pack the operator */
  /* This is the operator of the form (1/2)*[(1+mu) + (1-mu)*gamma_5*eps] */
  
  return new lovlapms(*Mact, state_, m_q,
		      numroot, coeffP, resP, rootQ, 
		      NEig, EigValFunc, state.getEigVec(),
		      MaxCGinner, RsdCGinner);
  
  END_CODE("Zolotarev4DLinOp::create");
}

//! Produce a linear operator for this action
/*!
 * The operator acts on the entire lattice
 *
 * \param state_	 gauge field state  	 (Read)
 */
const LinearOperator<LatticeFermion>* 
Zolotarev4DFermActBj::lMdagM(Handle<const ConnectState> state_, const Chirality& ichiral) const
{

  // If chirality is none, return traditional MdagM
  if ( ichiral == CH_NONE ) {
    cout << "Returning Normal MdagM" << endl << flush;
    return lMdagM(state_);
  }
  else { 
    const ZolotarevConnectState<LatticeFermion>& state = dynamic_cast<const ZolotarevConnectState<LatticeFermion>&>(*state_);
    
    if (state.getEigVec().size() != state.getEigVal().size())
      QDP_error_exit("Zolotarev4DLinOp: inconsistent sizes of eigenvectors and values");

    int NEigVal = state.getEigVal().size();

    /* The actual number of eigenvectors to project out.
       The highest of the valid low eigenmodes is not
       projected out. So we will put NEig = NEigVal - 1 */  
    int NEig;
    
    /* The number of residuals and poles */
    int numroot;

    /* The roots, i.e., the shifts in the partial fraction expansion */
    multi1d<Real> rootQ;
    
    /* The residuals in the partial fraction expansion */
    multi1d<Real> resP;
    
    /* This will be our alpha(0) which can be 0 depending on type */
    /* an even- or oddness of RatPolyDeg*/
    Real coeffP; 
    
    /* Array of values of the sign function evaluated on the eigenvectors of H */
    multi1d<Real> EigValFunc(NEigVal);
    
  // Common initialization
    init(numroot, coeffP, resP, rootQ, NEig, EigValFunc, state);

  
    /* Finally construct and pack the operator */
    /* This is the operator of the form (1/2)*[(1+mu) + (1-mu)*gamma_5*eps] */
    cout << "Returning lovddag" << endl << flush;
    return new lovddag(*Mact, state_, m_q,
		       numroot, coeffP, resP, rootQ, 
		       NEig, EigValFunc, state.getEigVec(),
		       MaxCGinner, RsdCGinner, ichiral);
  
  }
  END_CODE("Zolotarev4DlMdagM");
}

//! Produce a conventional MdagM operator for this action
/*!
 * This is the operator which you can always use
 * The operator acts on the entire lattice
 *
 * \param state_	 gauge field state  	 (Read)
 */
const LinearOperator<LatticeFermion>* 
Zolotarev4DFermActBj::lMdagM(Handle<const ConnectState> state_) const
{
  // linOp news the linear operator and gives back pointer, 
  // We call lmdagm with this pointer.
  // lmdagm is the only owner
  // No need to grab linOp with handle at this stage.
  return new lmdagm<LatticeFermion>( linOp(state_) );
}


const ZolotarevConnectState<LatticeFermion>*
Zolotarev4DFermActBj::createState(const multi1d<LatticeColorMatrix>& u_,
				  const Real& approxMin_) const 
{
  /* if ( ! ( approxMin_ > Real(0) )) { 
    ostringstream error_str;
    error_str << "Zolotarev4DFermActBj: approxMin_ has to be positive" << endl;
    throw error_str.str();
  }
  */

  // First put in the BC
  multi1d<LatticeColorMatrix> u_tmp = u_;
  getFermBC().modifyU(u_tmp);

  Real approxMax = Real(2)*Real(Nd);
  return new ZolotarevConnectState<LatticeFermion>(u_tmp, approxMin_, approxMax);
}


const ZolotarevConnectState<LatticeFermion>*
Zolotarev4DFermActBj::createState(const multi1d<LatticeColorMatrix>& u_,
				  const Real& approxMin_,
				  const Real& approxMax_) const
{
  ostringstream error_str;
  
  /*
  if ( !(approxMin_ > 0 )) { 
    error_str << "Zolotarev4DFermActBj: approxMin_ has to be positive" << endl;
    throw error_str.str();
  }

  if ( ! (approxMax_ > approxMin_) ) { 
    error_str << "Zolotarev4DFermActBj: approxMax_ has to be larger than approxMin_" << endl;
    throw error_str.str();
  }
  */

  // First put in the BC
  multi1d<LatticeColorMatrix> u_tmp = u_;
  getFermBC().modifyU(u_tmp);

  return new ZolotarevConnectState<LatticeFermion>(u_tmp, approxMin_, approxMax_);
}



const ZolotarevConnectState<LatticeFermion>*
Zolotarev4DFermActBj::createState(const multi1d<LatticeColorMatrix>& u_,
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

  return new ZolotarevConnectState<LatticeFermion>(u_tmp, lambda_lo_, evecs_lo_, lambda_hi_, approxMin, approxMax);
}

