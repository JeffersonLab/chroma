// $Id: zolotarev4d_linop_w.cc,v 1.9 2003-12-09 10:33:42 bjoo Exp $
/*! \file
 *  \brief 4D Zolotarev operator
 */

#include "chromabase.h"
#include "actions/ferm/linop/lmdagm_w.h"
#include "actions/ferm/linop/zolotarev4d_linop_w.h"
#include <zolotarev.h>

// #include "primitives.h"
//#include "common_declarations.h"

//! Internal creation routine
/*!
 */
void 
Zolotarev4DLinOp::init()
{
  START_CODE("Zolotarev4DLinOp::create");

  if (state.getEigVec().size() != state.getEigVal().size())
    QDP_error_exit("Zolotarev4DLinOp: inconsistent sizes of eigenvectors and values");

  NEigVal = state.getEigVal().size();

  /* The residual for the solutions of the multi-shift linear system */
  RsdCGinner = 1.0e-7;  // Hardwired the accuracy

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

  /* This will be our alpha(0) which can be 0 depending on type */
  /* an even- or oddness of RatPolyDeg*/
  Real coeffP; 
  /* The number of residuals and poles */
  int numroot;
  /* The roots, i.e., the shifts in the partial fraction expansion */
  multi1d<Real> rootQ(numroot);
  /* The residuals in the partial fraction expansion */
  multi1d<Real> resP(numroot);

  /* Array of values of the sign function evaluated on the eigenvectors of H */
  multi1d<Real> EigValFunc(NEigVal);

  /* The actual number of eigenvectors to project out.
     The highest of the valid low eigenmodes is not
     projected out. So we will put NEig = NEigVal - 1 */  
  int NEig;

  /* Auxiliary variables and counters */
  Real t;
  Real c;

    
  /* Accuracy to compute Q^(-1) */
  RsdCGinner;
  
  /* Hermitian 4D overlap operator 1/2 ( 1 + m_q + (1 - m_q) gamma5 * sgn(H)) 
     using a partial fraction expansion of the optimal rational function
     approximation to sgn. Here, H = 1/2 * gamma_5 * (1/kappa - D'). 
     The coefficients are computed by Zolotarev's formula. */

  /* The square M^dagger*M of the Wilson Dirac operators, used for
     solving the multi-shift linear system */
  /* H^2 = M^dag . M */
  LinearOperatorProxy<LatticeFermion>* MdagM = new lmdagm(M); 
  
  /* The operator gamma_5 * M with the M constructed here has its eigenvalues
     in the range m/(m + Nd) <= |gamma_5 * M| <= (m + 2*Nd)/(m + Nd) (in the 
     free case) where here m is arbitrary.
     So if we multiply M by a factor scale_fac = (m + Nd)/(m + 2*Nd) we have
     |gamma_5 * M| <= 1. */
  if(NEigVal == 0 ) 
  {
    scale_fac = (Real(Nd) - OverMass)/(Real(2)*Real(Nd) - OverMass);
    /* Choose some arbitrary value for the lower bound of the
       approximation interval <=> dangerous! */
    eps = Real(0.00089);
    NEig = 0;
  }
  else
  {
    scale_fac = Real(1.) / PolyRangUp ;
    eps = PolyRangLow / PolyRangUp ;
    /* The highest of the valid low eigenmodes is not
       projected out. */
    NEig = NEigVal - 1;
#if 0   
    push(nml_out,"Number of Wilson eigenmodes to be projected out");
    Write(nml_out, NEigVal);
    Write(nml_out, NEig);
    pop(nml_out);
#endif
  }
 

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
  maxerr = rdata -> Delta;
#if 0
  push(nml_out,"ConsZolotarev4DLinOp");
  Write(nml_out, eps);
  Write(nml_out, RatPolyDeg);
  Write(nml_out, type);
  Write(nml_out, maxerr);
  pop(nml_out);
#endif

  /* Allocate the roots and residua */
  numroot = rdata -> dd;
      
  /* Fill in alpha[0] = alpha[da] if it is not zero*/
  coeffP = 0;
  coeffP = rdata -> alpha[rdata -> da - 1];
  /* The coefficients from the partial fraction.
     Here, we write them out for the sake of bookkeeping. */
  resP = 0;
  rootQ = 0;
  for(int n=0; n < numroot; ++n)
  {
    resP[n] = rdata -> alpha[n];
    rootQ[n] = rdata -> ap[n];
    rootQ[n] = -rootQ[n];
  }
#if 0
  push(nml_out,"Original partial fraction Zolotarev coeff");
  Write(nml_out, scale_fac);
  Write(nml_out, coeffP);
  Write(nml_out, resP);
  Write(nml_out, rootQ);
  pop(nml_out);
#endif

  /* Now fill in the coefficients for real, i.e., taking the rescaling
     into account */
  /* Fill in alpha[0] = alpha[da] if it is not zero*/
  coeffP = rdata -> alpha[rdata -> da - 1] * scale_fac;
  /* Fill in the coefficients for the roots and the residua */
  /* Make sure that the smallest shift is in the last value rootQ(numroot-1)*/
  resP = 0;
  rootQ = 0;
  t = Real(1) / (scale_fac * scale_fac);
  for(int n=0; n < numroot; ++n)
  {
    resP[n] = rdata -> alpha[n] / scale_fac;
    rootQ[n] = rdata -> ap[n];
    rootQ[n] = -(rootQ[n] * t);
  }

#if 0
  /* Write them out into the namelist */
  push(nml_out,"Rescaled partial fraction Zolotarev coeff");
  Write(nml_out, scale_fac);
  Write(nml_out, coeffP);
  Write(nml_out, resP);
  Write(nml_out, rootQ);
  pop(nml_out);
#endif

  QDP_info("ZOLOTAREV_4d: Mass= %g n= %d scale= %g coeff= %g  Nwils= %d  m_q= %g  Rsd= %g",
	 toFloat(OverMass),RatPolyDeg,toFloat(scale_fac),toFloat(coeffP),NEigVal,
	   toReal(m_q),toReal(RsdCGinner));
  QDP_info("Auxiliary fermion action: OverAuxAct = %d",OverAuxAct);
  QDP_info("Approximation on [-1,-eps] U [eps,1] with eps = %g",toReal(eps));
  maxerr = rdata -> Delta;
  QDP_info("Maximum error |R(x) - sgn(x)| <= Delta = %g",maxerr);
  if(type == 0) {QDP_info("Approximation type %d with R(0) = 0",type);}
  else {QDP_info("Approximation type %d with R(0) =  infinity",type);}
  
  /* We will also compute the 'function' of the eigenvalues */
  /* for the Wilson vectors to be projected out. */
  if (NEig > 0)
  {
      
    for(int i = 0; i < NEigVal; i++)
    {
      if (toBool(EigVal[i] > 0.0))
	EigValFunc[i] = 1.0;
      else if (toBool(EigVal[i] < 0.0))
	EigValFunc[i] = -1.0;
      else
	EigValFunc[i] = 0.0;
    }
  }


  /* Finally construct and pack the operator */
  /* This is the operator of the form (1/2)*[(1+mu) + (1-mu)*gamma_5*eps] */
  LinearOperator<LatticeFermion>* A = new lovlapms(*MdagM, M, m_q,
						   numroot, coeffP, resP, rootQ, 
						   EigVec, EigValFunc, NEig, 
						   MaxCG, RsdCGinner);
  
  END_CODE("Zolotarev4DLinOp::create");
}



//! Destructor
/*!
 * \ingroup linop
 */
void Zolotarev4DLinOp::~Zolotarev()
{
//  delete MdagM;
//  delete A;
}



//! Apply the operator onto a source vector
/*! \ingroup linop
 *
 * The operator acts on the entire lattice
 *
 * \param chi 	  Pseudofermion field     	       (Write)
 * \param psi 	  Pseudofermion field     	       (Read)
 * \param isign   Flag ( PLUS | MINUS )   	       (Read)
 */
void Zolotarev4DLinOp::operator() (LatticeFermion& chi, const LatticeFermion& psi, 
				   enum PlusMinus isign) const
{
  START_CODE("Zolotarev4DLinOp");

  // Call underlying lovlapms operator
  (*A)(chi, psi, isign);

  END_CODE("Zolotarev4DLinOp");
}

