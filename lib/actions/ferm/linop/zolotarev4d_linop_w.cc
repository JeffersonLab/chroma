// $Id: zolotarev4d_linop_w.cc,v 1.1 2003-04-09 05:57:16 edwards Exp $
/*! \file
 *  \brief 4D Zolotarev operator
 */

#include "chromabase.h"
#include "actions/construct/conslinop_w.h"
#include "actions/linop/lmdagm_w.h"
#include "actions/linop/zolotarev4d_w.h"
#include "actions/linop/zolotarev.h"
#include "primitives.h"
#include "common_declarations.h"
#include "common_io.h"

//! Creation routine
/*!
 * \ingroup linop
 *
 * \param _u 	  gauge field  	       (Read)
 * \param _m_q    quark mass   	       (Read)
 */

void Zolotarev4D::create(const multi1d<LatticeColorMatrix>& _u, const Real& _m_q)
{
  START_CODE("Zolotarev4D::create");

  u = _u;
  m_q = _m_q;

  RsdCGinner = 1.0e-7;  // Hardwired the accuracy

  NEigVal = 0;

  make();
};

//! Creation routine
/*!
 * \ingroup linop
 *
 * \param _u 	  gauge field  	       (Read)
 * \param _m_q    quark mass   	       (Read)
 * \param _EigVec eigenvectors 	       (Read)
 * \param _EigVal eigenvalues 	       (Read)
 */

void Zolotarev4D::create(const multi1d<LatticeColorMatrix>& _u, const Real& _m_q,
			 const multi1d<LatticeFermion>& _EigVec, const multi1d<Real>& _EigVal)
{
  START_CODE("Zolotarev4D::create");

  u = _u;
  m_q = _m_q;

  RsdCGinner = 1.0e-7;  // Hardwired the accuracy

  NEigVal = _EigVal.size();
  EigVec = _EigVec;
  EigVal = _EigVal;

  if (EigVec.size() != EigVal.size())
    QDP_error_exit("Zolotarev4D::create: inconsistent sizes of eigenvectors and values");

  make();
};

//! Internal creation routine
/*!
 * \ingroup linop
 */

void Zolotarev4D::make()
{
  START_CODE("Zolotarev4D::create");

  /* The residual for the solutions of the multi-shift linear system */
  Real RsdCG;

  /* The kappa value of the Wilson Dirac operator corresponding to the
     negative mass*/
  Real Kappa_t;
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
  int n;
  int i;
  Real t;
  Real c;

  /* More local variables here perhaps ? */
  START_CODE("subroutine");;

    
  /* Accuracy to compute Q^(-1) */
  RsdCG = RsdCGinner;
  
  /* Hermitian 4D overlap operator 1/2 ( 1 + m_q + (1 - m_q) gamma5 * sgn(H)) 
     using a partial fraction expansion of the optimal rational function
     approximation to sgn. Here, H = 1/2 * gamma_5 * (1/kappa - D'). 
     The coefficients are computed by Zolotarev's formula. */

  /* Basic Wilson operator for H. Kappa_t is the kappa corresponding
     to OverMass, Kappa_t = 1/2 / (Nd - OverMass). For GW OverMass >
     0, corresponding to a negative shift, i.e. the negative sign of 
     OverMass is taken explicitely care of. Normally, one has 
     kappa = 0.5 / (Nd + mass). */
  Kappa_t = TO_REAL(0.5) / (TO_REAL(Nd) - OverMass);

  /* Construct the auxilliary (Wilson) operator for use within Neuberger's
     form. Here, M = (1 - kappa_t * D') */

  /* Basic Wilson operator for H */
  if (OverAuxAct == PLANAR_WILSON)
    Kappa_t = - OverMass;
  else
    Kappa_t = TO_REAL(0.5) / (TO_REAL(Nd) - OverMass);

  /* A linear operator containing the non-hermitian Wilson Dirac
     operator with negative mass */
  LinearOperator* M = ConsLinOp (u, Kappa_t, OverAuxAct);
  QDP_info("Auxiliary fermion action: OverAuxAct = %d",OverAuxAct);
  
  /* The square M^dagger*M of the Wilson Dirac operators, used for
     solving the multi-shift linear system */
  /* H^2 = M^dag . M */
  LinearOperator* MdagM = new lmdagm(*M); 
  
  /* The operator gamma_5 * M with the M constructed here has its eigenvalues
     in the range m/(m + Nd) <= |gamma_5 * M| <= (m + 2*Nd)/(m + Nd) (in the 
     free case) where here m is arbitrary.
     So if we multiply M by a factor scale_fac = (m + Nd)/(m + 2*Nd) we have
     |gamma_5 * M| <= 1. */
  if(NEigVal == 0 ) 
  {
    scale_fac = (TO_REAL(Nd) - OverMass)/(TO_REAL(2)*TO_REAL(Nd) - OverMass);
    /* Choose some arbitrary value for the lower bound of the
       approximation interval <=> dangerous! */
    eps = TO_REAL(0.00089);
    NEig = 0;
  }
  else
  {
    scale_fac = TO_REAL(1.) / PolyRangUp ;
    eps = PolyRangLow / PolyRangUp ;
    /* The highest of the valid low eigenmodes is not
       projected out. */
    NEig = NEigVal - 1;
    push(nml_out,"Number of Wilson eigenmodes to be projected out");
    Write(nml_out, NEigVal);
    Write(nml_out, NEig);
    pop(nml_out);
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
  push(nml_out,"ConsZolotarev4D");
  Write(nml_out, eps);
  Write(nml_out, RatPolyDeg);
  Write(nml_out, type);
  Write(nml_out, maxerr);
  pop(nml_out);

  /* Allocate the roots and residua */
  numroot = rdata -> dd;
      
  /* Fill in alpha[0] = alpha[da] if it is not zero*/
  coeffP = 0;
  coeffP = rdata -> alpha[rdata -> da - 1];
  /* The coefficients from the partial fraction.
     Here, we write them out for the sake of bookkeeping. */
  resP = 0;
  rootQ = 0;
  for(n=0; n < numroot; ++n)
  {
    resP[n] = rdata -> alpha[n];
    rootQ[n] = rdata -> ap[n];
    rootQ[n] = -rootQ[n];
  }
  push(nml_out,"Original partial fraction Zolotarev coeff");
  Write(nml_out, scale_fac);
  Write(nml_out, coeffP);
  Write(nml_out, resP);
  Write(nml_out, rootQ);
  pop(nml_out);

  /* Now fill in the coefficients for real, i.e., taking the rescaling
     into account */
  /* Fill in alpha[0] = alpha[da] if it is not zero*/
  coeffP = rdata -> alpha[rdata -> da - 1] * scale_fac;
  /* Fill in the coefficients for the roots and the residua */
  /* Make sure that the smallest shift is in the last value rootQ(numroot-1)*/
  resP = 0;
  rootQ = 0;
  t = TO_REAL(1) / (scale_fac * scale_fac);
  for(n=0; n < numroot; ++n)
  {
    resP[n] = rdata -> alpha[n] / scale_fac;
    rootQ[n] = rdata -> ap[n];
    rootQ[n] = -(rootQ[n] * t);
  }

  /* Write them out into the namelist */
  push(nml_out,"Rescaled partial fraction Zolotarev coeff");
  Write(nml_out, scale_fac);
  Write(nml_out, coeffP);
  Write(nml_out, resP);
  Write(nml_out, rootQ);
  pop(nml_out);

  QDP_info("ZOLOTAREV_4d: Mass= %g n= %d scale= %g coeff= %g  Nwils= %d  m_q= %g  Rsd= %g",
	 toFloat(OverMass),RatPolyDeg,toFloat(scale_fac),toFloat(coeffP),NEigVal,
	   toReal(m_q),toReal(RsdCG));
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
      
    for(i = 0; i < NEigVal; i++)
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
  LinearOperator* A = new lovlapms(*MdagM, M, m_q,
				   numroot, coeffP, resP, rootQ, 
				   EigVec, EigValFunc, NEig, RsdCG);
  
  END_CODE("Zolotarev4D::create");
}



//! Destructor
/*!
 * \ingroup linop
 */
void Zolotarev4D::~Zolotarev()
{
  delete MdagM;
  delete A;
}



//! Apply the operator onto a source vector
/*! \ingroup linop
 *
 * The operator acts on the entire lattice
 *
 * \param psi 	  Pseudofermion field     	       (Read)
 * \param isign   Flag ( PLUS | MINUS )   	       (Read)
 */
LatticeFermion Zolotarev4D::operator() (const LatticeFermion& psi, enum LinOpSign isign) const
{
  LatticeFermion chi;

  START_CODE("Zolotarev4D");

  // Call underlying lovlapms operator
  chi = (*A)(psi, isign);

  END_CODE("Zolotarev4D");

  return chi;
}

