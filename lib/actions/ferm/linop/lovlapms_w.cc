// $Id: lovlapms_w.cc,v 1.9 2003-12-17 13:22:59 bjoo Exp $
/*! \file
 *  \brief Overlap-pole operator
 */
#include <math.h>
#include "chromabase.h"
#include "actions/ferm/linop/lovlapms_w.h"


using namespace QDP;

//! Apply the operator onto a source vector
/*! \ingroup linop
 *
 * The operator acts on the entire lattice
 *
 * \param psi 	  Pseudofermion field     	       (Read)
 * \param isign   Flag ( PLUS | MINUS )   	       (Read)
 */
void lovlapms::operator() (LatticeFermion& chi, const LatticeFermion& psi, 
			   enum PlusMinus isign) const
{

  // Debugging section:
  QDPIO::cout << "m_q " << m_q << endl;
  QDPIO::cout << "numroot " << numroot << endl;
  QDPIO::cout << "constP " << constP <<endl;

  QDPIO::cout << "rootQ.size() " << rootQ.size() << endl;
  for(int qcount=0; qcount < rootQ.size(); qcount++) {
    QDPIO::cout << "rootQ["<<qcount<<"] " << rootQ[qcount] << endl;
  }
  
  QDPIO::cout << "resP.size() " << resP.size() << endl;
  for(int qcount=0; qcount < resP.size(); qcount++) {
    QDPIO::cout << "resP[" << qcount <<"] " << resP[qcount] << endl;
  }

  QDPIO::cout << "NEig " << NEig << endl;  
  QDPIO::cout << "EigValFunc.size() " << EigValFunc.size() << endl;
  for(int qcount=0; qcount < EigValFunc.size(); qcount++) {
    QDPIO::cout << "EigValFunc[" << qcount <<"] " << EigValFunc[qcount] << endl;
  }

  QDPIO::cout << "MaxCG " << MaxCG << endl;
  QDPIO::cout << "RsdCG " << RsdCG << endl;

  LatticeFermion tmp1, tmp2;
  int k, n_count;

  // Gamma_5 
  int G5 = Ns*Ns - 1;

  // Mass for shifted system
  Real mass = Real(1 + m_q) / Real(1 - m_q);
  

  switch (isign)
  {
  case PLUS:
    /*  Non-Dagger: psi is source and tmp1 */
    /*  chi  :=  gamma_5 * (gamma_5 * mass + eps(H)) * Psi  */
    tmp1 = psi;
    break;

  case MINUS:
    /*  Dagger: apply gamma_5 to source psi to make tmp1 */
    /*  chi  :=  (mass + eps(H) * gamma_5) * Psi  */
    tmp1 = Gamma(G5) * psi;
    break;

  default:
    QDP_error_exit("unknown isign value", isign);
  }


  chi = zero;

  /* Project out eigenvectors of source if desired */
  /* chi  +=  func(lambda) * EigVec * <EigVec, psi>  */
  /* Usually "func(.)" is eps(.); it is precomputed in EigValFunc. */

  QDPIO::cout << "NEig = " << NEig << endl;

  if (NEig > 0)
  {
    Complex cconsts;

    for(int i = 0; i < NEig; ++i)
    {
      cconsts = innerProduct(EigVec[i], tmp1);
      tmp1 -= EigVec[i] * cconsts;

      cconsts *= EigValFunc[i];
      chi += EigVec[i] * cconsts;
    }
  }

  /* tmp1 <- H * Projected psi */
  /*      <- gamma_5 * M * psi */

  (*M)(tmp2, tmp1, PLUS);
  tmp1 = Gamma(G5) * tmp2;
  
 

  Double c = norm2(tmp1);

  /* If exactly 0 norm, then solution must be 0 (for pos. def. operator) */
  if (toBool(c == 0.0))
  {
    chi = zero;

    return;
  }


  QDP_info("Doing multi-mass solve\n");
  // *******************************************************************
  // Solve  (MdagM + rootQ_n) chi_n = H * tmp1
  LatticeFermion Ap;
  LatticeFermion r;
  multi1d<LatticeFermion> p(numroot);

  Real a;              // Alpha for unshifted (isz) system
  Real as;             // alpha for current shifted system
  Real b;              // Beta for unshifted system
  Real bp;             // Beta previous for unshifted system

  Real ztmp;
  Double cp;
  Double d;
  Real z0;
  Real z1;   
  multi1d<Real> bs(numroot);  // beta for shifted system
  multi2d<Real> z(2,numroot); // zeta for shifted system
  Double chi_sq_new;
  Double chi_sq_diff;
  multi1d<bool> convsP(numroot);  // convergence mask for shifted system
  bool convP;                     // overall convergence mask
  
  int iz;
  int s;

  
  multi1d<LatticeFermion> x(numroot);

  
  Real rsdcg_sq = RsdCG * RsdCG;
  Real rsd_sq = c * rsdcg_sq;

  x = zero;

  // By default (could change), rootQ(isz) is considered the smallest shift 
  int isz = numroot-1;
  

  // chi[0] := mass*psi + c0*H*tmp1 + Eigvecs; 
  if (isign == PLUS)
  {
    //  chi  :=  gamma_5 * (gamma_5 * mass + eps(H)) * Psi 
    // Final mult by gamma_5 is at end 
    chi += Gamma(G5) * psi * mass;
  }
  else
  {
    // chi  :=  (mass + eps(H) * gamma_5) . Psi  
    chi += psi * mass;
  }

  chi += tmp1 * constP;


  // r[0] := p[0] := tmp1 
  r = tmp1;

  for(s = 0; s < numroot; ++s){
    p[s] = tmp1;
  }

  convsP = false;
  convP = false;

  iz = 1;  // z_index  z[ iz , s ] holds zeta(s)
           //          z[ 1-iz, s ] holds zeta_minus(s)

  z = 1;   // This fills both zeta and zeta_minus

  a = 0;   // Alpha for unshifted
  b = 1;   // beta for unshifted 


  for(k = 0; k <= MaxCG && ! convP ; ++k) {

    //  b[k] := | r[k] |**2 / < p[k], Ap[k] > ; 
    //  First compute  d  =  < p, A.p >  
    //  Ap = A . p 


    // This bit computes 
    // Ap = [  M^dag M + rootQ(isz)  ] p_isz
    

    (*MdagM)(Ap, p[isz], PLUS);
    Ap += p[isz] * rootQ[isz];

#if 000000
    // Project out eigenvectors
    if (k % 2 == 0)
      GramSchm (Ap, 1, EigVec, NEig);
#endif

    //  d =  < p, A.p >
    d = innerProductReal(p[isz], Ap);                       // 2 Nc Ns  flops 
    
    bp = b;                        // Store previous unshifted beta
    b = -Real(c/d);                // New unshifted beta

    // Compute the shifted bs and z 
    bs[isz] = b;

    // iz now points to previous z's
    iz = 1 - iz;
    
  
    for(s = 0; s < numroot; s++)
    {

      // Do this to avoid mitsmp compiler bug !!
      z0 = z[1-iz][s];  // The current z for the system under considerstion
      z1 = z[iz][s] ;   // The previous z for the system under consideration
  
      // We only compute beta and z factors if
      //   i) The system is not yet converged
      //   ii) The system is shifted
      if (s != isz &&  !convsP[s]) {


	// We write the new z-s in the place of the previous ones. 
	// Our ones will become previous next time around
	z[iz][s]  = z0*z1*bp ; 
	z[iz][s] /=  b*a*(z1 - z0) + z1*bp*(1 - (rootQ[s] - rootQ[isz])*b);
	bs[s] = b * z[iz][s] / z0;

      }

     
    }
    


    // r[k+1] += b[k] A . p[k] ; 
    r += b * Ap;	        // 2 Nc Ns  flops 


#if 000000
    // Project out eigenvectors 
    if (k % 2 == 0)
      GramSchm (r, 1, EigVec, NEig);
#endif
    
    //  chi[k+1] -= b[k] p[k] ; 

    // This is implemented as a sum. I don't know if this is kosher or not:
    //
    // Basically sgn(H) tmp1 
    // = sum_shifts resP(shift) x_shift 
    // = sum_shifts resP(shift) x_shift(iter)
    // = sum_shifts resP(shift) ( x_i(shift) - beta_
    ztmp = resP[isz] * b;
    tmp1 = p[isz] * ztmp;	// 2 Nc Ns  flops 

    
    for(s = 0; s < numroot; ++s) {
      if(s != isz  &&  !convsP[s]) {

	ztmp = bs[s] * resP[s];
	tmp1 += p[s] * ztmp;	// 2 Nc Ns  flops
      }
    }

    chi -= tmp1;                   // 2 Nc Ns  flops

    //  cp  =  | r[k] |**2 
    cp = c;

    //  c  =  | r[k] |**2 
    c = norm2(r);	                // 2 Nc Ns  flops

    //  a[k+1] := |r[k]|**2 / |r[k-1]|**2 ; 
    a = Real(c/cp);

    //  p[k+1] := r[k+1] + a[k+1] p[k];
    //  Compute the shifted as 
    //  ps[k+1] := zs[k+1] r[k+1] + a[k+1] ps[k]; 
 
    for(s = 0; s < numroot; ++s)
    {
      if (! convsP[s]) {

	if (s == isz) {
	  // Compute the actual solutions just for debugging
	  x[s] -= b * p[s];
    
	  p[s] *= a;	        // Nc Ns  flops 
	  p[s] += r;		// Nc Ns  flops 
	}
	else {
	  x[s] -= bs[s]*p[s];
	  
	  as = a * z[iz][s]*bs[s] / (z[1-iz][s]*b);
	  
	  p[s] *= as;	        // Nc Ns  flops 
	  p[s] += r * z[iz][s];	// Nc Ns  flops 
	  
	}
      }
    }

#if 0
    // Project out eigenvectors 
    if (k % 10 == 0)
      GramSchm (p, numroot, EigVec, NEig);
#endif
    
    // IF |r[k]| <= RsdCG |chi| THEN RETURN;
    convP = true;                          // Assume convergence and prove
                                           // otherwise

    cout << endl;
    cout << "Iter " << k << endl;
    for(s = 0; s < numroot; ++s) {
      if (! convsP[s]) {
	
        ztmp = Real(c) * z[iz][s]*z[iz][s];	
	cout << "Rsd(" << s << ") = " << ztmp << endl;
	bool btmp = toBool(ztmp < rsd_sq);
	convP = convP & btmp;
	convsP[s] = btmp;
      }
    }
    

  
    // eps(H) psi = ( c_0 + \sum_i { res_i / (z^2 + rootQ_i) } ) * H * psi 
    // Norm of diff of new and old soln unchanged by a potential gamma_5 
    // Check for convergence of chi 
    // Only converge if chi is converged. If vectors converge first, then error
    
    
    if (k > 0 &&  !convP) {
      chi_sq_new = rsdcg_sq * norm2(chi);
      chi_sq_diff = norm2(tmp1);      // the diff of old and new soln

      bool btmp = toBool(chi_sq_diff < chi_sq_new);

#if 0
      QDPIO::cout << "Lovlapms: k =" << k <<" chisq_new=" << chi_sq_new << " chisq_diff=" << Real(chi_sq_diff) << endl;
#endif

      if (! btmp && convP)
	QDP_error_exit("vectors converged but not final chi");

      /* convP = convP & btmp;  */
      convP = btmp;
    }
    
  }


  n_count = k;

  if (isign == PLUS)
  {
    // chi  :=  gamma_5 * (gamma_5 * mass + eps(H)) * Psi  
    tmp1 = Gamma(G5) * chi;
    chi = tmp1;
  }

  // Rescale to the correct normalization 
  chi *= 0.5 * (1 - m_q);


}

