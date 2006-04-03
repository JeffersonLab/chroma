// $Id: lovlap_double_pass_w.cc,v 3.0 2006-04-03 04:58:50 edwards Exp $
/*! \file
 *  \brief Overlap-pole operator
 */
#include <math.h>
#include "chromabase.h"
#include "actions/ferm/linop/lovlap_double_pass_w.h"
#include "meas/eig/gramschm.h"


namespace Chroma 
{ 
//! Apply the GW operator onto a source vector
/*! \ingroup linop
 *
 * This routine applies the 4D GW operator onto a source
 * vector. The coeffiecients for the approximation get 
 * wired into the class by the constructor and should
 * come fromt fermion action.
 *
 * The operator applied is:
 *       D       =    (1/2)[  (1+m) + (1-m)gamma_5 sgn(H_w) ] psi
 * or    D^{dag} =    (1/2)[  (1+m) + (1-m) sgn(H_w) gamma_5 psi
 * 
 * 
 * \param chi     result vector                              (Write)  
 * \param psi 	  source vector         	             (Read)
 * \param isign   Hermitian Conjugation Flag 
 *                ( PLUS = no dagger| MINUS = dagger )       (Read)
 */
void lovlap_double_pass::operator() (LatticeFermion& chi, 
				     const LatticeFermion& psi, 
				     enum PlusMinus isign) const
{
  operator()(chi, psi, isign, RsdCG);
}


//! Apply the GW operator onto a source vector
/*! \ingroup linop
 *
 * This routine applies the 4D GW operator onto a source
 * vector. The coeffiecients for the approximation get 
 * wired into the class by the constructor and should
 * come fromt fermion action.
 *
 * The operator applied is:
 *       D       =    (1/2)[  (1+m) + (1-m)gamma_5 sgn(H_w) ] psi
 * or    D^{dag} =    (1/2)[  (1+m) + (1-m) sgn(H_w) gamma_5 psi
 * 
 * 
 * \param chi     result vector                              (Write)  
 * \param psi 	  source vector         	             (Read)
 * \param isign   Hermitian Conjugation Flag 
 *                ( PLUS = no dagger| MINUS = dagger )       (Read)
 */
void lovlap_double_pass::operator() (LatticeFermion& chi, 
				     const LatticeFermion& psi, 
				     enum PlusMinus isign, Real epsilon) const
{
  START_CODE();

  LatticeFermion tmp1, tmp2;

  // Gamma_5 
  int G5 = Ns*Ns - 1;

  // Mass for shifted system
  Real mass = Real(1 + m_q) / Real(1 - m_q);
  

  switch (isign)
  {
  case PLUS:
    //  Non-Dagger: psi is source and tmp1 
    //  chi  :=  gamma_5 * (gamma_5 * mass + sgn(H)) * Psi  
    tmp1 = psi;
    break;

  case MINUS:
    // Dagger: apply gamma_5 to source psi to make tmp1 
    //  chi  :=  (mass + sgn(H) * gamma_5) * Psi  
    tmp1 = Gamma(G5) * psi;
    break;

  default:
    QDP_error_exit("unknown isign value", isign);
  }


  chi = zero;

  // Project out eigenvectors of source if desired 
  // chi  +=  func(lambda) * EigVec * <EigVec, psi>  
  // Usually "func(.)" is sgn(.); it is precomputed in EigValFunc. 
  // for all the eigenvalues
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

  // tmp1 <- H * Projected psi 
  //      <- gamma_5 * M * psi 
  (*M)(tmp2, tmp1, PLUS);
  tmp1 = Gamma(G5) * tmp2;
  
 
  // Effectively the multi mass solve is on tmp1 
  // *******************************************************************
  // Solve  (MdagM + rootQ_n) chi_n = H * tmp1

  Double c = norm2(tmp1);
  Double rsd_sq = norm2(psi) * epsilon*epsilon; // || psi ||^2 * epsilon^2

  // We are solving 2/(1-m) Q so rescale target residuum by ((1-m)/2)^2
  rsd_sq *=  (Real(1) - m_q)*(Real(1)-m_q)/Real(4);
  Double cp;
  Double d; // InnerProduct

  /* If exactly 0 norm, then solution must be 0 (for pos. def. operator) */
  if (toBool(c == 0))
  {
    chi = zero;
    END_CODE();
    return;
  }


  LatticeFermion Ap;
  LatticeFermion r;
  LatticeFermion p;

  multi1d<Double> a(MaxCG+1);              // Alpha for unshifted (isz) system
  multi1d<Double> b(MaxCG+1);              // Beta for unshifted system

  bool convP;
  int  iters_taken;

  // By default, rootQ(isz) is considered the smallest shift 
  int isz = numroot-1;             // isz identifies system with smalles shift
  

  // chi[0] := mass*psi + c0*H*tmp1 + Eigvecs; 
  if (isign == PLUS)
  {
    //  chi  :=  gamma_5 * (gamma_5 * mass + eps(H)) * Psi 
    // Final mult by gamma_5 is at end 
    tmp2 = Gamma(G5) * psi;

    // This will be an axpy
    chi += tmp2 * mass;
  }
  else
  {
    // chi  :=  (mass + eps(H) * gamma_5) . Psi  
    chi += psi * mass;
  }

  // Multiply in P(0) -- this may well be 0 for type 0 rational approximations
  chi += tmp1 * constP;


  // r[0] := p[0] := tmp1 
  r = tmp1;
  p = tmp1;


  cp = 0;
  convP = false;

  
  Double alpha_minus_one=1;
  b[0] = 0;         // b[0] -- do we ever use ?
  


  int k;
  multi1d<Double> c_iter(MaxCG+1);

  c_iter[0] = c;
  // Do the iterations proper 
  // first pass
  for(k = 0; k < MaxCG && ! convP ; ++k) {

    // Keep hold of the residuum
 
    (*MdagM)(Ap, p, PLUS);
    Ap += p * rootQ[isz];

    // Project out eigenvectors
    if (k % ReorthFreq == 0) {
      GramSchm(Ap, EigVec, NEig, all);
    }

    //  d =  < p, A.p >
    d = innerProductReal(p, Ap);                       // 2 Nc Ns  flops 
    

    // k+1 because a[0] corresponds to a_{-1}
    a[k] = c/d;
    r -= Real(a[k])*Ap;

    // Project out eigenvectors 
    if (k % ReorthFreq == 0) {
      GramSchm (r, EigVec, NEig, all);
    }
    
    cp = c;
    c = norm2(r);

    c_iter[k+1] = c;
    b[k+1] = c/cp;

    p = r + Real(b[k+1])*p;
 
    convP = toBool( c < rsd_sq );
  }

  int niters = k;

  // OK First pass done. I now need to compute gamma_j and c_j
  multi2d<Double> gamma(niters+1, numroot);
  for(int s = 0; s < numroot; s++) { 
    gamma[0][s] = 1;                     // Really gamma_0
  }
  

  multi1d<bool> convPs(numroot);
  multi1d<int>  convIter(numroot);
  convPs = false;

  for(int j=0; j < niters; j++) {
    for(int s =0; s < numroot; s++) {


      if( s != isz  && ! convPs[s]  ) { 
	Double a_minus;

	if( j == 0 ) { 
	  // alpha_minus_one -- no need to store. Only need alpha_0,...
	  // in second pass
	  a_minus = Double(1);
	}
	else {
	  a_minus = a[j-1];
	}

	Double ga = gamma[j][s];
	Double ga_minus;
	if( j == 0 ) { 
	  // gamma[-1][s] -- no need to store. Only need gamma[0][s]...
	  // in second pass
	  ga_minus =Double(1);
	}
	else{ 
	  ga_minus = gamma[j-1][s];
	}

	Double tmp_num = ga*ga_minus*a_minus;
	Double tmp_den = ga_minus*a_minus*(Double(1) + a[j]*(rootQ[s] - rootQ[isz]));
        tmp_den += a[j]*b[j]*(ga_minus - ga);

	gamma[j+1][s] = tmp_num/tmp_den;

	// If this system has converged at iter j+1, then dont update
	// gamma-s anymore and note the convergence point. Updating
	// ad infinitum causes underflow.
	if( toBool( gamma[j+1][s]*gamma[j+1][s]*c_iter[j+1] < rsd_sq ) ) {
	  convPs[s] = true;
	  convIter[s] = j+1;
	}

      }
      
	
    }

    
  } 

 
  multi1d<Double> sumC(niters+1);

  for(int j=0; j<=niters; j++) { 

    sumC[j] = 0;

    for(int m=0; m < niters-j; m++) { 
      
      Double qsum=resP[isz];
     
      for(int s =0; s < numroot; s++) { 

	if( s != isz ) { 
	  
	  // Only add gammas which are unconverged.
	  // Converged gamma-s are not updated beyond convergence
	  // which would cause problems (with underflows)
	  // gamma[convIters[s]] is last valid gamma
	  if( toBool( (j+m+1) <= convIter[s] ) ) {
	    qsum += resP[s]*gamma[m+j+1][s]*gamma[m+j][s]/gamma[j][s];
	  }
	  
	}
      }
      
          
      Double delta_m = Double(1);
      
      if( m > 0 ) { 
	for(int k=1; k <= m; k++) { 
	  delta_m *= b[j+k];
	}
      }
      
      
      sumC[j] += a[j+m]*delta_m*qsum;
    }
  }
    

  // Second pass Lanczos

  // Initialise r_0, p_0
  r=tmp1;
  p=tmp1;

  convP = false;

  for(k=0; k < niters && !convP ; k++) { 
    (*MdagM)(Ap, p, PLUS);
    Ap += p * rootQ[isz];

    // Project out eigenvectors
    if (k % ReorthFreq == 0) {
      GramSchm(Ap, EigVec, NEig, all);
    }

 
    // Update chi
    chi += Real(sumC[k])*r;

    // Update r
    r -= Real(a[k])*Ap;

    // Project out eigenvectors 
    if (k % ReorthFreq == 0) {
      GramSchm (r, EigVec, NEig, all);
    }

#if 0

    // This early termination criterion has to be abandoned when
    // using relaxed solevers
    Double chi_norm_new = epsilon*epsilon*norm2(chi);

    // Convergence criterion for total signum. Might be good enough
    // without running to full niters
    convP = toBool( c_iter[k+1]*sumC[k+1]*sumC[k+1] < chi_norm_new );
#endif
    p = r + Real(b[k+1])*p;
  }


  QDPIO::cout << "Overlap Inner Solve (lovlap_double_pass): " << k << " iterations " << endl;
  // End of MULTI SHIFTERY 

  // Now fix up the thing. Multiply in gamma5 if needed 
  // and then rescale to correct normalisation.

  if (isign == PLUS)
  {
    // chi  :=  gamma_5 * (gamma_5 * mass + eps(H)) * Psi  
    tmp1 = Gamma(G5) * chi;
    chi = tmp1;
  }

  // Rescale to the correct normalization 
  chi *= 0.5 * (1 - m_q);

  END_CODE();
}


}; // End Namespace Chroma


