// $Id: lovddag_double_pass_w.cc,v 3.0 2006-04-03 04:58:50 edwards Exp $
/*! \file
 *  \brief Overlap-pole operator
 */
#include <math.h>
#include "chromabase.h"
#include "actions/ferm/linop/lovddag_double_pass_w.h"
#include "meas/eig/gramschm.h"


namespace Chroma 
{ 
//! Apply the GW operator onto a source vector
/*! \ingroup linop
 *
 */

// We are evaluating 
//
//    init_fac*psi + (gamma_5 +/- 1) eps(H) psi_{+/-}
//
// Split psi  as  P_evec psi + (1 - P_evec) psi
// 
// and 
//
//    apply eps(H, P_evec psi) + eps(H, (1-P_evec)psi)
//
//     P_evec psi = Sum_i lambda_i | e_i >< e_i | psi>
// so  eps(H, P_evec psi ) = Sum_i eps(lambda_i) | e_i >< e_i | psi >
//
//   and ( 1 - P_evec ) psi = psi - P_evec_psi 
//                          = psi - Sum_i lambda_i | e_i >< e_i | psi >
//
// so  eps(H, (1-P_evec) psi )) = eps(H, Y) 
//                              = [ p_0 + sum_r p_r * (H^2 + q)^{-1} ] Y

void lovddag_double_pass::operator() (LatticeFermion& chi, 
				      const LatticeFermion& psi, 
				      enum PlusMinus isign) const
{
  operator()(chi,psi,isign,RsdCG);
}


//! Apply the GW operator onto a source vector
/*! \ingroup linop
 *
 */

// We are evaluating 
//
//    init_fac*psi + (gamma_5 +/- 1) eps(H) psi_{+/-}
//
// Split psi  as  P_evec psi + (1 - P_evec) psi
// 
// and 
//
//    apply eps(H, P_evec psi) + eps(H, (1-P_evec)psi)
//
//     P_evec psi = Sum_i lambda_i | e_i >< e_i | psi>
// so  eps(H, P_evec psi ) = Sum_i eps(lambda_i) | e_i >< e_i | psi >
//
//   and ( 1 - P_evec ) psi = psi - P_evec_psi 
//                          = psi - Sum_i lambda_i | e_i >< e_i | psi >
//
// so  eps(H, (1-P_evec) psi )) = eps(H, Y) 
//                              = [ p_0 + sum_r p_r * (H^2 + q)^{-1} ] Y

void lovddag_double_pass::operator() (LatticeFermion& chi, 
				      const LatticeFermion& psi, 
				      enum PlusMinus isign, 
				      Real epsilon) const
{
  START_CODE();

  LatticeFermion ltmp;
  LatticeFermion tmp1;
  Real ftmp;

  const int G5 = Ns * Ns - 1;
  int n;

  Real mass =  ( Real(1) + m_q ) / ( Real(1) - m_q );
  Real init_fac = (Real(1)/mass) + mass;

  LatticeFermion psi_proj = psi;

  // Accumulate the initial factor * psi this needs no projection
  // or any of the (gamma_5  + 1 ) tricks.
  chi = init_fac * psi;

 
  // Now we need to work out ( gamma_5 +/- 1 ) sgn(H) psi_{+/-}
  // Split psi into psi_proj and evecs
  //
  // evecs = sum lambda_i | e_i ><e_i | psi>

  // psi_proj = psi - sum lambda_i | e_i> <e_i | psi>
  //
  // Applying the ( gamma_5 +/-1 ) sgn(H) evecs
  // = (gamma_5 +/-1) sum sgn(lambda_i) |e_i><e_i | psi>
  // = (gamma_5 +/-1) EigValFunc[i] | e_i> <e_i|psi>
  // 
  // We accumulate this directly onto chi, and compute psi_proj 
  // in the process
  
  if (NEig > 0)
  {
    Complex cconsts;

    for(int i = 0; i < NEig; ++i)
    {
      cconsts = innerProduct(EigVec[i], psi_proj);
      
      psi_proj -= cconsts * EigVec[i];

   
      // Construct  tmp1 = (gamma_5 + ichiral) * EigVec 
      tmp1 = Gamma(G5)*EigVec[i];

      switch (ichiral) {
      case CH_PLUS:
	tmp1 += EigVec[i];
	break;
      case CH_MINUS:
	tmp1 -= EigVec[i];
	break;
      case CH_NONE:
	{
	  ostringstream error_message;
	  error_message << "Should not use lovddag_double_pass unless the chirality is either CH_PLUS or CH_MINUS" << endl;
	  throw error_message.str();
	}
	break;
      default:
	QDP_error_exit("ichiral is outside allowed range: %d\n", ichiral);
	break;
      }

      /* chi += "eps(lambda)" * (gamma_5 + ichiral) * EigVec  */
      cconsts *= EigValFunc[i];
      chi += cconsts * tmp1;

    }
    
  }
  
  // Chi Now holds the projected part + init_fac part
  // ie:
  // chi = init_fac*psi + (gamma_5 + ichiral) sum sgn(lambda) <e_i|psi> |e_i>

  // Now we have to apply (gamma_5 + 1 ) eps psi_proj
  //
  // We do this by doing the double pass algorithm
  // In this algorithm we do a multi_shift solve on 
  //  (H^2 + Q_i) x_i = H^psi_proj
  //
  // for all i, and then evaluate (gamma_5 + 1) ( constP H^psi_proj + 
  //                                             sum resP_i x_i
  //
  // The sum is carried out by the double pass
  //

  // tmp_1 will hold the source for the inversion
  // ie tmp_1 = H psi_proj
  //          = gamma_5 D_w psi_proj
  (*M)( ltmp, psi_proj, PLUS);
  tmp1 = Gamma(G5)*ltmp;
  Double c = norm2(tmp1);

  // If exactly 0 norm, then solution must be 0 (for pos. def. operator)
  // In which case we stop right now...
  if ( toBool(c == Real(0)) ) { 
    chi = zero;
    END_CODE();
    return;
  }



  // Now do the multi shiftery using the double pass
  // *******************************************************************
  // Solve  (MdagM + rootQ_n) chi_n = tmp1
  //
  // We put this part of the solution into ltmp (which we don't need
  // elsewhere

  Double rsd_sq = norm2(psi) * epsilon*epsilon; // || psi ||^2 * epsilon^2
  
  // Actually we are solving 4/( 1 - m_q^2 ) so need to rescale
  // target residue by (1-m_q^2)/4
  // (squared of course for rsd_sq)
  rsd_sq *= ( 1 - m_q*m_q)*(1 - m_q*m_q)/Real(16);

  Double cp;
  Double d; // InnerProduct
  
  LatticeFermion Ap;
  LatticeFermion r;
  LatticeFermion p;

  multi1d<Double> a(MaxCG+1);              // Alpha for unshifted (isz) system
  multi1d<Double> b(MaxCG+1);              // Beta for unshifted system

  bool convP;
  int  iters_taken;

  // By default, rootQ(isz) is considered the smallest shift 
  int isz = numroot-1;             // isz identifies system with smalles shift
  

  // Multiply in P(0) -- this may well be 0 for type 0 rational approximations
  ltmp = tmp1 * constP;


  // r[0] := p[0] := tmp1 
  r = tmp1;
  p = tmp1;
  

  cp = 0;
  convP = false;
 
  b[0] = 0;         
  int k;

  // Record || r_i ||^2 for each iteration
  multi1d<Double> c_iter(MaxCG+1);

  // We already know it for iter 0
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

  // Now work out the c_j constants for the second pass
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
    ltmp += Real(sumC[k])*r;

    // Update r
    r -= Real(a[k])*Ap;

    // Project out eigenvectors 
    if (k % ReorthFreq == 0) {
      GramSchm (r, EigVec, NEig, all);
    }

#if 0
    // I have to abandon this stopping criterion 
    // when using relaxed solver
    Double ltmp_norm_new = epsilon*epsilon*norm2(ltmp);

    // Convergence criterion for total signum. Might be good enough
    // without running to full niters
    convP = toBool( c_iter[k+1]*sumC[k+1]*sumC[k+1] < ltmp_norm_new );
#endif

    p = r + Real(b[k+1])*p;
  }

  // Now we have 
  // ltmp =  eps(H) psi_proj 

  // we have to add (gamma_5 + ichiral) eps(H) psi_proj
  // onto chi
  //
  tmp1 = Gamma(G5)*ltmp;
  switch (ichiral) {
  case CH_PLUS:
    tmp1 += ltmp;
    break;

  case CH_MINUS:
    tmp1 -= ltmp;
    break;

  default:
    QDP_error_exit("unsupported chirality", ichiral);
  }
  chi += tmp1;

  /* Rescale to the correct normalization */
  /*  chi <-  (1/4)*(1-m_q^2) * chi  */
  ftmp = Real(0.25) * (Real(1) - m_q*m_q);
  chi *= ftmp;	        /* 2 Nc Ns  flops */
  QDPIO::cout << "Overlap Inner Solve (lovddag_double_pass(" << ichiral << ")) = " << k << " iterations" << endl;

  END_CODE();
}

}; // End Namespace Chroma

