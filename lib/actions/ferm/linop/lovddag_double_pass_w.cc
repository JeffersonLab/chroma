// $Id: lovddag_double_pass_w.cc,v 1.2 2004-05-03 18:03:43 bjoo Exp $
/*! \file
 *  \brief Overlap-pole operator
 */
#include <math.h>
#include "chromabase.h"
#include "actions/ferm/linop/lovddag_double_pass_w.h"
#include "meas/eig/gramschm.h"

using namespace QDP;

//! Apply the GW operator onto a source vector
/*! \ingroup linop
 *
 */
void lovddag_double_pass::operator() (LatticeFermion& chi, const LatticeFermion& psi, 
			   enum PlusMinus isign) const
{
 

  LatticeFermion ltmp;
  LatticeFermion tmp1;
  Real ftmp;

  const int G5 = Ns * Ns - 1;
  int n;



  Real mass =  ( Real(1) + m_q ) / ( Real(1) - m_q );
  Real init_fac = (Real(1)/mass) + mass;

  LatticeFermion psi_proj = psi;

  chi = init_fac * psi;

 
  // Project out eigenvectors of source if desired 
  // chi  +=  eps(lambda) * (gamma_5 + ichiral) * EigVec  
  // Usually "func(.)" is eps(.); it is precomputed in EigValFunc.
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

  /* First part of application of D.psi */
  /* tmp1 <- H * psi */
  /*      <- gamma_5 * M * psi */
  (*M)( ltmp, psi_proj, PLUS);
  tmp1 = Gamma(G5)*ltmp;
  Double c = norm2(tmp1);

  /* If exactly 0 norm, then solution must be 0 (for pos. def. operator) */
  if ( toBool(c == Real(0)) ) { 
    chi = zero;
    return;
  }

  
  /*  chi  +=  constP*(gamma_5 + ichiral) * H * Psi  */
  switch (ichiral) {
  case CH_PLUS:
    ltmp += tmp1;
    break;

  case CH_MINUS:
    ltmp -= tmp1;
    break;

  default:
    QDP_error_exit("unsupported chirality", ichiral);
  }


  chi += constP * ltmp;

 // Effectively the multi mass solve is on tmp1 
  // *******************************************************************
  // Solve  (MdagM + rootQ_n) chi_n = H * tmp1

  Double rsd_sq = c * RsdCG*RsdCG; // || tmp 1 ||^2 * epsilon^2
  Double cp;
  Double d; // InnerProduct

  /* If exactly 0 norm, then solution must be 0 (for pos. def. operator) */
  if (toBool(c == 0))
  {
    chi = zero;
    return;
  }


  LatticeFermion Ap;
  LatticeFermion r;
  LatticeFermion p;

  Double a[MaxCG+1];              // Alpha for unshifted (isz) system
  Double b[MaxCG+1];              // Beta for unshifted system

  bool convP;
  int  iters_taken;

  // By default, rootQ(isz) is considered the smallest shift 
  int isz = numroot-1;             // isz identifies system with smalles shift
  

  // r[0] := p[0] := tmp1 
  r = tmp1;
  p = tmp1;


  cp = 0;
  convP = false;

  
  Double alpha_minus_one=1;
  b[0] = 0;         // b[0] -- do we ever use ?
  


  int k;
  Double c_iter[MaxCG+1];

  // Do the iterations proper 
  // first pass
  for(k = 0; k < MaxCG && ! convP ; ++k) {

    // Keep hold of the residuum
    c_iter[k] = c;

    (*MdagM)(Ap, p, PLUS);
    Ap += p * rootQ[isz];

    // Project out eigenvectors
    if (k % ReorthFreq == 0) {
      GramSchm(Ap, EigVec, NEig);
    }

    //  d =  < p, A.p >
    d = innerProductReal(p, Ap);                       // 2 Nc Ns  flops 
    

    // k+1 because a[0] corresponds to a_{-1}
    a[k] = c/d;
    r -= Real(a[k])*Ap;

    // Project out eigenvectors 
    if (k % ReorthFreq == 0) {
      GramSchm (r, EigVec, NEig);
    }
    
    cp = c;
    c = norm2(r);

    b[k+1] = c/cp;

    p = r + Real(b[k+1])*p;
 
    convP = toBool( c < rsd_sq );
  }


  int niters = k;

  // OK First pass done. I now need to compute gamma_j and c_j
  multi2d<Double> gamma(niters+1,numroot);
  for(int s = 0; s < numroot; s++) { 
    gamma[0][s] = 1;                     // Really gamma_0
  }
  
  for(int j=0; j < niters; j++) {
    for(int s =0; s < numroot; s++) {

      if( s != isz ) { 
	Double a_minus;

	if( j == 0 ) { 
	  // alpha_minus_one -- no need to store. Only need alpha_0,...
	  // in second pass
	  a_minus = alpha_minus_one;
	}
	else {
	  a_minus = a[j-1];
	}

	Double ga = gamma[j][s];
	Double ga_minus;
	if( j == 0 ) { 
	  // gamma[-1][s] -- no need to store. Only need gamma[0][s]...
	  // in second pass
	  ga_minus =1;
	}
	else{ 
	  ga_minus = gamma[j-1][s];
	}

	Double tmp_num = ga*ga_minus*a_minus;
	Double tmp_den = ga_minus*a_minus*(Double(1) + a[j]*(rootQ[s] - rootQ[isz]));
        tmp_den += a[j]*b[j]*(ga_minus - ga);

	gamma[j+1][s] = tmp_num/tmp_den;
 
      }
    }

  } 

  multi1d<Double> sumC(niters+1);

  for(int j=0; j<=niters; j++) { 

    sumC[j] = 0;

    for(int m=0; m <= niters-j; m++) { 
      
      Double qsum=resP[isz];

      for(int s =0; s < numroot; s++) { 

	if( s != isz ) { 
	  
	  qsum += resP[s]*gamma[m+j+1][s]*gamma[m+j][s]/gamma[j][s];
	  
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

  for(k=0; k < niters;k++) { 
    (*MdagM)(Ap, p, PLUS);
    Ap += p * rootQ[isz];

    // Project out eigenvectors
    if (k % ReorthFreq == 0) {
      GramSchm(Ap, EigVec, NEig);
    }

    // Update chi
    chi += Real(sumC[k])*r;

    // Update r
    r -= Real(a[k])*Ap;

    // Project out eigenvectors 
    if (k % ReorthFreq == 0) {
      GramSchm (r, EigVec, NEig);
    }

    p = r + Real(b[k+1])*p;
  }




  /* Rescale to the correct normalization */
  /*  chi <-  (1/4)*(1-m_q^2) * chi  */
  ftmp = Real(0.25) * (Real(1) - m_q*m_q);
  chi *= ftmp;	        /* 2 Nc Ns  flops */
  //  QDPIO::cout << "Overlap Inner Solve (lovddag_double_pass(" << ichiral << ")) = " << n_count << " iterations" << endl;
}

