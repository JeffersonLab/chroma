// $Id: lovddag_w.cc,v 3.0 2006-04-03 04:58:50 edwards Exp $
/*! \file
 *  \brief Overlap-pole operator
 */
#include <math.h>
#include "chromabase.h"
#include "actions/ferm/linop/lovddag_w.h"
#include "meas/eig/gramschm.h"


#undef LOVDDAG_RSD_CHK


namespace Chroma 
{ 
  //! Apply the GW operator onto a source vector
  /*! \ingroup linop
   *
   */
  void lovddag::operator() (LatticeFermion& chi, const LatticeFermion& psi, 
			    enum PlusMinus isign) const
  {
    operator()(chi, psi, isign, RsdCG );
  }

  //! Apply the GW operator onto a source vector
  //  epsilon specifies desired precision (absolute)
  /*! \ingroup linop
   *
   */
  void lovddag::operator() (LatticeFermion& chi, const LatticeFermion& psi, 
			    enum PlusMinus isign, Real epsilon) const
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
    chi = init_fac * psi;


    // Project out eigenvectors of source if desired 
    // chi  +=  eps(lambda) * (gamma_5 + ichiral) * EigVec  
    // Usually "func(.)" is eps(.); it is precomputed in EigValFunc.
    if (NEig > 0)
      {
	Complex cconsts;

	for(int i = 0; i < NEig; ++i)
	  {
      
	    cconsts = innerProduct(EigVec[i], psi);
      
	    psi_proj -= cconsts * EigVec[i];

   
	    // Construct  tmp1 = (gamma_5 +/- 1) * EigVec 
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
		error_message << "Should not use lovddag unless the chirality is either CH_PLUS or CH_MINUS" << endl;
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
      END_CODE();
      return;
    }

  
    /*  chi  +=  constP*(gamma_5 +/-1) * H * Psi  */
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


    /********************************************************************/
    /* Solve  (MdagM + rootQ_n) chi_n = H * tmp1 */
#ifdef LOVDDAG_RSD_CHK
    // DEBUG 
    LatticeFermion b_vec=tmp1;
    LatticeFermion x = zero;
#endif
    LatticeFermion Ap;
    LatticeFermion r; 
    multi1d<LatticeFermion> p(numroot);

    Real a;              // Alpha for unshifted (isz) system
    Real as;             // alpha for current shifted system
    Real b;              // Beta for unshifted system
    Real bp;             // Beta previous for unshifted system

    Double ztmp;           // Assorted reals (shifted residues)
    Double cp;
    Double d;
    Real z0;                    // temporary value for zeta previous
    Real z1;                    // temporary value for zeta current


    multi1d<Real> bs(numroot);  // beta for shifted system
    multi2d<Real> z(2,numroot); // zeta for shifted system


    Double chi_sq_new;              // sgn() convergence criterion
    Double chi_sq_diff;
    multi1d<bool> convsP(numroot);  // convergence mask for shifted system
    bool convP;                     // overall convergence mask
  
    int iz;                     //  Temporary index for getting at 
    //  zeta values in z array 

    int s;                      // Counter for loops over shifts

 
    // We are dealing with 4/(1-m_q)^2
    // so I should readjust the residua by that squared
    Real epsilon_normalise = epsilon*(Real(1)-m_q*m_q)/Real(4);

    // Now get the desired residuum to terminated the CG
    //
    // Which is epsilon/(2+epsilon) according to the wuppertal paper
    Real epsilon_target = epsilon_normalise/(Real(2)+epsilon_normalise);

    // Now square that up
    Real rsdcg_sq = epsilon_target*epsilon_target;

    Real rsd_sq = norm2(psi)*rsdcg_sq;      // Used for relative residue comparisons
    // r_t^2 * || r ||^2

    /* By default (could change), rootQ(isz) is considered the smallest shift */
    int isz = numroot-1;

          
    /* r[0] := p[0] := tmp1 */
    r = tmp1;
  
    for(s = 0; s < numroot; ++s) {
      p[s] = tmp1;
    }

    // Set convergence masks to fals
    convsP = false;
    convP = false;

    iz = 1;  // z_index  z[ iz , s ] holds zeta(s)
    //          z[ 1-iz, s ] holds zeta_minus(s)

    z = 1;   // This fills both zeta and zeta_minus

    a = 0;   // Alpha for unshifted
    b = 1;   // beta for unshifted 

    int k;
    for(k = 0; k <= MaxCG && ! convP ; ++k)
      {

	//  Unshifted beta value: b[k] -- k is iteration index
	//  b[k] := | r[k] |**2 / < p[k], Ap[k] > ; 
	//  First compute  d  =  < p, A.p >  
	//  Ap = A . p 

	// This bit computes 
	// Ap = [  M^dag M + rootQ(isz)  ] p_isz
	(*MdagM)(Ap, p[isz], PLUS);
	Ap += p[isz] * rootQ[isz];


	/* Project out eigenvectors */
	if (k % ReorthFreq  == 0){
	  GramSchm(Ap, EigVec, NEig, all);
	}

	/*  d =  < p, A.p >  */
	d = innerProductReal(p[isz], Ap);            /* 2 Nc Ns  flops */
    
	bp = b;                            // Store previous unshifted beta
	b = -Real(c/d);                    // New unshifted beta

	/* Compute the shifted bs and z */
	bs[isz] = b;

	// iz now points to previous z's
	iz = 1 - iz;
    
	// Compute new shifted beta and zeta values
	// as per Beat Jegerlehner's paper: hep-lat/9612014
	// eqns 2.42, 2.44 on page 7
	for(s = 0; s < numroot; s++) {

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
    
	// New residual of system with smallest shift
	// r[k+1] += b[k] A . p[k] ; 
	r += b * Ap;	        // 2 Nc Ns  flops 


	/* Project out eigenvectors */
	if (k % ReorthFreq == 0) {
	  GramSchm (r, EigVec, NEig, all);
	}
    
	// Work out new iterate for sgn(H).
	//
	// This in effect updates all x vectors and performs
	// immediately the linear sum.
	// However, the x's are never explicitly computed. Rather
	// the changes they would get get rolled onto the chi immediately
	//
	//  chi[k+1] -= sum_{shifts} resP_{shift} b_{shift}[k] p_{shift}[k] ; 
	//
	// since we are doing the linear combinations we multiply in the 
	// constant in the numerator too.

	// smallest shift first
	Real(rtmp);
	rtmp = resP[isz] * b;      
	ltmp = p[isz] * rtmp;	// 2 Nc Ns  flops 

	// Now the other shifts. 
	// Converged systems haven' changed, so we only add results
	// from the unconverged systems
	for(s = 0; s < numroot; ++s) {
	  if(s != isz  &&  !convsP[s]) {

	    rtmp = bs[s] * resP[s];
	    ltmp += p[s] * rtmp;	// 2 Nc Ns  flops
	  }
	}


    
	/*  diff :=  (gamma_5 + ichiral) * tmp1  */
	tmp1 = Gamma(G5)*ltmp;

	switch( ichiral ) { 
	case CH_PLUS: 
	  tmp1 += ltmp;
	  break;
	case CH_MINUS:
	  tmp1 -= ltmp;
	  break;
	default :
	  QDP_error_exit("lovddag: isign must be CH_PLUS or CH_MINUS\n");
	  break;
	}


	/* chi += diff . The minus comes from the CG method. */
	chi -= tmp1;                   /* 2 Nc Ns  flops */

#ifdef LOVDDAG_RSD_CHK
	/*  diff :=  (gamma_5 + ichiral) * tmp1  */
	x -= b*p[isz];
#endif

	// Store in cp the previous value of c
	// cp  =  | r[k] |**2 
	cp = c;

	// And compute the current norm of r into the (now saved) c
	//  c  =  | r[k] |**2 
	c = norm2(r);	                // 2 Nc Ns  flops

	// Work out the alpha factor for the system with smallest shift.
	//  a[k+1] := |r[k]|**2 / |r[k-1]|**2 ; 
	a = Real(c/cp);

	/*  p[k+1] := r[k+1] + a[k+1] p[k]; */
	/*  Compute the shifted as */
	/*  ps[k+1] := zs[k+1] r[k+1] + a[k+1] ps[k]; */
	for(s = 0; s < numroot; ++s) {
	  if (s == isz) {
	
	    p[s] *= a;	         // Nc Ns  flops 
	    p[s] += r;		 // Nc Ns  flops 

	
	  }
	  else {
	    if (! convsP[s]) {	  
	      as = a * z[iz][s]*bs[s] / (z[1-iz][s]*b);
	  
	      p[s] *= as;	/* Nc Ns  flops */
	      p[s] += r * z[iz][s];	/* Nc Ns  flops */
	  
	    }
	  }
	}

#if 0
	/* Project out eigenvectors */
	if (k % ReorthFreq == 0)
	  GramSchm (p, numroot, EigVec, NEig, all);
#endif

	// Convergence tests start here.
	//
	// These are two steps:
	//
	// i) Check that the shifted vectors have converged 
	//    by checking their accumulated residues
	//
	// ii) Check that the sign function itself has converged
	//
	// 
	// We set the global convergence flag to true, and then if 
	// any vectors are unconverged this will flip the global flag back
	// to false
	convP = true;                          // Assume convergence and prove
	// otherwise

	for(s = 0; s < numroot; ++s) {
      
	  // Only deal with unconverged systems
	  if (! convsP[s]) {
	
	    // Compute the shifted residuum squared
	    //
	    //  r_{shift} = || r || * zeta(shift)
	    //
	    //  hence || r_shift ||^2 = || r ||^2 * zeta_shift^2
	    //
	    // || r^2 || is already computed in c
	    //
	    // Store  || r_shift ||^2 in ztmp
	    ztmp = Real(c) * z[iz][s]*z[iz][s];	

	    // Check ztmp is smaller than the target residuum
	    bool btmp = toBool(ztmp < rsd_sq);
	    convP = convP & btmp;
	    convsP[s] = btmp;
	  }
	}
    
#ifdef LOVDDAG_RSD_CHK
	if(convP) {
	  LatticeFermion tmp_normcheck;
	  (*MdagM)(tmp_normcheck, x, PLUS);
	  tmp_normcheck += rootQ[isz]*x;
	  tmp_normcheck -= b_vec;
	  Double norm2check = norm2(tmp_normcheck);
	  Double check_ztmp = Real(c) * z[iz][isz]*z[iz][isz];
	  QDPIO::cout << "|| b - (Q_isz + MM)x || = " << norm2check << " accum = " << check_ztmp << endl;
	}
#endif

	// Now check convergence of the sgn() itself.
	// It was updated with 
	//               sum resP(shift) beta(shift) p_shift
	// and this quantity is still in tmp1.
	//
	// So norm tmp1 is like an abolute error in sgn() aka: Delta sgn()
	//
	// Here we ensure Delta sgn() < target r^2 || sgn H[k-1] ||
	//
	// ie that  the relative error in sgn() is smaller than the target.
	// sgn H is kept in chi
	//
	// Only converge if chi is converged. If vectors converge first, then error
    
#if 0  
	if (k > 0 &&  !convP) {

	  // Get target r^2 * || sgn (H) ||^2
	  chi_sq_new = rsd_sq;

	  // Get || Delta sgn() ||^2
	  chi_sq_diff = norm2(tmp1);      // the diff of old and new soln

	  // Check convergence
	  bool btmp = toBool(chi_sq_diff < chi_sq_new);

	  // If we havent converged globally but the vectors have then error
	  if (! btmp && convP) {
	    QDP_error_exit("vectors converged but not final chi");
	  }
      
	  // cnvP = convP & btmp; 
	  convP = btmp;
	}
#endif
      }

    int n_count = k;

    /* Rescale to the correct normalization */
    /*  chi <-  (1/4)*(1-m_q^2) * chi  */
    ftmp = Real(0.25) * (Real(1) - m_q*m_q);
    chi *= ftmp;	        /* 2 Nc Ns  flops */
    QDPIO::cout << "Overlap Inner Solve (lovddag(" << ichiral << ")) = " << n_count << " iterations" << endl;

    END_CODE();
  }

}; // End Namespace Chroma

