// $Id: lovlapms_w.cc,v 3.0 2006-04-03 04:58:50 edwards Exp $
/*! \file
 *  \brief Overlap-pole operator
 */
#include <math.h>
#include "chromabase.h"
#include "actions/ferm/linop/lovlapms_w.h"
#include "meas/eig/gramschm.h"


#undef LOVLAPMS_RSD_CHK

namespace Chroma 
{ 
void lovlapms::operator() (LatticeFermion& chi, const LatticeFermion& psi, 
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
void lovlapms::operator() (LatticeFermion& chi, const LatticeFermion& psi, 
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
  //
  //  Also we must bear in mind that if we want the dagger of the 
  // operator we must use gamma_5 psi instead of psi
  //
  // at this stage tmp1 holds either psi or gamma_5 psi as required
  // so we must project from tmp1

  if (NEig > 0)
  {
    Complex cconsts;
    LatticeFermion tmp3 = tmp1;

    for(int i = 0; i < NEig; ++i)
    {

      // BUG Should this not be innerProduct(EigVec[i], psi) ???
      //                     or innerProduct(EigVec[i], g5 psi) ????

      cconsts = innerProduct(EigVec[i], tmp3);
      tmp1 -= EigVec[i] * cconsts;

      cconsts *= EigValFunc[i];
      chi += EigVec[i] * cconsts;
    }

    // tmp3 should go out of scope here and be freed.
  }

  // tmp1 <- H * Projected tmp_1, where tmp1 = psi or gamma_5 psi as needed
  //      <- gamma_5 * M * tmp1
  (*M)(tmp2, tmp1, PLUS);
  tmp1 = Gamma(G5) * tmp2;
  
 

  Double c = norm2(tmp1);
  
  /* If exactly 0 norm, then solution must be 0 (for pos. def. operator) */
  if (toBool(c == 0))
  {
    chi = zero;
    END_CODE();
    return;
  }


  
  // *******************************************************************
  // Solve  (MdagM + rootQ_n) chi_n = tmp1 where
  //
  // tmp1 = H psi or H_gamma_5 psi
  //
#ifdef LOVLAPMS_RSD_CHK
  LatticeFermion b_vec = tmp1;
  LatticeFermion x;
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

  // We are solving with 2/(1-mu) D(mu) here so 
  // I should readjust the residuum by (1-mu)/2
  Real epsilon_normalise = epsilon*(Real(1)-m_q)/Real(2);

  // Real target for sign function -- from Wuppertal paper
  Real epsilon_target = epsilon_normalise/(Real(2) + epsilon_normalise);

  // Square it up
  Real rsdcg_sq = epsilon_target*epsilon_target;   // Target residuum squared

  // Get relative target
  Real rsd_sq = norm2(psi)*rsdcg_sq;      // Used for relative residue comparisons
                                          // r_t^2 * || r ||^2


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

  // Initialise search vectors for shifted systems
  for(s = 0; s < numroot; ++s){
    p[s] = tmp1;
  }

  // Set convergence masks to false
  convsP = false;
  convP = false;

  
  iz = 1;  // z_index  z[ iz , s ] holds zeta(s)
           //          z[ 1-iz, s ] holds zeta_minus(s)

  z = 1;   // This fills both zeta and zeta_minus

  a = 0;   // Alpha for unshifted
  b = 1;   // beta for unshifted 

  int k;
  // Do the iterations proper 
  for(k = 0; k <= MaxCG && ! convP ; ++k) {

    //  Unshifted beta value: b[k] -- k is iteration index
    //  b[k] := | r[k] |**2 / < p[k], Ap[k] > ; 
    //  First compute  d  =  < p, A.p >  
    //  Ap = A . p 

    // This bit computes 
    // Ap = [  M^dag M + rootQ(isz)  ] p_isz
    

    (*MdagM)(Ap, p[isz], PLUS);
    Ap += p[isz] * rootQ[isz];

    // Project out eigenvectors
    if (k % ReorthFreq == 0) {
      GramSchm(Ap, EigVec, NEig, all);
    }

    //  d =  < p, A.p >
    d = innerProductReal(p[isz], Ap);                       // 2 Nc Ns  flops 
    
    bp = b;                        // Store previous unshifted beta
    b = -Real(c/d);                // New unshifted beta

    // Compute the shifted bs and z 
    bs[isz] = b;


    // iz now points to previous z's
    iz = 1 - iz;
    
    // Compute new shifted beta and zeta values
    // as per Beat Jegerlehner's paper: hep-lat/9612014
    // eqns 2.42, 2.44 on page 7
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


    // New residual of system with smallest shift
    // r[k+1] += b[k] A . p[k] ; 
    r += b * Ap;	        // 2 Nc Ns  flops 

#ifdef LOVLAPMS_RSD_CHK
    x -= b * p[isz];
#endif

    // Project out eigenvectors 
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
    tmp1 = p[isz] * rtmp;	// 2 Nc Ns  flops 

    // Now the other shifts. 
    // Converged systems haven' changed, so we only add results
    // from the unconverged systems
    for(s = 0; s < numroot; ++s) {
      if(s != isz  &&  !convsP[s]) {

	rtmp = bs[s] * resP[s];
	tmp1 += p[s] * rtmp;	// 2 Nc Ns  flops
      }
    }

    // Now update the sgn(H) with the above accumulated linear sum
    chi -= tmp1;                   // 2 Nc Ns  flops

    // Store in cp the previous value of c
    // cp  =  | r[k] |**2 
    cp = c;

    // And compute the current norm of r into the (now saved) c
    //  c  =  | r[k] |**2 
    c = norm2(r);	                // 2 Nc Ns  flops

    // Work out the alpha factor for the system with smallest shift.
    //  a[k+1] := |r[k]|**2 / |r[k-1]|**2 ; 
    a = Real(c/cp);

    // Now update search vectors p_{shift}[k]
    // the system with smallest shift gets updated as
    //
    //  p[k+1] := r[k+1] + a[k+1] p[k];

    // The updated systems get updated as
    //  
    //  ps[k+1] := zs[k+1] r[k+1] + as[k+1] ps[k]; 
    //
    // where the as[]-s are the shifted versions of alpha. 
    // we must first computed these as per Beat's paper hep-lat/9612014
    // eq 2.43 on page 7.
    // 
    // As usual we only update the unconverged systems
    for(s = 0; s < numroot; ++s)
    {

      if (s == isz) {
	// Smallest shift 	  
	// p[k+1] = r[k+1] + a[k+1] p[k]
	// 
	// k is iteration index
	p[s] *= a;	        // Nc Ns  flops 
	p[s] += r;		// Nc Ns  flops 
	
      }
      else {
	if (! convsP[s]) {
	  // Unshifted systems
	  // First compute shifted alpha
	  as = a * z[iz][s]*bs[s] / (z[1-iz][s]*b);
	  
	  // Then update
	  //    ps[k+1] := zs[k+1] r[k+1] + as[k+1] ps[k]; 
	  //
	  // k is iteration index
	  p[s] *= as;	        // Nc Ns  flops 
	  p[s] += r * z[iz][s];	// Nc Ns  flops 
	  
	}
      }
    }

#if 0
    // Project out eigenvectors 
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


#ifdef LOVLAPMS_RSD_CHK
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
      chi_sq_new = rsd_sq * norm2(chi);

      // Get || Delta sgn()
      chi_sq_diff = norm2(tmp1);      // the diff of old and new soln

#if 0
      QDPIO::cout << "Iter " << k << " || delta Sgn() || " << sqrt(chi_sq_diff) << endl;
#endif 
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

  QDPIO::cout << "Overlap Inner Solve (lovlapms): " << k << " iterations " << endl;
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

