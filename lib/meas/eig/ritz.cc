/*! \file
 *  \brief Ritz code for eigenvalues
 */

#include "chromabase.h"
#include "meas/eig/ritz.h"
#include "meas/eig/gramschm.h"

namespace Chroma {

//! Minimizes the Ritz functional with a CG based algorithm
/*!
 * \ingroup eig
 *
 * This subroutine minimizes the Ritz functional with a CG based
 * algorithm to find the n-th lowest eigenvalue of a hermitian A

 * lambda_n = min_z <z|Az>/<z|z>

 * In the subspace orthogonal to the (n-1) lower eigenvectors of A

 * This routine can handle the original version or be a part of a 
 * Kalkreuter Simma algorithm.
 * 
 * The algorithm has been modified that there is now only one target
 * residue -- the relative one.
 *
 *  so the original stopping criterion is now:
 * 
 *         || g || < Rsd_r * lambda
 *
 * The two key features for Kalkreuter - Simma over the normal procedure
 * are: 
 *   i) Adaptive termination of the procedure (controlled from caller)
 * 
 *   The iteration will stop if the gradient has decreased by at least
 * a given factor of gamma^{-1} which according to the paper is O(10)
 *
 *  ii) The initial error estimate is now the "delta_cycle" criterion.
 *  basically this scales the norm of || g^2 ||. So the iteration will
 *  also stop if   delta_cycle || gamma^2 || < Rsd_r mu.
 *
 * To run in non K-S mode, set gamma = 1, delta_cycle = 1.
 *
 * For details of the Kalkreuter Simma algorithm see 
 * hep-lat/9507023
 *
 * Algorithm:

 *   Apsi[0] :=  A . Psi[0] ;
 *   mu[0]    :=  < Psi[0] | A . Psi[0] > ; 	Initial Ritz value
 *             Note: we assume  < Psi[0] |Psi[0] > = 1
 *   p[1] = g[0]   :=  ( A - mu[0] ) Psi[0] ;	Initial direction
 *   if converged then return
 *   FOR k FROM 1 TO MaxCG DO			CG iterations
 *       s1 = (mu[k-1] + < p[k] | A p[k] >/p[k]|^2) / 2;
 *       s2 = (mu[k-1] - < p[k] | A p[k] >/p[k]|^2) / 2;
 *       s3 = |g[k-1]|^2 / |p[k]|;
 *       Compute a and cos(theta), sin(theta)
 *          (see DESY internal report, September 1994, by Bunk, Jansen,
 *           Luescher and Simma)
 *       lambda = mu[k] = mu[k-1] - 2 a sin^2(theta);
 *       Psi[k] = cos(theta) Psi[k-1] + sin(theta) p[k]/|p[k]|;
 *       Apsi[k] = cos(theta) Apsi[k-1] + sin(theta) A p[k]/|p[k]|;
 *       g[k] = Apsi[k] - mu[k] Psi[k]
 *      
 *       if now converged, then exit
 * 
 *       b[k+1] := cos(theta) |g[k]|^2 / |g[k-1]|^2;
 *       p[k+1] := g[k] + b[k+1] (p[k] - Psi[k] < Psi[k] | p[k] >);	New direction

 * Arguments:

 *  A		The hermitian Matrix as a lin op	(Read)
 *  lambda      Eigenvalue to find              (Write)
 *  Psi_all	Eigenvectors			(Modify)
 *  N_eig	Current Eigenvalue number 	(Read)
 *  RsdR_r	(relative) residue		(Read)
 *  n_renorm	Renormalize every n_renorm iter.(Read)
 *  n_min       Minimum no of CG iters to do    (Read)
 *  n_max       Maximum no of CG iters to do    (Read)
 *  n_count	Number of CG iteration	(done)	(Write)
 *  Kalk_Sim	Use Kalkreuter-Simma criterion  (Read)
 *  delta_cycle The initial error estimate      (Read)
 *  gamma	"Convergence factor" gamma required	(Read)
 *  ProjApsiP	flag for projecting A.psi	(Read)

  * Local Variables:

 *  psi			New eigenstd::vector
 *  p			Direction std::vector
 *  Apsi		Temporary for  A.psi
 *  Ap			Temporary for  A.p, and other
 *  mu			Ritz functional value
 *  g2			| g[k] |**2
 *  p2			| p[k] |**2
 *  k			CG iteration counter
 *  b                   beta[k+1]
 *  and some others...
 */

template < typename T >
void Ritz_t(const LinearOperator<T>& A, // Herm Pos Def
	    Real& lambda,               // Current E-value
	    multi1d<T>& psi_all,        // E-std::vector array
	    int N_eig,                  // Current evec index
	    const Real& Rsd_r,          // Target relative residue
	    const Real& Rsd_a,          // Target absolute residue
	    const Real& zero_cutoff,    // if ev-slips below this 
	    // we consider it to be a zero
	    int n_renorm,               // Renormalise frequency
	    int n_min, int n_max,       // minimum / maximum no of iters to do
		 int MaxCG,                  // Max iters after which we bomb
	    bool ProjApsiP,             // Project A (?) -- user option
	    int& n_count,               // No of iters actually taken
	    Real& final_grad,           // || g || at the end
	    bool Kalk_Sim,              // Are we in Kalk Simma mode?
	    const Real& delta_cycle,    // Initial error estimate (KS mode)
	    const Real& gamma_factor)   // Convergence factor Gamma
{
  START_CODE();
  
  T psi;
  T p;
  T Apsi;
  T Ap;

  const Subset& s = A.subset(); // Subset over which A acts

  //  Make Psi_all(N_eig-1) orthogonal to the previous eigenvectors  */
  int N_eig_index = N_eig - 1;
  int N_eig_minus_one = N_eig_index;

  // For now equate worry about concrete subsetting later
  psi[s] = psi_all[N_eig_index];

  // Project out subspace of previous 
  if( N_eig_minus_one > 0 )
    GramSchm(psi, psi_all, N_eig_minus_one, s);

  // Normalize
  Double d = sqrt(norm2(psi,s));
  psi[s] /= d;

  // Now we can start
  //  Apsi[0]   :=  A . Psi[0]
  A(Apsi, psi, PLUS);

  // Project to orthogonal subspace, if wanted 
  // Should not be necessary, following Kalkreuter-Simma **/
  if( N_eig_minus_one > 0 ) {
    if (ProjApsiP) {
      GramSchm(Apsi, psi_all, N_eig_index, s);
    }
  }
  
  //  mu  := < Psi[0] | A Psi[0] > 
  Double mu = innerProductReal(psi, Apsi, s);

  //  p[0] = g[0]   :=  ( A - mu[0] ) Psi[0]  =  Apsi - mu psi */
  lambda = Real(mu);


  // Psi and Apsi are both projected so p should also be (?)
  p[s] = Apsi - lambda * psi;
    
  //  g2 = p2 = |g[0]|^2 = |p[0]|^2
  Double g2;
  Double p2;
  Double g2_0;
  Double g2_p;

  g2 = norm2(p,s);
  
  // Keep hold of initial g2
  g2_0 = g2;
  p2 = g2;

#if 0
  // Debugging
  QDPIO::cout << "Starting Ritz: N_eig=" << N_eig << ", mu = " << mu
	      << ", g2_0 = " << g2_0 << std::endl;
#endif

  // Check whether we have converged
  bool convP;
  bool minItersDoneP = n_min <= 0;
  bool maxItersDoneP = n_max <= 0;
  bool zeroFoundP;
  bool CGConvP;
  bool KSConvP;
  bool deltaCycleConvP;

  Double rsd_r_sq;
  Double rsd_a_sq;
 
  // Make up the stopping criteria
  rsd_a_sq = Rsd_a * Rsd_a;                // Absolute
  rsd_r_sq = Rsd_r * Rsd_r * mu * mu;      // Relative
  
  // Converged if g2 < Max( absolute_target, relative-target) 
  // This is necessary for low lying non-zero modes, where 
  // the relative error can reach below machine precision
  // (ie g2 never gets smaller than it )
  CGConvP = toBool( g2 < rsd_r_sq  ) || toBool(g2 < rsd_a_sq );

  // Zero Ev-s can reach the cutoff value quickly, also 
  // the absolute error doesn't apply to them.
  // I could try and use the absolute error as a zero cutoff...
  // But lets allow this flexibility 
  zeroFoundP = toBool( fabs(mu) < zero_cutoff );

  // error based on delta_cycle
  // NOT USED AT THIS TIME BUT POSSIBLY IN THE FUTURE
  Double delta_cycle_err = delta_cycle;

  if( Kalk_Sim == false ) { 

    // Normal CG stopping criterion
    // The stopping criterion is satisfied AND 
    // we have performed the minimum number of iterations
    convP = minItersDoneP && ( CGConvP || zeroFoundP );
  }
  else { 
    // Kalk-Sim stopping criterion
    // The iteration terminates if
    //    
    //    i) The minimum no of iters is done and max iters is also done
    //    (this is iteration 0 so I ignore max iters at this pint
    // 
    // AND either the following is true
    //
    //      i) The delta cycle bound is reached  (*)
    //     ii) If the normal CG criterion is satisfied
    //    iii) If mu has reached the zero cutoff
    //
    //  (*) KS delta cycle bount is not in use at this time.
    //KSConvP = toBool( delta_cycle_err < rsd );
    KSConvP = false;         // Not in use

    convP = minItersDoneP && ( KSConvP || CGConvP || zeroFoundP );
  }

  if ( convP ) {
    n_count = 0;
    END_CODE();
    return;
  }


  Double s1, s2, s3, a, b ;
  const Double one = Double(1);
  const Double two = Double(2);
  const Double half = Double(0.5);
  Complex xp;

  Double st, ct; // Sin theta, cos theta
  
  // FOR k FROM 1 TO MaxCG do 
  // MaxCG is an absolute max after which we bomb
  for(int k = 1; k <= MaxCG; ++k)
  {
    //  Ap = A * p  */
    A(Ap, p, PLUS);

    //  Project to orthogonal subspace, if wanted 
    if( N_eig_minus_one > 0 ) { 
      if (ProjApsiP) {
	GramSchm(Ap, psi_all, N_eig_minus_one, s);
      }
    }

    //  d = < p | A p > 
    d = innerProductReal(p, Ap, s);

    // Minimise Ritz functional along circle.
    // Work out cos theta and sin theta
    // See internal report cited above for details
    d = d / p2;
    s1 = half * (mu+d);
    s2 = half * (mu-d);
    p2 = sqrt(p2);
    p2 = one / p2;
    s3 = g2 * p2;
    a = fabs(s2);

    if( toBool (a >= s3) )
    {
      d = s3 / s2;
      d = one + d*d;
      d = sqrt(d);
      a *= d;
    }
    else
    {
      d = s2 / s3;
      d = one + d*d;
      d = sqrt(d);
      a = s3 * d;
    }

    s2 /= a;			// Now s2 is cos(delta) 
    s3 /= a;			// Now s3 is sin(delta) 

    if( toBool( s2 > 0)  )
    {
      s2 = half * (one+s2);
      d  = -sqrt(s2);		// Now d is sin(theta) 
      s2 = -half * s3 / d;	// Now s2 is cos(theta)
    }
    else
    {
      s2 = half * (one-s2);
      s2 = sqrt(s2);		// Now s2 is cos(theta)
      d = -half * s3 / s2;	// Now d is sin(theta)
    }

    // mu[k] = mu[k-1] - 2 a d^2
    mu -= (two * a * d * d);
    lambda = Real(mu);

    // del_lamb is the change in lambda 
    // used in an old stopping criterion, which 
    // we may have to revisit
    // del_lam = Real(s1);

    st = d*p2;		       // Now st is sin(theta)/|p|
    ct = s2;	               // Now ct is cos(theta)

    //  Psi[k] = ct Psi[k-1] + st p[k-1] 
    //  Apsi[k] = ct Apsi[k-1] + st Ap
    psi[s]  *= Real(ct);
    psi[s]  += p * Real(st);
    Apsi[s] *= Real(ct);
    Apsi[s] += Ap * Real(st);
    

    //  Ap = g[k] = Apsi[k] - mu[k] Psi[k]
    Ap[s] = Apsi - psi*lambda;
    

    //  g2  =  |g[k]|**2 = |Ap|**2
    s1 = g2;			// Now s1 is |g[k-1]|^2
    g2_p = g2;                  // g2_p is g2 "previous" 
    g2 = norm2(Ap,s);

    // Convergence test
    // Check for min_KS / max_KS iters done
    minItersDoneP = ( k >= n_min );
    maxItersDoneP = ( k >= n_max );
    
    // Conservative CG test
    rsd_a_sq = Rsd_a*Rsd_a;
    rsd_r_sq = Rsd_r*Rsd_r*mu*mu;

    // Converged if g2 < Max( absolute_target, relative-target) 
    // This is necessary for low lying non-zero modes, where 
    // the relative error can reach below machine precision
    // (ie g2 never gets smaller than it )
    CGConvP = toBool( g2 < rsd_a_sq ) || toBool( g2 < rsd_r_sq);

    // Zero Ev-s can reach the cutoff value quickly, also 
    // the absolute error doesn't apply to them.
    // I could try and use the absolute error as a zero cutoff...
    // But lets allow this flexibility 
    zeroFoundP = toBool( fabs(mu) < zero_cutoff );
    
    if( Kalk_Sim == false ) {

      // Non Kalk Sim criteria. We have done the minimum 
      // Number of iterations, and we have either done the maximum
      // number too or we have converged
      convP = minItersDoneP && ( maxItersDoneP || CGConvP || zeroFoundP );
    }
    else { 
      // Kalk-Sim stopping criterion
      // The iteration terminates if
      //    
      //    i) The minimum no of iters is done
      // AND either the following is true
      //      i) KS_max iters is also done
      //     ii) The delta cycle bound is reached  (*)
      //    iii) If the normal CG criterion is satisfied
      //     vi) If mu has reached the zero cutoff
      //
      //  (*) KS delta cycle bount is not in use at this time.
      
      /* The following commented lines were a stab at the delta_cycle bound
	 
      // work out decrease in || g^2 ||.
      Double g_decr_factor = g2 / g2_p;
      
      // "Extrapolate the delta_cycle error with the decrease in g"
      delta_cycle_err *= fabs(g_decr_factor);
      
      deltaCycleConvP = toBool( delta_cycle_err < rsd );
      */

      deltaCycleConvP = false;
      Double g2_g0_ratio = g2 / g2_0;
      KSConvP = deltaCycleConvP || toBool( g2_g0_ratio <= Double(gamma_factor)) ;
      convP = minItersDoneP && ( maxItersDoneP || CGConvP || KSConvP || zeroFoundP );

    }

    
    if( convP ) { 
      n_count = k;

      // Project onto orthog subspace
      if( N_eig_minus_one > 0 ) { 
	GramSchm(psi, psi_all, N_eig_minus_one, s);
      }

      // Renormalise
      d = sqrt(norm2(psi,s));
      psi[s] /= d;
      d -= one;
      

      // Print out info about convergence 
#if 0
      QDPIO::cout << "Converged at iter=" << k << ", lambda = " << lambda
		  << ",  rsd | mu | = " << rsd << ",  || g || = "
		  << sqrt(g2) << " || x || - 1 = " << d << std::endl;
      
      if(Kalk_Sim) 
      { 
	// Extra info for KalkSimma Mode
	QDPIO::cout << "KS: gamma = "<< gamma_factor << ",  || g ||^2/|| g_0 ||^2="
		    << g2/g2_0 
		    << ",  delta_cycle_err=" << delta_cycle_err << std::endl;
	QDPIO::cout << "KS: CGConvP = " << CGConvP << ",  KSConvP = " << KSConvP << " deltaCycleConvP = " << deltaCycleConvP << std::endl;
      }
#endif
      // Recompute lambda
      // Apsi = A psi 
      A(Apsi, psi, PLUS);
      s1 = mu;
      mu = innerProductReal(psi, Apsi, s);
      lambda = Real(mu);
#if 0
      QDPIO::cout << "Mu-s at convergence: old " << s1 << " vs " << mu << std::endl;

#endif
      // Copy std::vector back into psi array.
      psi_all[N_eig_index][s] = psi;

      // Work out final gradient 
      Ap[s] = Apsi - psi*lambda;
      final_grad = sqrt(norm2(Ap,s));
      END_CODE();
      return;
    }
     

    // Work out new beta
    /* b = beta[k] = cos(theta) |g[k]|^2 / |g[k-1]|^2 */
    b = ct * g2/g2_p;

    ct *= 0.05 * b;
    d = sqrt(g2);
    ct /= p2*d;

    // If beta gets very big, 
    if( toBool( ct > Double(1)) ) 
    {
      /* Restart: p[k] = g[k] = Ap */
      QDPIO::cout << "Restart at iter " << k << " since beta = " << b << std::endl;
      p[s] = Ap;
    }
    else    {
      /* xp = < Psi[k] | p[k-1] > */
      xp = innerProduct(psi, p, s);

      /* p[k] = g[k] + b (p[k-1] - xp psi[k]) */
      p[s] -= psi * xp;
      p[s] *= b;
      p[s] += Ap;
    }

    if( k % n_renorm == 0 )
    {
      /* Renormalize, and re-orthogonalize */
      /*  Project to orthogonal subspace  */
      if( N_eig_minus_one > 0 ) {
	GramSchm (psi, psi_all, N_eig_minus_one, s);
      }
      /* Normalize */
      d = sqrt(norm2(psi,s));

      ct = Double(1)/d;
      psi *= ct;
      d -= one;

      //  Apsi  :=  A . Psi 
      A(Apsi, psi, PLUS);

      // Project to orthogonal subspace, if wanted 
      // Should not be necessary, following Kalkreuter-Simma
      if (ProjApsiP) {
	GramSchm (Apsi, psi_all, N_eig_index, s);
      }

      //  mu  := < Psi | A Psi > 
      mu = innerProductReal(psi, Apsi, s);

      /*  g[k] = Ap = ( A - mu ) Psi  */
      lambda = Real(mu);
      Ap[s] = Apsi - psi * lambda;
      
      /*  g2  =  |g[k]|**2 = |Ap|**2 */
      g2 = norm2(Ap,s);

      /*  Project p[k] to orthogonal subspace  */
      if( N_eig_minus_one > 0 ) { 
	GramSchm(p, psi_all, N_eig_minus_one, s);
      }
      
      /*  Make p[k] orthogonal to Psi[k]  */
      GramSchm(p, psi, s);

      /*  Make < g[k] | p[k] > = |g[k]|**2: p[k] = p_old[k] + xp g[k],
	  xp = (g2 - < g | p_old >) / g2; g[k] = Ap */
      p2 = 0;
      ct = Double(1) / g2;
      xp = Real(ct)*(cmplx(Real(g2),Real(0)) - innerProduct(Ap, p, s));

      p[s] += Ap * xp;
    }
    else if (ProjApsiP && N_eig_minus_one > 0)
    {
      /*  Project psi and p to orthogonal subspace  */
      if( N_eig_minus_one > 0 ) { 
	GramSchm (psi,  psi_all, N_eig_minus_one, s);
	GramSchm (p, psi_all, N_eig_minus_one, s);
      }
    }

    /*  p2  =  |p[k]|**2 */
    p2 = norm2(p, s);
  }

  /*  Copy psi back into psi_all(N_eig-1)  */
 
  psi_all[N_eig_index][s] = psi;
  n_count = MaxCG;
  final_grad = sqrt(g2);
  QDPIO::cerr << "too many CG/Ritz iterations: n_count=" << n_count
	      << ", rsd_r =" << sqrt(rsd_r_sq) <<" rsd_a=" << Rsd_a << ", ||g||=" << sqrt(g2) << ", p2=" << p2
	      << ", lambda" << lambda << std::endl;
  QDP_abort(1);
  END_CODE();
}




void Ritz(const LinearOperator<LatticeFermion>& A,   // Herm Pos Def
	  Real& lambda,                            // Current E-value
	  multi1d<LatticeFermion>& psi_all,        // E-std::vector array
	  int N_eig,                  // Current evec index
	  const Real& Rsd_r,          // Target relative residue
	  const Real& Rsd_a,          // Absolute target residue
	  const Real& zero_cutoff,    // If ev-slips below this we consider it
	                              // to be zero
	  int n_renorm,               // Renormalise frequency
	  int n_min, int n_max,       // minimum / maximum no of iters to do
	  int MaxCG,                  // Max no of iters after which we bomb
	  bool ProjApsiP,             // Project A (?) -- user option
	  int& n_count,               // No of iters actually taken
	  Real& final_grad,           // Final gradient at the end.
	  bool Kalk_Sim,              // Are we in Kalk Simma mode?
	  const Real& delta_cycle,    // Initial error estimate (KS mode)
	  const Real& gamma_factor)   // Convergence factor Gamma
{
  Ritz_t(A, lambda, psi_all, N_eig, Rsd_r, Rsd_a, zero_cutoff, 
	 n_renorm, n_min, n_max, MaxCG,
	 ProjApsiP, n_count, final_grad, 
	 Kalk_Sim, delta_cycle, gamma_factor);
}

}  // end namespace Chroma
