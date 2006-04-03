// $Id: ritz_array.h,v 3.0 2006-04-03 04:58:57 edwards Exp $
#ifndef __ritz_array_h__
#define __ritz_array_h__

#include "chromabase.h"
#include "linearop.h"

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

 *  psi			New eigenvector
 *  p			Direction vector
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
void Ritz_t(const LinearOperatorArray<T>& A, // Herm Pos Def
	    Real& lambda,               // Current E-value
	    multi2d<T>& psi_all,        // E-vector array
	    int N_eig,                  // Current evec index
	    const Real& Rsd_r,          // Target relative residue
	    const Real& Rsd_a,          // Target absolute residue
	    const Real& zero_cutoff,    // If ev slips below this we consider
	                                // if to be zero
	    int n_renorm,               // Renormalise frequency
	    int n_min, int n_max,       // minimum / maximum no of iters to do
	    int MaxCG,                  // Max iters after which we bomb
	    bool ProjApsiP,             // Project A (?) -- user option
	    int& n_count,               // No of iters actually taken
	    Real& final_grad,           // Final gradient
	    bool Kalk_Sim,              // Are we in Kalk Simma mode?
	    const Real& delta_cycle,    // Initial error estimate (KS mode)
	    const Real& gamma_factor) ;  // Convergence factor Gamma


void Ritz(const LinearOperatorArray<LatticeFermion>& A,   // Herm Pos Def
	  Real& lambda,                            // Current E-value
	  multi2d<LatticeFermion>& psi_all,        // E-vector array
	  int N_eig,                  // Current evec index
	  const Real& Rsd_r,          // Target relative residue
	  const Real& Rsd_a,          // Target absolute residue
	  const Real& zero_cutoff,    // if ev slips below this we consider 
	                              // if to be zero.
	  int n_renorm,               // Renormalise frequency
	  int n_min, int n_max,       // minimum / maximum no of iters to do
	  int MaxCG,                  // Maximum no if iters after which we bomb
	  bool ProjApsiP,             // Project A (?) -- user option
	  int& n_count,               // No of iters actually taken
	  Real& final_grad,           // Final Gradient
	  bool Kalk_Sim,              // Are we in Kalk Simma mode?
	  const Real& delta_cycle,    // Initial error estimate (KS mode)
	  const Real& gamma_factor) ;  // Convergence factor Gamma

}  // end namespace Chroma

#endif
