// $Id: ritz.cc,v 1.1 2004-01-04 21:56:04 edwards Exp $
/*! \file
 *  \brief Ritz code for eigenvalues
 */

#error "CONVERSION NOT COMPLETE: NEED TO MAKE APPROPRIATE SUBCLASSING"

#include "chromabase.h"
#include "meas/eig/ritz.h"

using namespace QDP;

//! Minimizes the Ritz functional with a CG based algorithm
/*!
 * \ingroup eig
 *
 * This subroutine minimizes the Ritz functional with a CG based
 * algorithm to find the n-th lowest eigenvalue of a hermitian A

 * lambda_n = min_z <z|Az>/<z|z>

 * In the subspace orthogonal to the (n-1) lower eigenvectors of A

 * This version includes the "early" termination criterion of
 * Kalkreuter and Simma.

 * Algorithm:

 *  Apsi[0] :=  A . Psi[0] ;
 *  mu[0]    :=  < Psi[0] | A . Psi[0] > ; 	Initial Ritz value
 * Note: we assume  < Psi[0] |Psi[0] > = 1
 *  p[1] = g[0]   :=  ( A - mu[0] ) Psi[0] ;	Initial direction
 *  IF |g[0]| <= max(RsdR_a, RsdR_r |mu[0]|) THEN RETURN;	Converged?
 *  FOR k FROM 1 TO MaxCG DO			CG iterations
 *      s1 = (mu[k-1] + < p[k] | A p[k] >/p[k]|^2) / 2;
 *      s2 = (mu[k-1] - < p[k] | A p[k] >/p[k]|^2) / 2;
 *      s3 = |g[k-1]|^2 / |p[k]|;
 *      Compute a and cos(theta), sin(theta)
 *      (see DESY internal report, September 1994, by Bunk, Jansen,
 *       Luescher and Simma)
 *      lambda = mu[k] = mu[k-1] - 2 a sin^2(theta);
 *      Psi[k] = cos(theta) Psi[k-1] + sin(theta) p[k]/|p[k]|;
 *      Apsi[k] = cos(theta) Apsi[k-1] + sin(theta) A p[k]/|p[k]|;
 *      g[k] = Apsi[k] - mu[k] Psi[k]
 *      IF |g[k]| <= max(RsdR_a, RsdR_r |mu[k]|)
 *         && |del_lam| < Rsdlam |lambda| THEN RETURN;	Converged?
 *      b[k+1] := cos(theta) |g[k]|^2 / |g[k-1]|^2;
 *      p[k+1] := g[k] + b[k+1] (p[k] - Psi[k] < Psi[k] | p[k] >);	New direction

 * Arguments:

 *  A		The hermitian Matrix as a lin op	(Read)
 *  Psi_all	Eigenvectors			(Modify)
 *  N_eig	Eigenvalue number 		(Read)
 *  lambda	N_eig-th Eigenvalue		(Write)
 *  RsdR_a	(absolute) residue		(Read)
 *  RsdR_r	(relative) residue		(Read)
 *  Rsdlam	relative accuracy of lambda	(Read)
 *  Cutl_zero	Zero cut-off for lambda		(Read)
 *  n_renorm	Renormalize every n_renorm iter.	(Read)
 *  Ncb		Number of sublattices		(Read)
 *  N_Count	Number of CG iteration		(Write)
 *  Kalk_Sim	Use Kalkreuter-Simma criterion  (Read)
 *  N_min	Minimal CG iterations		(Read)
 *  N_max	Maximal CG iterations		(Read)
 *  Cv_fact	"Convergence factor" required	(Read)
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

 * Global Variables:

 *  MaxCG       Maximum number of CG iterations allowed

 * Subroutines:

 *  A		Apply matrix A to vector 
 */


void Ritz(const LinearOperator<LatticeFermion>& A, 
	  multi2d<LatticeFermion>& psi_all(),
	  int N_eig, const Real& lambda, const Real& RsdR_a, const Real& RsdR_r, 
	  const Real& Rsdlam, const Real& Cutl_zero,
	  int n_renorm, int& n_count, bool Kalk_Sim, int N_min, int N_max, 
	  const Real& Cv_fact, bool ProjApsiP)
{
  START_CODE("Ritz");
  
  int NLayers = 1;
  multi1d<LatticeFermion> psi(NLayers);
  multi1d<LatticeFermion> p(NLayers);
  multi1d<LatticeFermion> Apsi(NLayers);
  multi1d<LatticeFermion> Ap(NLayers);
  LatticeComplex lc;
  LatticeComplex lct;
  LatticeReal lr;
  LatticeReal lrt;
  Complex xp;
  Double mu;
  Double p2;
  Double g2;
  Double a;
  Double d;
  Double s1;
  Double s2;
  Double s3;
  Double one;
  Double two;
  Double half;
  Real g2_0;
  Real del_lam;
  Real rsd;
  Real rsda_sq;
  Real rsdr_sq;
  Real ct;
  Real st;
  Real b;
  Real acc;
  int cb;
  int n;
  
  one = 1;
  two = 2;
  half = 0.5;
  acc = pow(fuzz/5.0, 2);
  rsda_sq = RsdR_a * RsdR_a;
  rsdr_sq = RsdR_r * RsdR_r;
  rsda_sq = max(rsda_sq, acc);
                    

  /*  Make Psi_all(N_eig-1) orthogonal to the previous eigenvectors  */
  int NN_eig = N_eig - 1;

  /*  First copy it into a the local psi  */
  for(cb = 0; cb < NLayers; ++cb)
    psi[cb] = psi_all[NN_eig][cb];

  GramSchm (psi, 1, psi_all, NN_eig, NLayers);

  /* Normalize */
  d = sqrt(norm2(psi));

  ct = Real(1) / Real(d);
  for(cb = 0; cb < NLayers; ++cb)
    psi[cb] *= ct;

  /* Now we can start */
  /*  Apsi[0]   :=  A . Psi[0]  */
  A(Apsi, psi, PLUS);

  /*  Project to orthogonal subspace, if wanted  */
  /** Should not be necessary, following Kalkreuter-Simma **/
  if (ProjApsiP)
    GramSchm(Apsi, 1, psi_all, NN_eig, NLayers);

  /*  mu  := < Psi[0] | A Psi[0] >  */
  mu = innerProductReal(psi, Apsi);

  /*  p[0] = g[0]   :=  ( A - mu[0] ) Psi[0]  =  Apsi - mu psi */
  lambda = Real(mu);
  for(cb = 0; cb < NLayers; ++cb)
  {
    p[cb] = Apsi[cb];
    p[cb] -= psi[cb] * lambda;
  }
  
  /*  g2 = p2 = |g[0]|^2 = |p[0]|^2 */
  g2 = norm2(p);
  p2 = g2;
  g2_0 = max(Cv_fact * Real(g2), acc);

#if 1
  if (ProjApsiP)
    QDPIO::cout << "NN_eig= " << NN_eig << "  lambda= " << mu
		<< "  g2=" << g2
		<< "  g2_0=" << g2_0 << endl;
#endif

  /*  IF |g[0]| <= max(RsdR_a, RsdR_r |mu[0]|) THEN RETURN; */
  rsd = rsdr_sq * lambda * lambda;
  rsd = max(rsda_sq, rsd);
  /* if ( ! Kalk_Sim && Real(g2)  <=  rsd ) */
  if ( ! Kalk_Sim && N_min <= 0 && Real(g2)  <=  rsd )
  {
    /*  Copy psi back into psi_all(N_eig-1)  */
    for(cb = 0; cb < NLayers; ++cb)
      psi_all[NN_eig][cb] = psi[cb];

    n_count = 0;
    END_CODE("subroutine");;
    return;
  }

  /*  FOR k FROM 1 TO MaxCG DO */
  for(int k = 1; k <= MaxCG; ++k)
  {
    /*  Ap = A * p  */
    A(Ap, p, PLUS);

    /*  Project to orthogonal subspace, if wanted  */
    /** Should not be necessary, following Kalkreuter-Simma **/
    if (ProjApsiP)
      GramSchm(Ap, 1, psi_all, NN_eig, NLayers);

    /*  d = < p | A p >  */
    d = innerProductReal(p, Ap);

    d = d / p2;
    s1 = half * (mu+d);
    s2 = half * (mu-d);
    p2 = sqrt(p2);
    p2 = one / p2;
    s3 = g2 * p2;
    a = fabs(s2);
    if( a >= s3 )
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

    s2 /= a;			/* Now s2 is cos(delta) */
    s3 /= a;			/* Now s3 is sin(delta) */
    if( s2 > 0 )
    {
      s2 = half * (one+s2);
      d  = -sqrt(s2);		/* Now d is sin(theta) */
      s2 = -half * s3 / d;	/* Now s2 is cos(theta) */
    }
    else
    {
      s2 = half * (one-s2);
      s2 = sqrt(s2);		/* Now s2 is cos(theta) */
      d = -half * s3 / s2;	/* Now d is sin(theta) */
    }

    /* mu[k] = mu[k-1] - 2 a d^2 */
    s1 = two * a * d * d;
    mu -= s1;
    lambda = Real(mu);
    del_lam = Real(s1);

    st = Real(d*p2);		/* Now st is sin(theta)/|p| */
    ct = Real(s2);		/* Now ct is cos(theta) */

    /*  Psi[k] = ct Psi[k-1] + st p[k-1] */
    /*  Apsi[k] = ct Apsi[k-1] + st Ap */
    for(cb = 0; cb < NLayers; ++cb)
    {
      psi[cb] *= ct;
      psi[cb] += p[cb] * st;
      Apsi[cb] *= ct;
      Apsi[cb] += Ap[cb] * st;
    }

    /*  Ap = g[k] = Apsi[k] - mu[k] Psi[k] */
    for(cb = 0; cb < NLayers; ++cb)
    {
      Ap[cb] = Apsi[cb] - psi[cb]*lambda;
    }

    /*  g2  =  |g[k]|**2 = |Ap|**2 */
    s1 = g2;			/* Now s1 is |g[k-1]|^2 */
    g2 = norm2(Ap);

    /*+  */
    /*  IF |g[k]| <= max(RsdR_a, RsdR_r |mu[k]|)
	&& |del_lam| <= Rsdlam*|lambda|  THEN RETURN; */
    rsd = rsdr_sq * lambda * lambda;
    rsd = max(rsda_sq, rsd);
    st = Rsdlam * fabs(mu);	/* old value of st is no longer needed */

/*    if ( Real(g2) <= rsd && del_lam <= st ) */
    if ( (! Kalk_Sim && k >= N_min && (fabs(mu) < Cutl_zero ||
				       (Real(g2) <=  rsd && del_lam <= st ) ) ) ||
	 (Kalk_Sim && ( del_lam <= st || fabs(mu) < Cutl_zero ||
			(k >= N_min && (Real(g2) < g2_0 || k >= N_max ) ) ) ) )
    {
      QDPIO::cout << "Converged at iter " << k << ", lambda = " << lambda
		  << ", del_lam = " << del_lam << endl
		  << "                      rsd = " << rsd 
		  << ", g2 = " << g2 << ", g2_0 = " << g2_0 << endl;
      n_count = k;

      /* Renormalize and recompute lambda */
      /*  Project to orthogonal subspace  */
      GramSchm (psi, 1, psi_all, NN_eig, NLayers);

      /* Normalize */
      d = norm2(psi);

      d = sqrt(d);
      ct = Real(1) / Real(d);
      for(cb = 0; cb < NLayers; ++cb)
	psi[cb] *= ct;
      d -= one;
      QDPIO::cout << "Deviation at convergence: " << d << endl;

      /*  Apsi  :=  A . Psi  */
      A(Apsi, psi, PLUS);

      /*  mu  := < Psi | A Psi >  */
      s1 = mu;
      mu = innerProductReal(psi, Apsi);
      lambda = mu;
      QDPIO::cout << "Mu-s at convergence: old " << s1 << " vs. NN_eig " << mu << endl;

      /*  Copy psi back into psi_all(N_eig-1)  */
      psi_all[NN_eig] = psi;

      END_CODE("Ritz");
      return;
    }

    /* b = beta[k] = cos(theta) |g[k]|^2 / |g[k-1]|^2 */
    b = ct * Real(g2/s1);
    ct *= 0.05 * b;
    d = sqrt(g2);
    ct /= Real(p2*d);
    if( ct > 1.0 )
    {
      /* Restart: p[k] = g[k] = Ap */
      QDPIO::cout << "Restart at iter " << k << " since beta = " << b << endl;
      p = Ap;
    }
    else
    {
      /* xp = < Psi[k] | p[k-1] > */
      xp = innerProduct(psi, p);

      /* p[k] = g[k] + b (p[k-1] - xp psi[k]) */
      for(cb = 0; cb < NLayers; ++cb)
      {
        p[cb] -= psi[cb] * xp;
	p[cb] *= b;
	p[cb] += Ap[cb];
      }
    }

    if( k % n_renorm == 0 )
    {
      /* Renormalize, and re-orthogonalize */
      /*  Project to orthogonal subspace  */
      GramSchm (psi, 1, psi_all, NN_eig, NLayers);

      /* Normalize */
      d = norm2(psi);

      d = sqrt(d);
      ct = Real(1) / Real(d);
      for(cb = 0; cb < NLayers; ++cb)
	psi[cb] *= ct;
      d -= one;
#if 1
      if (ProjApsiP)
	QDPIO::cout << "Deviation at iter " << k << " : " << d << endl;
#endif

      /*  Apsi  :=  A . Psi  */
      A(Apsi, psi, PLUS);

      /*  Project to orthogonal subspace, if wanted  */
      /** Should not be necessary, following Kalkreuter-Simma **/
      if (ProjApsiP)
	GramSchm (Apsi, 1, psi_all, NN_eig, NLayers);

      /*  mu  := < Psi | A Psi >  */
      mu = innerProductReal(psi, Apsi);
#if 1
      if (ProjApsiP)
	QDPIO::cout << "Mu-s at iter " << k ": old " << s1 << " vs. NN_eig " << mu << endl;
#endif

      /*  g[k] = Ap = ( A - mu ) Psi  */
      lambda = Real(mu);
      for(cb = 0; cb < NLayers; ++cb)
      {
	Ap[cb] = Apsi[cb];
	Ap[cb] -= psi[cb] * lambda;
      }

      /*  g2  =  |g[k]|**2 = |Ap|**2 */
      g2 = norm2(Ap);

      /*  Project p[k] to orthogonal subspace  */
      GramSchm (p, 1, psi_all, NN_eig, NLayers);

      /*  Make p[k] orthogonal to Psi[k]  */
      GramSchm (p, 1, psi, 1, NLayers);

      /*  Make < g[k] | p[k] > = |g[k]|**2: p[k] = p_old[k] + xp g[k],
	  xp = (g2 - < g | p_old >) / g2; g[k] = Ap */
      p2 = 0;
      ct = Real(1) / Real(g2);
      xp = ct*(cmplx(g2,p2) - innerProduct(Ap, p));

      for(cb = 0; cb < NLayers; ++cb)
	p[cb] += Ap[cb] * xp;
    }
    else if (ProjApsiP && NN_eig > 0)
    {
      /*  Project psi and p to orthogonal subspace  */
      GramSchm (psi, 1, psi_all, NN_eig, NLayers);
      GramSchm (p, 1, psi_all, NN_eig, NLayers);
    }

    /*  p2  =  |p[k]|**2 */
    p2 = 0;
    for(cb = 0; cb < NLayers; ++cb)
      p2 += norm2(p[cb]);	                /* 2 Nc Ns  flops */

#if 1
    /* if( (2*k)%n_renorm == 0 ) */
    if (ProjApsiP)
      QDPIO::cout << "At iter " << k << ", lambda = " << lambda 
		  << ", del_lam = " << del_lam
		  << ", g2 = " << g2 << endl;
#endif
  }

  /*  Copy psi back into psi_all(N_eig-1)  */
  for(cb = 0; cb < NLayers; ++cb)
    psi_all[NN_eig][cb] = psi[cb];

  n_count = MaxCG;
  QDPIO::cerr << "too many CG/Ritz iterations: n_count=" << n_count
	      << ", rsd=" << rsd << ", g2" << g2 << ", p2" << p2
	      << ", lambda" << lambda << endl;
  QDP_abort(1);
  END_CODE("Ritz");
}
