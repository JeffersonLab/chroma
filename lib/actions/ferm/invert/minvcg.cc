// $Id: minvcg.cc,v 1.1 2003-04-08 21:36:13 edwards Exp $

/*! \file
 *  \brief Multishift Conjugate-Gradient algorithm for a Linear Operator
 */

#include "chroma.h"
#include "actions/ferm/invert/minvcg.h"

using namespace QDP;

//! Multishift Conjugate-Gradient (CG1) algorithm for a  Linear Operator
/*! \ingroup invert
 *
 * This subroutine uses the Conjugate Gradient (CG) algorithm to find
 * the solution of the set of linear equations
 * Method used is described in  Jegerlehner, hep-lat/9708029

 * We are searching in a subspace orthogonal to the eigenvectors EigVec
 * of A. The source chi is assumed to already be orthogonal!

 *   	    Chi  =  A . Psi

 * Algorithm:

 *  Psi[0] :=  0;      	                      Zeroed
 *  r[0]   :=  Chi;                           Initial residual
 *  p[1]   :=  Chi ;	       	       	      Initial direction
 *  b[0]   := |r[0]|**2 / <p[0],Ap[0]> ;
 *  z[0]   := 1 / (1 - (m - m(0))*b) 
 *  bs[0]  := b[0] * z[0]  
 *  r[1] += b[k] A . p[0] ; 	       	      New residual
 *  Psi[1] = - b[k] p[k] ;   	       	      Starting solution vector
 *  IF |r[0]| <= RsdCG |Chi| THEN RETURN;        Converged?
 *  FOR k FROM 1 TO MaxCG DO    	       	       CG iterations
 *      a[k] := |r[k]|**2 / |r[k-1]|**2 ;
 *      p[k] := r[k] + a[k] p[k-1];   	       New direction
 *      b[k+1] := |r[k]|**2 / <p[k],Ap[k]> ;
 *      r[k+1] += b[k+1] A . p[k] ; 	       	       New residual
 *      Psi[k+1] -= b[k+1] p[k] ;   	       	       New solution vector
 *      IF |[k+1]| <= RsdCG |Chi| THEN RETURN;    Converged?

 * Arguments:

 *  A	        Hermitian linear operator      (Read)
 *  Chi	        Source   	               (Read)
 *  Psi	        array of solutions    	       (Write)
 *  Chi_Norm    |Chi|	    	    	       (Read)
 *  mass        shifts of form  A + mass       (Read)
 *  Nmass       Number of shifts               (Read)
 *  isz         which mass in list is smallest (Read)
 *  RsdCG       residual accuracy              (Read)
 *  N_Count     Number of CG iteration	       (Write)
 *  EigVec      Eigenvectors, to which we want orthogonality	(Read)
 *  NEig        Number of eigenvectors         (Read)

  * Local Variables:

 *  p   	       Direction vector
 *  r   	       Residual vector
 *  cp  	       | r[k] |**2
 *  c   	       | r[k-1] |**2
 *  k   	       CG iteration counter
 *  a   	       a[k]
 *  b   	       b[k+1]
 *  d   	       < p[k], A.p[k] >
 *  Ap  	       Temporary for  M.p

 *  MaxCG       Maximum number of CG iterations allowed

 * Subroutines:
 *  A	       Apply matrix hermitian A to vector 
 */

void MInvCG(const LinearOperator& A, const LatticeFermion& chi, LatticeFermion& psi,
	    const Double& chi_norm, 
	    const multi1d<Real>& mass, int Nmass, int isz, 
	    const multi1d<Real>& RsdCG, int n_count&,
	    const multi1d<LatticeFermion>& EigVec, int NEig)

/*  Vec_orth    Vectors to which we want orthogonality	(Read) */
/*  NVec        Number of vectors Vec_orth     (Read) */
{ /* Local Variables */
  multi1d<LatticeFermion> p(Nmass);
  LatticeFermion r;
  LatticeFermion Ap;
  LatticeFermion ltmp;
  Real a;
  Real as;
  Real b;
  Real bp;
  multi1d<Real> bs(Nmass);
  multi2d<Real> z(2, Nmass);
  Real css;
  Real ztmp;
  Double c;
  Double cs;
  Double d;
  Double cp;
  Boolean convP;
  multi1d<Boolean> convsP(Nmass);
  int iz;
  int k;
  int s;
  multi1d<Real> rsd_sq(Nmass);
  multi1d<Real> rsdcg_sq(Nmass);
  
  START_CODE("MInvCGm");
  
  /* If exactly 0 norm, then solution must be 0 (for pos. def. operator) */
  if (chi_norm == 0.0)
  {
    n_count = 0;
    psi = 0;
    END_CODE("MInvCGm");
    return;
  }

    
  cp = chi_norm * chi_norm;
  for(s = 0; s < Nmass; ++s)
  {
    rsdcg_sq[s] = RsdCG[s] * RsdCG[s];
    rsd_sq[s] = Real(cp) * rsdcg_sq[s];
  }

              
  /* Psi[0] := 0; */
  /* r[0] := p[0] := Chi */
  r = chi;
  for(s = 0; s < Nmass; ++s)
    p[s] = chi;

  /*  b[0] := - | r[0] |**2 / < p[0], Ap[0] > ; */
  /*  First compute  d  =  < p, A.p >  */
  /*  Ap = A . p  */
  ltmp = p[isz];
  A (ltmp, Ap, 0);

#if 1
  /* Project out eigenvectors */
  GramSchm (Ap, 1, EigVec, NEig);
#endif

  Ap += p[isz] * mass[isz];

  /*  d =  < p, A.p >  */
  d = sum(real(trace(adj(Ap) * ltmp)));                        /* 2 Nc Ns  flops */
  
  b = -Real(cp/d);

  /* Compute the shifted bs and z */
  z[isz][0] = 1;
  z[isz][1] = 1;
  bs[isz] = b;
  iz = 1;
  for(s = 0; s < Nmass; ++s)
  {
    if (s == isz) continue;

    z[s][1-iz] = 1;
    z[iz][s] = Real(1) / (1 - (mass[s]-mass[isz])*b);
    bs[s] = b * z[s][iz];
  }

  /*  r[1] += b[0] A . p[0]; */
  r += Ap * b;	   /* 2 Nc Ns  flops */

  /*  Psi[1] -= b[0] p[0] = - b[0] chi; */
  for(s = 0; s < Nmass; ++s)
    psi[s] = -(chi * bs[s]); /* 2 Nc Ns  flops */

  /*  c = |r[1]|^2 */
  c = norm2(r);   	       	   /* 2 Nc Ns  flops */

  /*  Check convergence of first solution */
  convsP = false;

  css = Real(c);
  convP = css < rsd_sq[isz];

  PRINTF("MInvCG: k = %d  r = %g\n",0,sqrt(c));

  /*  FOR k FROM 1 TO MaxCG DO */
  /*  IF |psi[k+1] - psi[k]| <= RsdCG |psi[k+1]| THEN RETURN; */
  for(k = 1; k <= MaxCG && ! convP ; ++k)
  {
    /*  a[k+1] := |r[k]|**2 / |r[k-1]|**2 ; */
    a = Real(c/cp);

    /*  p[k+1] := r[k+1] + a[k+1] p[k]; */
    /*  Compute the shifted as */
    /*  ps[k+1] := zs[k+1] r[k+1] + a[k+1] ps[k]; */
    for(s = 0; s < Nmass; ++s)
    {
      if (convsP[s]) continue;

      if (s == isz)
      {
	p[s] *= a;	/* Nc Ns  flops */
	p[s] += r;	/* Nc Ns  flops */
	}
      }
      else
      {
	as = a * z[iz][s]*bs[s] / (z[1-iz][s]*b);

	p[s] *= as;	/* Nc Ns  flops */
	p[s] += r * z[s][iz];	/* Nc Ns  flops */
      }
    }

#if 0
    /* Project out eigenvectors */
    if (k % 10 == 0)
      GramSchm (p, Nmass, EigVec, NEig);
#else
#if 1
    /* Project out eigenvectors */
    if (k % 10 == 0)
    {
      ltmp = p[isz];
      GramSchm (ltmp, 1, EigVec, NEig);
      p[isz] = ltmp;
    }
#endif
#endif

    /*  cp  =  | r[k] |**2 */
    cp = c;

    /*  b[k] := | r[k] |**2 / < p[k], Ap[k] > ; */
    /*  First compute  d  =  < p, A.p >  */
    /*  Ap = A . p  */
    ltmp = p[isz];
    A (ltmp, Ap, 0);

#if 1
    /* Project out eigenvectors */
    if (k % 10 == 0)
      GramSchm (Ap, 1, EigVec, NEig);
#endif

    Ap += p[isz] * mass[isz];

    /*  d =  < p, A.p >  */
    d = sum(real(trace(adj[Ap] * ltmp)));       /* 2 Nc Ns  flops */
    
    bp = b;
    b = -Real(cp/d);

    /* Compute the shifted bs and z */
    bs[isz] = b;
    iz = 1 - iz;
    for(s = 0; s < Nmass; ++s)
    {
      if (s == isz || VALUE(convsP(s))) continue;

      ztmp = VALUE(z(1-iz,s))*VALUE(z(iz,s))*bp / 
	(b*a*(VALUE(z(iz,s)) - VALUE(z(1-iz,s))) + 
	 VALUE(z(iz,s))*bp*(1 - (VALUE(mass(s)) - VALUE(mass(isz)))*b));
      VALUE(bs(s)) = b * ztmp / VALUE(z(1-iz,s));
      z[s][iz] = ztmp;
    }

    /*  r[k+1] += b[k] A . p[k] ; */
    r += Ap * b;	        /* 2 Nc Ns  flops */

#if 1
    /* Project out eigenvectors */
    if (k % 10 == 0)
      GramSchm (r, 1, EigVec, NEig);
#endif

    /*  Psi[k+1] -= b[k] p[k] ; */
    for(s = 0; s < Nmass; ++s)
    {
      if (VALUE(convsP(s))) continue;

      psi[s] -= p[s] * bs[s];	/* 2 Nc Ns  flops */
    }

    /*  c  =  | r[k] |**2 */
    c = norm2(r);	                /* 2 Nc Ns  flops */

    /*    IF |psi[k+1] - psi[k]| <= RsdCG |psi[k+1]| THEN RETURN; */
    /* or IF |r[k+1]| <= RsdCG |chi| THEN RETURN; */
    for(s = 0; s < Nmass; ++s)
    {
      if (VALUE(convsP(s))) continue;

      /* Convergence methods */
#if 1
      /* Check norm of shifted residuals */
      css = Real(c) * VALUE(z(iz,s))*VALUE(z(iz,s));

      convsP[s] = css < rsd_sq[s];

      if (s == isz)
	PRINTF("MInvCG: k = %d  s = %d  r = %g\n",k,s,sqrt(css));
#else
      /* Check relative error of solution */
      cs = norm2(p[s]);         	        /* 2 Nc Ns  flops */
      d = norm2(psi[s]);         	/* 2 Nc Ns  flops */

      cs *= VALUE(bs(s))*VALUE(bs(s));
      d *= VALUE(rsdcg_sq(s));
      convsP[s] = cs < d;

      if (s == isz)
	PRINTF("MInvCG: k = %d  s = %d  r = %g  d = %g\n",k,s,sqrt(cs),sqrt(d));
#endif
    }
    /* Final convergence is determined by smallest shift */
    convP = convsP[isz];

    n_count = k;
  }

#if 1
  GramSchm (psi, Nmass, EigVec, NEig);
#endif

#if 0
  /* HACK **/
  for(s = 0; s < Nmass; ++s)
  {
    ltmp = psi[s];
    A (ltmp, Ap, 0);
    Ap += psi[s] * mass[s];

    Ap -= chi;

    c = norm2(Ap);	                /* 2 Nc Ns  flops */

    PRINTF("MInvCG (conv): s = %d  r = %g  m = %g\n",s,sqrt(c),VALUE(mass(s)));
  }
  PRINTF("\n");
  /* end */
#endif

  if (n_count == MaxCG)
    QDP_error_exit("too many CG iterations", n_count, rsd_sq[isz], d, c, cp, a, b);
  END_CODE("subroutine");;
  return;
}
