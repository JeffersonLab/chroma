// $Id: lovddag_w.cc,v 1.3 2003-11-09 22:35:19 edwards Exp $
/*! \file
 *  \brief Internal Overlap-pole operator for D^dag.D
 */

#include "chromabase.h"
#include "actions/ferm/linop/lovlapms_w.h"

using namespace QDP;

//! Apply the operator onto a chiral source vector
/*! \ingroup linop
 *
 * The operator acts on the entire lattice
 *
 * \param psi 	  Pseudofermion field     	       (Read)
 * \param isign   Flag ( PLUS | MINUS )   	       (Read)
 */
LatticeFermion lovddag::operator() (const LatticeFermion& psi, enum PlusMinus isign) const
{
  LatticeFermion chi;

  Real m_q;
  int numroot;
  Real constP;
  multi1d<Real> resP(numroot);
  multi1d<Real> rootQ(numroot);
  Real RsdCG;
  int NEig;
  multi1d<Real> EigValFunc(NEig);
  multi2d<LatticeFermion> EigVec(NEig);    /* Nsubl and Ncb better agree */
  int ichiral;                    /* source is up (+1) or down (-1) */

  Real mass;

  LatticeFermion Ap;
  LatticeFermion ltmp;
  multi1d<LatticeFermion> p(numroot);
  LatticeFermion r;
  LatticeFermion tmp1;

  Real a;
  Real as;
  Real b;
  Real bp;
  Real ztmp;
  Real z0;
  Real z1;   
  Double cp;
  Double c;
  Double d;
  multi1d<Real> bs(numroot);
  multi2d<Real> z(2, numroot);
  Real rsd_sq;
  Real rsdcg_sq;
  Double chi_sq_new;
  Double chi_sq_diff;

  multi1d<Boolean> convsP(numroot);
  Boolean convP;
  Boolean btmp;

  int iz;
  int n_count;
  int n;
  int k;
  int i;
  int s;
  int isz;
  int cb;

  START_CODE("lovddag");

  int G5 = Ns*Ns - 1;

  /* Extract arguments */
  UNPACK_LINEAR_OPERATOR(A, B, ichiral);

  UNPACK_LINEAR_OPERATOR(B, MdagM, M, m_q, 
			 numroot, constP, resP, rootQ, 
			 EigVec, EigValFunc, NEig, RsdCG);

  Real mass = (1 + m_q) / (1 - m_q);

        
  LatticeFermion psi_proj = psi;

  /* Initially  chi := (1/mass + mass) * psi  */
  chi = (mass + 1 / mass) * psi;

  /* Project out eigenvectors of source if desired */
  /* chi  +=  eps(lambda) * (gamma_5 + ichiral) * EigVec  */
  /* Usually "func(.)" is eps(.); it is precomputed in EigValFunc. */
  if (NEig > 0)
  {
    Complex cconsts;

    for(i = 0; i < NEig; ++i)
    {
      cconsts = innerProduct(EigVec[i], psi_proj);
      psi_proj -= EigVec[i] * cconsts;

      /* Construct  tmp1 = (gamma_5 + ichiral) * EigVec */
      tmp1 = Gamma(G5) * EigVec[i];

      if (ichiral == PLUS)
      {
	tmp1 += EigVec[i];
      }
      else
      {
	tmp1 -= EigVec[i];
      }

      /* chi += "eps(lambda)" * (gamma_5 + ichiral) * EigVec  */
      cconsts *= EigValFunc[i];
      chi += tmp1 * cconsts;
    }
    
  }

  /* First part of application of D.psi */
  /* tmp1 <- H * psi */
  /*      <- gamma_5 * M * psi */
  ltmp = M(psi_proj, PLUS);
  tmp1 = Gamma(G5) * ltmp;

  c = norm2(tmp1);

  /* If exactly 0 norm, then solution must be 0 (for pos. def. operator) */
  if (c == 0.0)
  {
    chi = 0;
            
    END_CODE("lovddag");
    return;
  }

  /*  chi  +=  constP*(gamma_5 + ichiral) * H * Psi  */
  switch (ichiral)
  {
  case PLUS:
    ltmp += tmp1;
    break;

  case MINUS:
    ltmp -= tmp1;
    break;

  default:
    QDP_error_exit("unsupported chirality", ichiral);
  }

  chi += ltmp * constP;

  /********************************************************************/
  /* Solve  (MdagM + rootQ_n) chi_n = H * tmp1 */
  rsdcg_sq = RsdCG * RsdCG;
  rsd_sq = c * rsdcg_sq;

  /* By default (could change), rootQ(isz) is considered the smallest shift */
  isz = numroot-1;
          
  /* r[0] := p[0] := tmp1 */
  r = tmp1;
  for(s = 0; s < numroot; ++s)
    p[s] = r;

  z = 1;
  a = 0;
  b = 1;
  convsP = false;
  convP = false;
  iz = 1;

  for(k = 0; k <= MaxCG && ! convP ; ++k)
  {
    n_count = k;

    /*  b[k] := | r[k] |**2 / < p[k], Ap[k] > ; */
    /*  First compute  d  =  < p, A.p >  */
    /*  Ap = A . p  */
    ltmp = p[isz];
    Ap = MdagM(ltmp, PLUS);
    
    Ap += p[isz] * rootQ[isz];

#if 1
    /* Project out eigenvectors */
    if (k % 10 == 0)
      GramSchm (Ap, 1, EigVec, NEig);
#endif

    /*  d =  < p, A.p >  */
    d = innerProduct(Ap, ltmp);   /* 2 Nc Ns  flops */
    
    bp = b;
    b = -Real(c/d);

    /* Compute the shifted bs and z */
    bs[isz] = b;
    iz = 1 - iz;
    for(s = 0; s < numroot; ++s)
    {
      /* Do this to avoid mitsmp compiler bug !! */
      z0 = z[s][1-iz];
      z1 = z[s][iz];

      if (s == isz || convsP[s]) continue;

      z[iz][s] = z0*z1*bp / 
	(b*a*(z1 - z0) + z1*bp*(1 - (rootQ[s] - rootQ[isz])*b));
      bs[s] = b * z[iz][s] / z0;
    }

    /*  r[k+1] += b[k] A . p[k] ; */
    r += Ap * b;	        /* 2 Nc Ns  flops */

#if 1
    /* Project out eigenvectors */
    if (k % 10 == 0)
      GramSchm (r, 1, EigVec, NEig, Ncb);
#endif
    
    /*  soln[k+1] -= b[k] p[k] ; */
    ztmp = b * resP[isz];
    ltmp = p[isz] * ztmp;	/* 2 Nc Ns  flops */

    for(s = 0; s < numroot; ++s)
    {
      if (s == isz || convsP[s]) continue;

      ltmp += p[s] * bs[s] * resP[s];	/* 2 Nc Ns  flops */
    }
    
    /*  diff :=  (gamma_5 + ichiral) * tmp1  */
    tmp1 = Gamma(G5) * ltmp;
    if (ichiral == PLUS)
    {
      tmp1 += ltmp;
    }
    else
    {
      tmp1 -= ltmp;
    }

    /* chi += diff . The minus comes from the CG method. */
    chi -= tmp1;                   /* 2 Nc Ns  flops */


    /*  cp  =  | r[k] |**2 */
    cp = c;

    /*  c  =  | r[k] |**2 */
    c = norm2(r);	                /* 2 Nc Ns  flops */

    /*  a[k+1] := |r[k]|**2 / |r[k-1]|**2 ; */
    a = Real(c/cp);

    /*  p[k+1] := r[k+1] + a[k+1] p[k]; */
    /*  Compute the shifted as */
    /*  ps[k+1] := zs[k+1] r[k+1] + a[k+1] ps[k]; */
    for(s = 0; s < numroot; ++s)
    {
      if (convsP[s]) continue;

      if (s == isz)
      {
	p[s] *= a;	/* Nc Ns  flops */
	p[s] += r;	/* Nc Ns  flops */
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
      GramSchm (p, numroot, EigVec, NEig);
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

    /*  IF |r[k]| <= RsdCG |chi| THEN RETURN; */
    convP = true;
    for(s = 0; s < numroot; ++s)
    {
      if (convsP[s]) continue;

      ztmp = Real(c) * z[iz][s]*z[iz][s];

      btmp = ztmp < rsd_sq;
      convP = convP & btmp;
      convsP[s] = btmp;

#if 0
      if (k % 10 == 0)
	QDP_info("Lovddag: k = %d  s = %d  r = %g  rsd = %g",k,s,sqrt(ztmp),sqrt(rsd_sq));
#endif
    }


    /* soln = ( \sum_i { res_i / (z^2 + rootQ_i) } ) * H * psi  */
    /* Check for convergence of chi */
    /* Only converge if chi is converged. If vectors converge first, then error */
    if (k > 0 && ! convP)
    {
      chi_sq_new = 0;
      chi_sq_diff = 0;

      chi_sq_new += norm2(chi);
      chi_sq_diff += norm2(tmp1);

      chi_sq_new *= Double(rsdcg_sq);
      btmp = chi_sq_diff < chi_sq_new;

#if 0
      if (k % 10 == 0)
	QDP_info("Lovddag (chi): k = %d  diff = %g  new = %g  rsdcg = %g",
	       k,sqrt(chi_sq_diff),sqrt(chi_sq_new),sqrt(rsdcg_sq));
#endif

      if (! btmp && convP)
	QDP_error_exit("vectors converged but not final chi");

      convP = btmp;
    }
  }

  /* Rescale to the correct normalization */
  /*  chi <-  (1/4)*(1-m_q^2) * chi  */
  chi *= 0.25 * (1 - m_q*m_q);   /* 2 Nc Ns  flops */

#if 0
  QDP_info("Lovddag (chi): k = %d  diff = %g  new = %g  rsdcg = %g",
	 n_count,sqrt(chi_sq_diff),sqrt(chi_sq_new),sqrt(rsdcg_sq));
#endif

#if 1
  QDP_info("ovddag(%d): %d", ichiral, n_count);
#endif

  END_CODE("lovddag");

  return chi;
}

