// $Id: lovlapms_w.cc,v 1.6 2003-12-15 17:52:51 bjoo Exp $
/*! \file
 *  \brief Overlap-pole operator
 */
#include <math.h>
#include "chromabase.h"
#include "actions/ferm/linop/lovlapms_w.h"

using namespace QDP;

//! Apply the operator onto a source vector
/*! \ingroup linop
 *
 * The operator acts on the entire lattice
 *
 * \param psi 	  Pseudofermion field     	       (Read)
 * \param isign   Flag ( PLUS | MINUS )   	       (Read)
 */
void lovlapms::operator() (LatticeFermion& chi, const LatticeFermion& psi, 
			   enum PlusMinus isign) const
{
  START_CODE("lovlapms");
  LatticeFermion Ap,tmp1, tmp2;
  LatticeFermion r;
  multi1d<LatticeFermion> p(numroot);

  Real a;
  Real as;
  Real b;
  Real bp;
  Real ztmp;
  Double cp;
  Double d;
  Real z0;
  Real z1;   
  multi1d<Real> bs(numroot);
  multi2d<Real> z(2, numroot);
  Double chi_sq_new;
  Double chi_sq_diff;


  int G5 = Ns*Ns - 1;

  Real mass = Real(1 + m_q) / Real(1 - m_q);
  


  switch (isign)
  {
  case PLUS:
    /*  Non-Dagger: psi is source and tmp1 */
    /*  chi  :=  gamma_5 * (gamma_5 * mass + eps(H)) * Psi  */
    tmp1 = psi;
    break;

  case MINUS:
    /*  Dagger: apply gamma_5 to source psi to make tmp1 */
    /*  chi  :=  (mass + eps(H) * gamma_5) * Psi  */
    tmp1 = Gamma(G5) * psi;
    break;

  default:
    QDP_error_exit("unknown isign value", isign);
  }

  chi = zero;

  /* Project out eigenvectors of source if desired */
  /* chi  +=  func(lambda) * EigVec * <EigVec, psi>  */
  /* Usually "func(.)" is eps(.); it is precomputed in EigValFunc. */

  QDPIO::cout << "NEig = " << NEig << endl;

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

  /* tmp1 <- H * Projected psi */
  /*      <- gamma_5 * M * psi */

  M(tmp2, tmp1, PLUS);
  tmp1 = Gamma(G5) * tmp2;
  


  Double c = norm2(tmp1);

  /* If exactly 0 norm, then solution must be 0 (for pos. def. operator) */
  if (toBool(c == 0.0))
  {
    chi = zero;

    END_CODE("lovlapms");
    return;
  }


  QDP_info("Doing multi-mass solve\n");
  /********************************************************************/
  /* Solve  (MdagM + rootQ_n) chi_n = H * tmp1 */

  Real rsdcg_sq = RsdCG * RsdCG;
  Real rsd_sq = c * rsdcg_sq;


  cout << "Rsd_sq = " << rsd_sq << endl;

  /* By default (could change), rootQ(isz) is considered the smallest shift */
  int isz = numroot-1;

  cout << "isz = " << isz << endl;
  /* chi[0] := mass*psi + c0*H*tmp1 + Eigvecs; */
  if (isign == PLUS)
  {
    /*  chi  :=  gamma_5 * (gamma_5 * mass + eps(H)) * Psi  */
    /* Final mult by gamma_5 is at end */
    chi += Gamma(G5) * psi * mass;
  }
  else
  {
    /*  chi  :=  (mass + eps(H) * gamma_5) . Psi  */
    chi += psi * mass;
  }

  chi += tmp1 * constP;

  /* r[0] := p[0] := tmp1 */
  r = tmp1;

  for(int s = 0; s < numroot; ++s)
    p[s] = tmp1;

  multi1d<bool> convsP(numroot);
  convsP = false;
  bool convP = false;
  int iz = 1;



  z = 1;
  a = 0;
  b = 1;

  int n_count;
  for(int k = 0; k <= MaxCG && ! convP ; ++k)
  {
    n_count = k;

    /*  b[k] := | r[k] |**2 / < p[k], Ap[k] > ; */
    /*  First compute  d  =  < p, A.p >  */
    /*  Ap = A . p  */


    /* This bit computes 
       Ap = [  M^dag M + rootQ(isz)  ] p_isz
    */

    MdagM(Ap, p[isz], PLUS);
    Ap += p[isz] * rootQ[isz];

#if 000000
    /* Project out eigenvectors */
    if (k % 2 == 0)
      GramSchm (Ap, 1, EigVec, NEig);
#endif

    /*  d =  < p, A.p >  */
    d = innerProductReal(p[isz], Ap);                       /* 2 Nc Ns  flops */
    
    bp = b;
    b = -Real(c/d);

    /* Compute the shifted bs and z */
    bs[isz] = b;
    iz = 1 - iz;
    for(int s = 0; s < numroot; ++s)
    {
      /* Do this to avoid mitsmp compiler bug !! */
      z0 = z[s][1-iz];
      z1 = z[s][iz];

      if (s == isz || convsP[s]) continue;

      z[s][iz] = z0*z1*bp / 
	(b*a*(z1 - z0) + z1*bp*(1 - (rootQ[s] - rootQ[isz])*b));
      bs[s] = b * z[s][iz] / z0;
    }

    /*  r[k+1] += b[k] A . p[k] ; */
    r += Ap * b;	        /* 2 Nc Ns  flops */

#if 000000
    /* Project out eigenvectors */
    if (k % 2 == 0)
      GramSchm (r, 1, EigVec, NEig);
#endif
    
    /*  chi[k+1] -= b[k] p[k] ; */
    ztmp = b * resP[isz];
    tmp1 = p[isz] * ztmp;	/* 2 Nc Ns  flops */

    for(int s = 0; s < numroot; ++s)
    {
      if (s == isz || convsP[s]) continue;

      ztmp = bs[s] * resP[s];
      tmp1 += p[s] * ztmp;	/* 2 Nc Ns  flops */
    }

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
    for(int s = 0; s < numroot; ++s)
    {
      if (convsP[s]) continue;

      if (s == isz)
      {
	p[s] *= a;	        /* Nc Ns  flops */
	p[s] += r;		/* Nc Ns  flops */
      }
      else
      {
	as = a * z[s][iz]*bs[s] / (z[s][1-iz]*b);

	p[s] *= as;	        /* Nc Ns  flops */
	p[s] += r * z[s][iz];	/* Nc Ns  flops */
      }
    }

#if 0
    /* Project out eigenvectors */
    if (k % 10 == 0)
      GramSchm (p, numroot, EigVec, NEig);
#endif
    
    /*  IF |r[k]| <= RsdCG |chi| THEN RETURN; */
    convP = true;
    for(int s = 0; s < numroot; ++s)
    {
      if (convsP[s]) continue;

      ztmp = Real(c) * z[s][iz]*z[s][iz];

      bool btmp = toBool(ztmp < rsd_sq);
      convP = convP & btmp;
      convsP[s] = btmp;
    }


    /* eps(H) psi = ( c_0 + \sum_i { res_i / (z^2 + rootQ_i) } ) * H * psi  */
    /* Norm of diff of new and old soln unchanged by a potential gamma_5 */
    /* Check for convergence of chi */
    /* Only converge if chi is converged. If vectors converge first, then error */
    if (k > 0 && ! convP)
    {
      chi_sq_new = rsdcg_sq * norm2(chi);
      chi_sq_diff = norm2(tmp1);      /* the diff of old and new soln */

      bool btmp = toBool(chi_sq_diff < chi_sq_new);


      QDPIO::cout << "Lovlapms: k =" << k <<" c= " << Real(c) << endl;


      if (! btmp && convP)
	QDP_error_exit("vectors converged but not final chi");

/*      convP = convP & btmp; */
      convP = btmp;
    }
  }

  if (isign == PLUS)
  {
    /*  chi  :=  gamma_5 * (gamma_5 * mass + eps(H)) * Psi  */
    tmp1 = Gamma(G5) * chi;
    chi = tmp1;
  }

  /* Rescale to the correct normalization */
  chi *= 0.5 * (1 - m_q);

  QDP_info("ovlapms: %d", n_count);
                
  END_CODE("lovlapms");
}

