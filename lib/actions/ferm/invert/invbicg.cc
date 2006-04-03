// $Id: invbicg.cc,v 3.0 2006-04-03 04:58:49 edwards Exp $

// THIS ROUTINE IS NOT FUNCTIONAL YET - SIMPLE CLEANUP EDITING IS REQUIRED

namespace Chroma {

/* INVBICG */

/* This subroutine uses the stabilized BiConjugate Gradient algorithm to */
/* determine the solution of the set of linear equations */

/*   	    Chi  =  A . Psi */

/* Algorithm: */

/*  Psi[0]                                      Argument */
/*  r[0]    :=  Chi  -  A . Psi[0] ;    	Initial residual */
/*  rv      :=  r[0] ;                          for second Krylov space */
/*  p[1]    :=  r[0] ;                          initial direction */
/*  IF |r[0]| <= RsdBiCG |Chi| THEN RETURN;     Converged? */
/*  FOR k FROM 1 TO MaxCG DO                    BiCG iterations */
/*      mp[k]   := A.p[k] ; */
/*      a[k]    := <rv,r[k-1]> / <rv,mp[k]> ; */
/*      s[k]    := r[k-1] - a[k] mp[k] ; */
/*      ms[k]      := A.s[k] ; */
/*      omega[k]:= <ms[k],s[k]> / <ms[k],ms[k]> ; */
/*      Psi[k]  += omega[k] s[k] + a[k] p[k] ;  New solution vector */
/*      r[k]    := s[k] - omega[k] ms[k] ;      New residual */
/*      IF |r[k]| <= RsdBiCG |Chi| THEN RETURN; Converged? */
/*      b[k+1]  := (<rv,r[k]> / <rv,r[k-1]>) * (a[k]/omega[k]) ; */
/*      p[k+1]  := r[k] + b[k+1] * (p[k] - omega[k] mp[k]) */

/* Arguments: */

/*  A	       Linear Operator                 (Read) */
/*  Chi	       Pseudofermion Source   	       (Read) */
/*  Psi	       Solution    	    	       (Modify) */
/*  Chi_Norm   |Chi|	    	    	       (Read) */
/*  RsdBiCG    BiCG residue                    (Read) */
/*  Ncb        Number of checkerboards         (Write) */
/*  N_Count    Number of BiCG iteration	       (Write) */

/* Local Variables: */

/*  r          Residual vector */
/*  rv         For second Krylov space */
/*  p          Direction vector */
/*  cp         | r[k] |**2 */
/*  c          | r[k-1] |**2 */
/*  rvr        < rv, r[k-1] > */
/*  k          BiCG iteration counter */
/*  a          a[k] */
/*  omega      omega[k] */
/*  ms2        | Ms[k] |**2 */
/*  d          Temporary for < rv, A.p[k] >, < A.s[k], s[k] > and rvr_old */
/*  Mp 	       Temporary for  M.p */
/*  s  	       Temporary */
/*  Ms 	       Temporary for  M.s */

/* Global Variables: */

/*  MaxCG       Maximum number of BiCG iterations allowed */

/* Subroutines: */

/*  A    Apply matrix linear operator to vector */

/* Operations: */

/*  ???? */

include(types.mh)

SUBROUTINE(InvBiCGm, A,chi,psi,chi_norm,RsdBiCG,Ncb,n_count)

LINEAR_OPERATOR(A);
multi1d<LatticeFermion> psi(Ncb);
multi1d<LatticeFermion> chi(Ncb);
int n_count;
Double chi_norm;
int Ncb;
Real RsdBiCG;

{ /* Local variables */
  include(COMMON_DECLARATIONS)

  PROTOTYPE(A, DESC, DESC, DESC, VAL, VAL)

  multi1d<LatticeFermion> r(Ncb);
  multi1d<LatticeFermion> rv(Ncb);
  multi1d<LatticeFermion> p(Ncb);
  multi1d<LatticeFermion> mp(Ncb);
  multi1d<LatticeFermion> s(Ncb);
  multi1d<LatticeFermion> ms(Ncb);
  multi1d<LatticeComplex> lc(Ncb);
  Complex a;
  Complex b;
  Complex omega;
  DComplex rvr;
  DComplex d;
  DComplex ct1;
  DComplex ct2;
  Complex ct3;
  Double re_rvr;
  Double im_rvr;
  Double one;
  Double zero;
  Double cp;
  Double c;
  Double ms2;
  Double rt1;
  Double rt2;
  Real rt3;
  Real rt4;
  Real one_s;
  Real re_a;
  Real im_a;
  Real re_b;
  Real im_b;
  int k;
  int cb;
  Real rsd_sq;

  START_CODE();

  rsd_sq = (RsdBiCG * chi_norm)* (RsdBiCG * chi_norm);

  one = 1;
  one_s = 1;
  zero = 0;

  
  /*  r[0]  :=  Chi - A . Psi[0] */
  /*  r  :=  A . Psi  */
  A (A, psi, r, Ncb, PLUS);

  /*  r  :=  ( Chi - R )  =  [ Chi  -  M(u)  . psi ] */
  r = -r;      	       	   /* Nc Ns  flops */
  r += chi;       	       	   /* Nc Ns  flops */

  /*  Cp = |r[0]|^2 */
  cp = 0;
  for(cb = 0; cb < Ncb; ++cb)
    cp += norm2(r[cb]);   	       	   /* 2 Nc Ns  flops */

  /*  IF |r[0]| <= RsdBiCG |Chi| THEN RETURN; */
/*  printf("ClovBiCG: k = 0  cp = %g\n",cp); */

  if ( TO_REAL(cp)  <=  rsd_sq )
  {
    n_count = 0;
    END_CODE();
    return;
  }

                          
  /*  p[1]  :=  r[0] */
  p = r;
  /*  rv    :=  r[0] */
  rv = r;

  /*  rvr = < rv, r[0] > = |r[0]|^2 = cp */
  rvr = cmplx(cp,zero);

  /*  FOR k FROM 1 TO MaxCG DO */
/*  printf("MaxCG = %d\n",MaxCG); */
/*  printf("(RsdBiCG * chi_norm)**2 = %g\n", rsd_sq; */
  for(k = 1; k <= MaxCG; ++k)
  {
    /*  c  =  | r[k-1] |**2 */
    c = cp;

    /*  mp = M(u) * p  */
    A (A, p, mp, Ncb, PLUS);

    /*  a[k] := < rv, r[k-1] > / < rv, mp > ; */
    /*  d = < rv, mp > */
        lc = trace(adj[rv] * mp);
    d = 0;
    for(cb = 0; cb < Ncb; ++cb)
      d += sum(lc[cb]);        	   /* ?? flops */
    
    /*  a = rvr / d */
    /*#ct1 = rvr / d; */
    /* Hack instead: */
    rt1 = localNorm2(d);
    ct2 = adj(d) * rvr;
    rt2 = one / rt1;
    ct1 = ct2 * rt2;
    a = FLOAT(ct1);

    /*  s = r - a mp */
    s = r;
    for(cb = 0; cb < Ncb; ++cb)
      s[cb] -= mp[cb] * a;  	   /* 2 Nc Ns  flops */

    /*  ms = A . s  */
    A (A, s, ms, Ncb, PLUS);

    /*  omega[k] := < ms[k], s[k] > / < ms[k], ms[k] > ; */
    /*  d = < ms, s > */
        lc = trace(adj[ms] * s);
    d = 0;
    for(cb = 0; cb < Ncb; ++cb)
      d += sum(lc[cb]);        	   /* ?? flops */
    
    /*  ms2 = | Ms | ** 2  */
    ms2 = 0;
    for(cb = 0; cb < Ncb; ++cb)
      ms2 += norm2(ms[cb]);       	   /* 2 Nc Ns  flops */

    /*  omega = d / ms2 */
    /*#ct1 = d / ms2; */
    /* Hack instead: */
    rt1 = one / ms2;
    ct1 = d * rt1;
    omega = FLOAT(ct1);

    /*  Psi[k] += omega[k] s[k] + a[k] p[k] ; */
    for(cb = 0; cb < Ncb; ++cb)
    {
      psi[cb] += s[cb] * omega;   /* 2 Nc Ns  flops */
      psi[cb] += p[cb] * a;  	   /* 2 Nc Ns  flops */
    }

    /*  r[k] := s[k] - omega[k] ms[k] ; */
    r = s;
    for(cb = 0; cb < Ncb; ++cb)
      r[cb] -= ms[cb] * omega;  	   /* 2 Nc Ns  flops */

    /*  cp  =  | r[k] |**2 */
    cp = 0;
    for(cb = 0; cb < Ncb; ++cb)
      cp += norm2(r[cb]);   	   /* 2 Nc Ns  flops */

/*    printf("ClovBiCG: k = %d  cp = %g\n",k, cp); */

    if ( TO_REAL(cp)  <=  rsd_sq )
    {
      n_count = k;
      END_CODE();
      return;
    }

    /*  b[k+1] := (<rv,r[k]> / <rv,r[k-1]>) * (a[k]/omega[k]) ; */
    /*  d = < rv, r[k-1] >  */
    d = rvr;

    /*  rvr = < rv, r[k] > */
        lc = trace(adj[rv] * r);
    rvr = 0;
    for(cb = 0; cb < Ncb; ++cb)
      rvr += sum(lc[cb]);        	   /* ?? flops */
    
    /*  b = rvr / d */
    /*#ct1 = rvr / d; */
    /* Hack instead: */
    rt1 = localNorm2(d);
    ct2 = adj(d) * rvr;
    rt2 = one / rt1;
    ct1 = ct2 * rt2;
    b = FLOAT(ct1);

    /*  b = b * a / omega */
    ct3 = b * a;
    /*#b = ct3 / omega; */
    /* Hack instead: */
    rt3 = localNorm2(omega);
    a = adj(omega) * ct3;
    rt4 = one_s / rt3;
    b = ct3 * rt4;

    /*  p[k+1] := r[k] + b[k+1] * (p[k] - omega[k] mp[k]) */
    for(cb = 0; cb < Ncb; ++cb)
    {
      p[cb] -= mp[cb] * omega;   /* 2 Nc Ns  flops */
      mp[cb] = p[cb] * b;        /* 2 Nc Ns  flops */
      p[cb] = mp[cb];
      p[cb] += r[cb];                   /* Nc Ns  flops */
    }
  }
  n_count = MaxCG;
  re_rvr = real(rvr);
  im_rvr = imag(rvr);
  re_a = real(a);
  im_a = imag(a);
  re_b = real(b);
  im_b = imag(b);

                            
  QDP_error_exit("too many BiCG iterations", n_count, rsd_sq, cp, c, re_rvr, im_rvr, re_a, im_a, re_b, im_b);
  END_CODE();
}

}  // end namespace Chroma
