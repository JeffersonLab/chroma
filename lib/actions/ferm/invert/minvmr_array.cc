// $Id: minvmr_array.cc,v 1.1 2007-04-11 03:40:15 edwards Exp $
/*! \file
 *  \brief Multishift Minimal-residual algorithm for a Linear Operator
 */

#error "THIS ROUTINE IS NOT FULLY CONVERTED"

#include "linearop.h"
#include "actions/ferm/invert/minvmr.h"

namespace Chroma 
{

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
   *  z[0]   := 1 / (1 - (shift - shift(0))*b) 
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
   *  shifts        shifts of form  A + mass       (Read)
   *  RsdCG       residual accuracy              (Read/Write)
   *  n_count     Number of CG iteration	       (Write)

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
   *
   * @{
   */

  template<typename T>
  void MInvMR_a(const LinearOperator<T>& A, 
		const T& chi, 
		multi1d<T>& psi,
		const multi1d<Real>& shifts, 
		const multi1d<Real>& RsdMR, 
		int MaxMR,
		int& n_count)
  {
    START_CODE();

    const Subset& sub = A.subset();

    if (shifts.size() != RsdMR.size()) 
    {
      QDPIO::cerr << "MInvMR: number of shifts and residuals must match" << endl;
      QDP_abort(1);
    }

    int n_shift = shifts.size();

    if (n_shift == 0) 
    {
      QDPIO::cerr << "MInvMR: You must supply at least 1 mass: mass.size() = " 
		  << n_shift << endl;
      QDP_abort(1);
    }

    multi1d<LatticeFermion> ltmp(Ncb);
    multi1d<LatticeFermion> r(Ncb);
    multi1d<LatticeFermion> Ar(Ncb);
    multi1d<LatticeComplex> lc(Ncb);
    Complex a;
    Complex as;
    Complex atmp;
    Complex atmp2;
    DComplex b;
    DComplex c;
    Double d;
    Double one;
    Double dd;
    Double cp;
    Double cs;
    Real ones;
    Real asq;
    multi1d<Real> ass(Nmass);
    multi1d<Complex> rhos(Nmass);
    Boolean btmp;
    Boolean convP;
    multi1d<Boolean> convsP(Nmass);
    int iz;
    int k;
    int cb;
    int s;
    multi1d<Real> rsdmr_sq(Nmass);
/*  Real rsd_sq; */

    one = 1;
    ones = 1;
  
    /* If exactly 0 norm, then solution must be 0 (for pos. def. operator) */
    if (chi_norm == 0.0)
    {
      n_count = 0;
      psi = zero;
      END_CODE();
      return;
    }

    cp = chi_norm * chi_norm;
    rsdmr_sq = localNorm2(RsdMR);
/*  rsd_sq = cp * rsdmr_sq; */

              
    /* Psi[0] := 0; */
    /* r[0] := Chi */
    psi = zero;
    r = chi;
    FILL(convP,FALSE);
    FILL(convsP,FALSE);
    rhos = 1;

    /*  FOR k FROM 1 TO MaxMR DO */
    /*  IF |psi[k+1] - psi[k]| <= RsdMR |psi[k+1]| THEN RETURN; */
    for(k = 1; k <= MaxMR && ! convP ; ++k)
    {
      /*  cp = |r[k-1]|^2 */
      cp = zero;
      for(cb = 0; cb < Ncb; ++cb)
	cp += norm2(r[cb]);      	           /* 2 Nc Ns  flops */

      /*  a[k-1] := < A.r[k-1], r[k-1] >/ < A.r[k-1], A.r[k-1] > ; */
      /*  Ar = A * r  */
      A (A, r, Ar, Ncb, PLUS);

      for(cb = 0; cb < Ncb; ++cb)
	Ar[cb] += r[cb] * mass[isz];

      /*  c = < A.r, r > */
      lc = trace(adj[Ar] * r);
      c = zero;
      for(cb = 0; cb < Ncb; ++cb)
	c += sum(lc[cb]);		        /* ?? flops */

      /*  d = | A.r | ** 2  */
      d = zero;
      for(cb = 0; cb < Ncb; ++cb)
	d += norm2(Ar[cb]);	                /* 2 Nc Ns  flops */

      /*  a = c / d */
      /*#b = c / d; */
      dd = one / d;
      b = c * dd;
      a = FLOAT(b);
        
      /* Compute the shifted as and rhos */
      /*  Psi[k+1] += a[k] r[k] ; */
      for(s = 0; s < Nmass; ++s)
      {
	if (convsP[s]) continue;

	atmp = 1;
	asq = mass[s] - mass[isz];
	atmp += adj(a) * asq;
	asq = localNorm2(atmp);
	asq = ones / asq;
	atmp = atmp * asq;

	atmp2 = a * atmp;
	as = rhos[s] * atmp2;
	atmp2 = rhos[s] * atmp;
	rhos[s] = atmp2;

	ass[s] = localNorm2(as);

	for(cb = 0; cb < Ncb; ++cb)
	  psi[s][cb] += r[cb] * as;	/* 2 Nc Ns  flops */
      }
            
      /*  r[k] -= a[k-1] M . r[k-1] ; */
      for(cb = 0; cb < Ncb; ++cb)
	r[cb] -= Ar[cb] * a;	/* 2 Nc Ns  flops */

#if 0
      /* Project out eigenvectors */
      GramSchm (r, 1, OperEigVec, NOperEig, Ncb);
      GramSchm (psi, Nmass, OperEigVec, NOperEig, Ncb);
#endif
    
      /*  IF |psi[k+1] - psi[k]| <= RsdMR |psi[k+1]| THEN RETURN; */
      for(s = 0; s < Nmass; ++s)
      {
	if (convsP[s]) continue;

	d = zero;
	for(cb = 0; cb < Ncb; ++cb)
	{
	  d += norm2(psi[s][cb]);         	/* 2 Nc Ns  flops */
	}

	d *= rsdmr_sq[s];
	cs = FLOAT(ass[s]);
	cs = cs * cp;
	convsP[s] = cs < d;

#if 0
	PRINTF("MInvMR: k = %d  s = %d  r = %g  d = %g\n",k,s,sqrt(cs),sqrt(d));
	FLUSH_WRITE_NAMELIST(stdout);
#endif
      }
      /* Final convergence is determined by smallest shift */
      convP = convsP[isz];

      n_count = k;
    }

#if 0
    /* HACK **/
    for(s = 0; s < Nmass; ++s)
    {
      for(cb = 0; cb < Ncb; ++cb)
	ltmp[cb] = psi[s][cb];
      A (A, ltmp, Ar, Ncb, PLUS);
      for(cb = 0; cb < Ncb; ++cb)
	Ar[cb] += psi[s][cb] * mass[s];

      Ar -= chi;

      cp = zero;
      for(cb = 0; cb < Ncb; ++cb)
	cp += norm2(Ar[cb]);	                /* 2 Nc Ns  flops */

      PRINTF("MInvMR (conv): s = %d  r = %g  m = %g\n",s,sqrt(cp),mass[s]);
    }
    PRINTF("\n");
    /* end */
#endif

    if (n_count == MaxMR)
      QDP_error_exit("too many MR iterations", n_count, rsdmr_sq[isz]);
    END_CODE();
    return;
  }



  /*! \ingroup invert */
  template<>
  void MInvMR(const LinearOperator<LatticeFermion>& M,
	      const LatticeFermion& chi, 
	      multi1d<LatticeFermion>& psi, 
	      const multi1d<Real>& shifts,
	      const multi1d<Real>& RsdMR, 
	      int MaxMR,
	      int &n_count)
  {
    MInvMR_a(M, chi, psi, shifts, RsdMR, MaxMR, n_count);
  }


  /*! \ingroup invert */
  template<>
  void MInvMR(const DiffLinearOperator<LatticeFermion,
	                               multi1d<LatticeColorMatrix>,
	                               multi1d<LatticeColorMatrix> >& M,
	      const LatticeFermion& chi, 
	      multi1d<LatticeFermion>& psi, 
	      const multi1d<Real>& shifts,
	      const multi1d<Real>& RsdMR, 
	      int MaxMR,
	      int &n_count)
  {
    MInvMR_a(M, chi, psi, shifts, RsdMR, MaxMR, n_count);
  }

  /*! @} */  // end of group invert

}  // end namespace Chroma
