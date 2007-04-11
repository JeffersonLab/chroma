// $Id: minvmr.cc,v 1.1 2007-04-11 03:40:15 edwards Exp $
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

    Complex a;
    Complex as;
    Complex atmp;
    Complex atmp2;
    Double d;
    Real asq;
    multi1d<Real> ass(Nmass);
    multi1d<Complex> rhos(Nmass);
    Boolean btmp;
    Boolean convP;
    multi1d<Boolean> convsP(Nmass);
    multi1d<Real> rsdmr_sq(Nmass);
/*  Real rsd_sq; */
 
    /* If exactly 0 norm, then solution must be 0 (for pos. def. operator) */
    if (chi_norm == 0.0)
    {
      n_count = 0;
      psi = zero;
      END_CODE();
      return;
    }

    T Ar;                    moveToFastMemoryHint(Ar);
    T chi_internal;          moveToFastMemoryHint(chi_internal);
    chi_internal[sub] = chi;
    moveToFastMemoryHint(psi,true);

    FlopCounter flopcount;
    flopcount.reset();
    StopWatch swatch;
    swatch.reset();
    swatch.start();

    Double cp = chi_norm * chi_norm;
    rsdmr_sq = localNorm2(RsdMR);
/*  rsd_sq = cp * rsdmr_sq; */

    // For this algorithm, all the psi have to be 0 to start
    // Only is that way the initial residuum r = chi
    for(int i= 0; i < n_shift; ++i) { 
      psi[i][sub] = zero;
    }
  
    /* Psi[0] := 0; */
    /* r[0] := Chi */
    T r = chi;
    r[sub] = chi_internal;
    FILL(convP,FALSE);
    FILL(convsP,FALSE);
    rhos = 1;

    /*  FOR k FROM 1 TO MaxMR DO */
    /*  IF |psi[k+1] - psi[k]| <= RsdMR |psi[k+1]| THEN RETURN; */
    for(k = 1; k <= MaxMR && ! convP ; ++k)
    {
      /*  cp = |r[k-1]|^2 */
      cp = norm2(r, sub);   	       	         flopcount.addSiteFlops(4*Nc*Ns,sub);

      /*  a[k-1] := < A.r[k-1], r[k-1] >/ < A.r[k-1], A.r[k-1] > ; */
      /*  Ar = A * r  */
      A(Ar, r, PLUS);                            flopcount.addFlops(A.nFlops());
      Ar[sub] += r * mass[isz];                  flopcount.addSiteFlops(4*Nc*Ns,sub);

      /*  c = < A.r, r > */
      DComplex c = innerProduct(Ar, r, sub);    flopcount.addSiteFlops(4*Nc*Ns,sub);

      /*  d = | A.r | ** 2  */
      d = norm2(Ar, sub);   	       	        flopcount.addSiteFlops(4*Nc*Ns,sub);

      /*  a = c / d */
      a = c / d;
        
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
	Real cs = Real(ass[s]) * cp;
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
    for(int s = 0; s < Nmass; ++s)
    {
      A(Ar, psi[s], PLUS);
      Ar[sub] += psi[s] * mass[s];
      Ar[sub] -= chi;

      cp = norm2(Ar, sub);   	       	         flopcount.addSiteFlops(4*Nc*Ns,sub);
      QDPIO::cout << "MInvMR (conv): s = " << s << "  r = " << sqrt(cp) << "  m = " << mass[s] << endl;
    }
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
