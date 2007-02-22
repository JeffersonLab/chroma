// $Id: minvcg.cc,v 3.6 2007-02-22 21:11:46 bjoo Exp $

/*! \file
 *  \brief Multishift Conjugate-Gradient algorithm for a Linear Operator
 */

#include "linearop.h"
#include "actions/ferm/invert/minvcg.h"

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
   */

  template<typename T>
  void MInvCG_a(const LinearOperator<T>& A, 
		const T& chi, 
		multi1d<T>& psi,
		const multi1d<Real>& shifts, 
		const multi1d<Real>& RsdCG, 
		int MaxCG,
		int& n_count)
  {
    START_CODE();

    const Subset& sub = A.subset();

    if (shifts.size() != RsdCG.size()) 
    {
      QDPIO::cerr << "MInvCG: number of shifts and residuals must match" << endl;
      QDP_abort(1);
    }

    int n_shift = shifts.size();

    if (n_shift == 0) 
    {
      QDPIO::cerr << "MInvCG: You must supply at least 1 mass: mass.size() = " 
		  << n_shift << endl;
      QDP_abort(1);
    }

    /* Now find the smallest mass */
    int isz = 0;
    for(int findit=1; findit < n_shift; ++findit) {
      if ( toBool( shifts[findit] < shifts[isz])  ) { 
	isz = findit;
      }
    }

#if 0 
    QDPIO::cout << "n_shift = " << n_shift << " isz = " << isz << " shift = " << shifts[0] << endl;
#endif

    // We need to make sure, that psi is at least as big as the number
    // of shifts. We resize it if it is not big enough.
    // However, it is allowed to be bigger.
    if( psi.size() <  n_shift ) { 
      psi.resize(n_shift);
    }

    // For this algorithm, all the psi have to be 0 to start
    // Only is that way the initial residuum r = chi
    for(int i= 0; i < n_shift; ++i) { 
      psi[i][sub] = zero;
    }
  
    T chi_internal;      moveToFastMemoryHint(chi_internal);
    chi_internal[sub] = chi;
    moveToFastMemoryHint(psi,true);

    FlopCounter flopcount;
    flopcount.reset();
    StopWatch swatch;
    swatch.reset();
    swatch.start();

    // If chi has zero norm then the result is zero
    Double chi_norm_sq = norm2(chi_internal,sub);    flopcount.addSiteFlops(4*Nc*Ns,sub);
    Double chi_norm = sqrt(chi_norm_sq);

    if( toBool( chi_norm < fuzz )) 
    {
      swatch.stop();

      n_count = 0;

      QDPIO::cout << "MInvCG: " << n_count << " iterations" << endl;
      flopcount.report("minvcg", swatch.getTimeInSeconds());
      revertFromFastMemoryHint(psi,true);

      // The psi are all zero anyway at this point
      // for(int i=0; i < n_shift; i++) { psi[i] = zero; }
      END_CODE();
      return;
    }

    multi1d<Double> rsd_sq(n_shift);
    multi1d<Double> rsdcg_sq(n_shift);

    Double cp = chi_norm_sq;
    int s;
    for(s = 0; s < n_shift; ++s)  {
      rsdcg_sq[s] = RsdCG[s] * RsdCG[s];  // RsdCG^2
      rsd_sq[s] = Real(cp) * rsdcg_sq[s]; // || chi ||^2 RsdCG^2
    }

  
    // r[0] := p[0] := Chi 
    T r;                     moveToFastMemoryHint(r);
    r[sub] = chi_internal;                          // no flops

    // Psi[0] := 0;
    multi1d<T> p(n_shift);   moveToFastMemoryHint(p);
    for(s = 0; s < n_shift; ++s) {
      p[s][sub] = chi_internal;                     // no flops
    }


    //  b[0] := - | r[0] |**2 / < p[0], Ap[0] > ;/
    //  First compute  d  =  < p, A.p > 
    //  Ap = A . p  */
    T Ap;                    moveToFastMemoryHint(Ap);
    A(Ap, p[isz], PLUS);                            flopcount.addFlops(A.nFlops());
    Ap[sub] += shifts[isz]*p[isz];                  flopcount.addSiteFlops(4*Nc*Ns,sub);

    /*  d =  < p, A.p >  */
    Double d = innerProductReal(p[isz], Ap, sub);   flopcount.addSiteFlops(4*Nc*Ns,sub);


    Double b = -cp/d;

    /* Compute the shifted bs and z */
    multi1d<Double> bs(n_shift);
    multi2d<Double> z(2, n_shift);
    int iz;

    z[0][isz] = Double(1);
    z[1][isz] = Double(1);
    bs[isz] = b;
    iz = 1;

    for(s = 0; s < n_shift; ++s)
    {
      if( s != isz ) {
	z[1-iz][s] = Double(1);
	z[iz][s] = Double(1) / (Double(1) - (Double(shifts[s])-Double(shifts[isz]))*b);
	bs[s] = b * z[iz][s];
      }
    }

    //  r[1] += b[0] A . p[0]; 
    r[sub] += Real(b)*Ap;                        flopcount.addSiteFlops(4*Nc*Ns,sub);

    //  Psi[1] -= b[0] p[0] = - b[0] chi;
    for(s = 0; s < n_shift; ++s) {
      psi[s][sub] = - Real(bs[s])*chi_internal;  flopcount.addSiteFlops(2*Nc*Ns,sub);
    }
  
    //  c = |r[1]|^2   
    Double c = norm2(r,sub);   	       	         flopcount.addSiteFlops(4*Nc*Ns,sub);

    // Check convergence of first solution
    multi1d<bool> convsP(n_shift);
    for(s = 0; s < n_shift; ++s) {
      convsP[s] = false;
    }

    bool convP = toBool( c < rsd_sq[isz] );

#if 0 
    QDPIO::cout << "MInvCG: k = 0  r = " << sqrt(c) << endl;
#endif

    //  FOR k FROM 1 TO MaxCG DO
    //  IF |psi[k+1] - psi[k]| <= RsdCG |psi[k+1]| THEN RETURN; 
    Double z0, z1;
    Double ztmp;
    Double cs;
    Double a;
    Double as;
    Double  bp;
    int k;
  
    for(k = 1; k <= MaxCG && !convP ; ++k)
    {
      //  a[k+1] := |r[k]|**2 / |r[k-1]|**2 ; 
      a = c/cp;

      //  p[k+1] := r[k+1] + a[k+1] p[k]; 
      //  Compute the shifted as */
      //  ps[k+1] := zs[k+1] r[k+1] + a[k+1] ps[k];
      for(s = 0; s < n_shift; ++s) 
      {
	// Always update p[isz] even if isz is converged
	// since the other p-s depend on it.
	if (s == isz) 
	{
	  /* p[s][sub] *= Real(a);
	     p[s][sub] += r;
	  */
	  p[s][sub] = r + Real(a)*p[s];                    flopcount.addSiteFlops(4*Nc*Ns,sub);
	}
	else 
	{
	  // Don't update other p-s if converged.
	  if( ! convsP[s] ) 
	  {
	    as = a * z[iz][s]*bs[s] / (z[1-iz][s]*b);
	    /*
	      p[s][sub] *= Real(as);
	      p[s][sub] += Real(z[iz][s])*r;
	    */
	    p[s][sub] = Real(z[iz][s])*r + Real(as)*p[s];  flopcount.addSiteFlops(6*Nc*Ns,sub);
	  }
	}
      }

      //  cp  =  | r[k] |**2 
      cp = c;

      //  b[k] := | r[k] |**2 / < p[k], Ap[k] > ;
      //  First compute  d  =  < p, A.p >  
      //  Ap = A . p 
      A(Ap, p[isz], PLUS);                                 flopcount.addFlops(A.nFlops());
      Ap[sub] += shifts[isz]*p[isz];                       flopcount.addSiteFlops(4*Nc*Ns,sub);

      /*  d =  < p, A.p >  */
      d = innerProductReal(p[isz], Ap, sub);               flopcount.addSiteFlops(4*Nc*Ns,sub);
    
      bp = b;
      b = -cp/d;

      // Compute the shifted bs and z 
      bs[isz] = b;
      iz = 1 - iz;
      for(s = 0; s < n_shift; s++) 
      {
      	if (s != isz && !convsP[s] ) 
	{
	  z0 = z[1-iz][s];
	  z1 = z[iz][s];
	  z[iz][s] = z0*z1*bp;
	  z[iz][s] /= b*a*(z1-z0) + z1*bp*(Double(1) - (shifts[s] - shifts[isz])*b);
	  bs[s] = b*z[iz][s]/z0;
	}
      }

      //  r[k+1] += b[k] A . p[k] ; 
      r[sub] += Real(b)*Ap;                                flopcount.addSiteFlops(4*Nc*Ns,sub);


      //  Psi[k+1] -= b[k] p[k] ; 
      for(s = 0; s < n_shift; ++s) 
      {
	if (! convsP[s] ) 
	{
	  psi[s][sub] -= Real(bs[s])*p[s];                 flopcount.addSiteFlops(2*Nc*Ns,sub);
	}
      }

      //  c  =  | r[k] |**2 
      c = norm2(r,sub);	                                   flopcount.addSiteFlops(4*Nc*Ns,sub);

      //    IF |psi[k+1] - psi[k]| <= RsdCG |psi[k+1]| THEN RETURN;
      // or IF |r[k+1]| <= RsdCG |chi| THEN RETURN;
      convP = true;
      for(s = 0; s < n_shift; s++) 
      {
	if (! convsP[s] ) 
	{
	  // Convergence methods 
	  // Check norm of shifted residuals 
	  Double css = c * z[iz][s]* z[iz][s];

#if 0	
	  QDPIO::cout << "MInvCG (shift=" << s << ") k = " << k <<"  r =  " 
		      << css << " rsd_sq["<<s<<"] = " << rsd_sq[s] << endl;
#endif 

	  convsP[s] = toBool( css < rsd_sq[s] );

#if 0
     
	  // 
	  // Check relative error of solution 

	  // cs holds | beta_s p_s |^2 = | psi_next |^2
	  cs = norm2(p[s],sub);	                           flopcount.addSiteFlops(4*Nc*Ns,sub);
	  cs *= bs[s]*bs[s];

	  // d holds | psi |^2 * epsilon^2
	  d = norm2(psi[s],sub);                           flopcount.addSiteFlops(4*Nc*Ns,sub);
	  d *= rsdcg_sq[s];


	  // Terminate if | psi |^2/|psi_next|^2 < epsilon^2
	  convsP[s] = toBool( cs < d );

#if 0
	  QDPIO::cout  << "MInvCG (shift=" << s << ") k = " << k << " cs = " 
		       << cs << " d = " << d << endl;
#endif
#endif

	}
	convP &= convsP[s];
      }

      n_count = k;
    }

    swatch.stop();

#if 0
    // This is paranoia - do not need

    // Expicitly check the ALL solutions
    for(s = 0; s < n_shift; ++s)
    {
      A(Ap, psi[s], PLUS);
      Ap[sub] += shifts[s]*psi[s];

      Ap[sub] -= chi_internal;

      c = norm2(Ap,sub);	                /* 2 Nc Ns  flops */

      QDPIO::cout << "MInvCG (conv): s = " << s 
		  << " shift = " << shifts[s]
		  << " r = " <<  Real(sqrt(c)/chi_norm) << endl;
		
    }
    /* end */
#endif

    QDPIO::cout << "MInvCG: " << n_count << " iterations" << endl;
    flopcount.report("minvcg", swatch.getTimeInSeconds());
    revertFromFastMemoryHint(psi,true);

    if (n_count == MaxCG) {
      QDP_error_exit("too many CG iterationns: %d\n", n_count);
    }

    END_CODE();
    return;
  }


  /*! \ingroup invert */
  template<>
  void MInvCG(const LinearOperator<LatticeFermion>& M,
	      const LatticeFermion& chi, 
	      multi1d<LatticeFermion>& psi, 
	      const multi1d<Real>& shifts,
	      const multi1d<Real>& RsdCG, 
	      int MaxCG,
	      int &n_count)
  {
    MInvCG_a(M, chi, psi, shifts, RsdCG, MaxCG, n_count);
  }


  /*! \ingroup invert */
  template<>
  void MInvCG(const DiffLinearOperator<LatticeFermion,
	                               multi1d<LatticeColorMatrix>,
	                               multi1d<LatticeColorMatrix> >& M,
	      const LatticeFermion& chi, 
	      multi1d<LatticeFermion>& psi, 
	      const multi1d<Real>& shifts,
	      const multi1d<Real>& RsdCG, 
	      int MaxCG,
	      int &n_count)
  {
    MInvCG_a(M, chi, psi, shifts, RsdCG, MaxCG, n_count);
  }

}  // end namespace Chroma
