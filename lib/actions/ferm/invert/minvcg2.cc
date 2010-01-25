// $Id: minvcg2.cc,v 3.3 2009-02-13 20:17:46 bjoo Exp $

/*! \file
 *  \brief Multishift Conjugate-Gradient algorithm for a Linear Operator
 */

#include "linearop.h"
#include "actions/ferm/invert/minvcg2.h"
#undef PAT
#ifdef PAT
#include <pat_api.h>
#endif
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

   *   	    Chi  =  M^\dag M . Psi

   * Algorithm:

   *  Psi[0] :=  0;      	                      Zeroed
   *  r[0]   :=  Chi;                           Initial residual
   *  p[1]   :=  Chi ;	       	       	      Initial direction
   *  b[0]   := |r[0]|**2 / <M p[0], M p[0]> ;
   *  z[0]   := 1 / (1 - (shift - shift(0))*b) 
   *  bs[0]  := b[0] * z[0]  
   *  r[1] += b[k] M^\dag M . p[0] ; 	       	      New residual
   *  Psi[1] = - b[k] p[k] ;   	       	      Starting solution vector
   *  IF |r[0]| <= RsdCG |Chi| THEN RETURN;        Converged?
   *  FOR k FROM 1 TO MaxCG DO    	       	       CG iterations
   *      a[k] := |r[k]|**2 / |r[k-1]|**2 ;
   *      p[k] := r[k] + a[k] p[k-1];   	       New direction
   *      b[k+1] := |r[k]|**2 / <M p[k], Mp[k]> ;
   *      r[k+1] += b[k+1] M^\dag M . p[k] ; 	       	       New residual
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

  template<typename T, typename R>
  void MInvCG2_a(const LinearOperator<T>& M, 
		 const T& chi, 
		 multi1d<T>& psi,
		 const multi1d<R>& shifts, 
		 const multi1d<R>& RsdCG, 
		 int MaxCG,
		 int& n_count)
  {
    START_CODE();

    const Subset& sub = M.subset();

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

      QDPIO::cout << "MInvCG2: " << n_count << " iterations" << endl;
      flopcount.report("minvcg2", swatch.getTimeInSeconds());
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


    T p_0;                   moveToFastMemoryHint(p_0);
    p_0[sub] = chi_internal;


    // Psi[0] := 0;
    multi1d<T> p(n_shift);   moveToFastMemoryHint(p);
    for(s = 0; s < n_shift; ++s) {
      p[s][sub] = chi_internal;                     // no flops
    }


    //  b[0] := - | r[0] |**2 / < p[0], Ap[0] > ;/
    //  First compute  d  =  < p, A.p > 
    //  Ap = A . p  */
    T Mp, MMp;                    moveToFastMemoryHint(Mp); moveToFastMemoryHint(MMp);
    M(Mp, p_0, PLUS);                            flopcount.addFlops(M.nFlops());

    /*  d =  < M p, M.p >  */
    Double d = norm2(Mp, sub);   flopcount.addSiteFlops(4*Nc*Ns,sub);

    M(MMp, Mp, MINUS); flopcount.addFlops(M.nFlops());

    Double b = -cp/d;

    //  r[1] += b[0] A . p[0]; 
    R b_r = b;
    r[sub] += b_r*MMp;                        flopcount.addSiteFlops(4*Nc*Ns,sub);

    /* Compute the shifted bs and z */
    multi1d<Double> bs(n_shift);
    multi2d<Double> z(2, n_shift);
    int iz;

    /* -- These are no longer special... */
    /*
    z[0][isz] = Double(1);
    z[1][isz] = Double(1);
    bs[isz] = b;
    */

    iz = 1;

    for(s = 0; s < n_shift; ++s)
    {
      z[1-iz][s] = Double(1);
      // z[iz][s] = Double(1) / (Double(1) - (Double(shifts[s])-Double(shifts[isz]))*b);
      z[iz][s] = Double(1) / (Double(1) - Double(shifts[s])*b);
      bs[s] = b * z[iz][s];
    }


    //  Psi[1] -= b[0] p[0] = - b[0] chi;
    for(s = 0; s < n_shift; ++s) {
      R bs_r = bs[s];
      psi[s][sub] = - bs_r*chi_internal;  flopcount.addSiteFlops(2*Nc*Ns,sub);
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

      // Update p
      R a_r = a;
      p_0[sub] = r + a_r*p_0;                        flopcount.addSiteFlops(4*Nc*Ns,sub);
      //  p[k+1] := r[k+1] + a[k+1] p[k]; 
      //  Compute the shifted as */
      //  ps[k+1] := zs[k+1] r[k+1] + a[k+1] ps[k];
      for(s = 0; s < n_shift; ++s) {

	  // Don't update other p-s if converged.
	if( ! convsP[s] ) {

	  as = a * z[iz][s]*bs[s] / (z[1-iz][s]*b);
	  R zizs= z[iz][s];
	  R as_r = as;
	  p[s][sub] = zizs*r + as_r*p[s];  flopcount.addSiteFlops(6*Nc*Ns,sub);
	}
      }

      //  cp  =  | r[k] |**2 
      cp = c;

      //  b[k] := | r[k] |**2 / < p[k], Ap[k] > ;
      //  First compute  d  =  < p, A.p >  
      //  Ap = A . p 
      M(Mp, p_0, PLUS);                                 flopcount.addFlops(M.nFlops());

      /*  d =  < p, A.p >  */
      d = norm2(Mp, sub);                           flopcount.addSiteFlops(4*Nc*Ns,sub);

      M(MMp, Mp, MINUS);                                 flopcount.addFlops(M.nFlops());

      bp = b;
      b = -cp/d;
      //  r[k+1] += b[k] A . p[k] ; 
      b_r = b;
      r[sub] += b_r*MMp;                                flopcount.addSiteFlops(4*Nc*Ns,sub);
      //  c  =  | r[k] |**2 
      c = norm2(r,sub);	                                   flopcount.addSiteFlops(4*Nc*Ns,sub);

      // Compute the shifted bs and z 
      iz = 1 - iz;
      for(s = 0; s < n_shift; s++) {
      	if ( !convsP[s] ) {
	  z0 = z[1-iz][s];
	  z1 = z[iz][s];
	  z[iz][s] = z0*z1*bp;
	  z[iz][s] /= b*a*(z1-z0) + z1*bp*(Double(1) - shifts[s]*b);
	  bs[s] = b*z[iz][s]/z0;
	}
      }



      //  Psi[k+1] -= b[k] p[k] ; 
      for(s = 0; s < n_shift; ++s) 
      {
	if (! convsP[s] ) 
	{
	  R bs_r = bs[s];

	  psi[s][sub] -= bs_r*p[s];                 flopcount.addSiteFlops(2*Nc*Ns,sub);
	}
      }


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



	  convsP[s] = toBool( css < rsd_sq[s] );



	}
	convP &= convsP[s];
      }

      n_count = k;
    }

    swatch.stop();

#if 0
    // Check answers 
    for(int s=0; s < n_shift; s++) { 
      r[sub] = chi;
      M(Mp, psi[s], PLUS);
      M(MMp, Mp, MINUS);
      MMp[sub] += shifts[s]*psi[s];
      r[sub] -= MMp;
      Double rnorm2 = norm2(r, sub);
      Double chinorm2= norm2(chi,sub);
      QDPIO::cout << "shift("<<s<<"): || r || / || chi || = "<<sqrt(rnorm2/chinorm2)<< endl;
    }
#endif
    QDPIO::cout << "MInvCG2: " << n_count << " iterations" << endl;
    flopcount.report("minvcg", swatch.getTimeInSeconds());
    revertFromFastMemoryHint(psi,true);

    if (n_count == MaxCG) {
      QDP_error_exit("too many CG iterationns: %d\n", n_count);
    }

    END_CODE();
    return;
  }



  /*! \ingroup invert */

  void MInvCG2(const LinearOperator<LatticeFermionF>& M,
	      const LatticeFermionF& chi, 
	      multi1d<LatticeFermionF>& psi, 
	      const multi1d<RealF>& shifts,
	      const multi1d<RealF>& RsdCG, 
	      int MaxCG,
	      int &n_count)
  {
#ifdef PAT
    int ierr=PAT_region_begin(22, "MInvCG2LinOp");
#endif
    MInvCG2_a(M, chi, psi, shifts, RsdCG, MaxCG, n_count);
#ifdef PAT
    ierr=PAT_region_end(22);
#endif
  }

  /*! \ingroup invert */
  void MInvCG2(const LinearOperator<LatticeFermionD>& M,
	      const LatticeFermionD& chi, 
	      multi1d<LatticeFermionD>& psi, 
	      const multi1d<RealD>& shifts,
	      const multi1d<RealD>& RsdCG, 
	      int MaxCG,
	      int &n_count)
  {
#ifdef PAT
    int ierr=PAT_region_begin(22, "MInvCG2LinOp");
#endif
    MInvCG2_a(M, chi, psi, shifts, RsdCG, MaxCG, n_count);
#ifdef PAT
    ierr=PAT_region_end(22);
#endif
  }


  /*! \ingroup invert */

  void MInvCG2(const DiffLinearOperator<LatticeFermionF,
	                               multi1d<LatticeColorMatrixF>,
	                               multi1d<LatticeColorMatrixF> >& M,
	      const LatticeFermionF& chi, 
	      multi1d<LatticeFermionF>& psi, 
	      const multi1d<RealF>& shifts,
	      const multi1d<RealF>& RsdCG, 
	      int MaxCG,
	      int &n_count)
  {
#ifdef PAT
     int ierr=PAT_region_begin(23,"MInvCG2DiffLinOp");
#endif
    MInvCG2_a(M, chi, psi, shifts, RsdCG, MaxCG, n_count);
#ifdef PAT
     ierr=PAT_region_end(23);
#endif
  }


  void MInvCG2(const DiffLinearOperator<LatticeFermionD,
	                               multi1d<LatticeColorMatrixD>,
	                               multi1d<LatticeColorMatrixD> >& M,
	      const LatticeFermionD& chi, 
	      multi1d<LatticeFermionD>& psi, 
	      const multi1d<RealD>& shifts,
	      const multi1d<RealD>& RsdCG, 
	      int MaxCG,
	      int &n_count)
  {
#ifdef PAT
     int ierr=PAT_region_begin(23,"MInvCG2DiffLinOp");
#endif
    MInvCG2_a(M, chi, psi, shifts, RsdCG, MaxCG, n_count);
#ifdef PAT
     ierr=PAT_region_end(23);
#endif
  }

}  // end namespace Chroma
