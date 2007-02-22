// $Id: minvcg_array.cc,v 3.2 2007-02-22 21:11:46 bjoo Exp $

/*! \file
 *  \brief Multishift Conjugate-Gradient algorithm for a Linear Operator
 */

#include "linearop.h"
#include "actions/ferm/invert/minvcg_array.h"

using namespace QDP::Hints;

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

  template<typename T, typename C>
  void MInvCG_a(const C& A, 
		const multi1d<T>& chi, 
		multi1d< multi1d<T> >& psi,
		const multi1d<Real>& shifts, 
		const multi1d<Real>& RsdCG, 
		int MaxCG,
		int& n_count)
  {
    START_CODE();

    // Setup
    const Subset& sub = A.subset();
    int n_shift = shifts.size();
    int N = A.size();
    FlopCounter flopcount;
    StopWatch swatch;

    if (shifts.size() != RsdCG.size()) 
    {
      QDPIO::cerr << "MinvCG: number of shifts and residuals must match" << endl;
      QDP_abort(1);
    }

    if (n_shift == 0) {
      QDPIO::cerr << "MinvCG: You must supply at least 1 shift. shift.size() = " << n_shift << endl;
      QDP_abort(1);
    }

    /* Now find the smallest mass */
    int isz = 0;
    for(int findit=1; findit < n_shift; ++findit) {
      if ( toBool( shifts[findit] < shifts[isz])  ) { 
	isz = findit;
      }
    }

    // Now arrange the memory
    // The outer size of the multi1d is n_shift
    psi.resize(n_shift);

    // For this algorithm, all the psi have to be 0 to start
    // Only is that way the initial residuum r = chi
    for(int i= 0; i < n_shift; ++i) { 
      psi[i].resize(N);
    }

    multi1d<T> r(N);
    multi1d< multi1d<T> > p(n_shift);

    for(int s = 0; s < n_shift; ++s) {
      p[s].resize(N); 
    }

    multi1d<T> Ap(N);
    multi1d<T> chi_internal(N);

    // Now arrange the memory:
    // These guys get most hits so locate them first
    moveToFastMemoryHint(Ap);
    moveToFastMemoryHint(p[isz]); 
    moveToFastMemoryHint(r);
    moveToFastMemoryHint(psi[isz]);
    moveToFastMemoryHint(chi_internal);

    { 
      multi1d<T> blocker1(N); moveToFastMemoryHint(blocker1);
      multi1d<T> blocker2(N); moveToFastMemoryHint(blocker2);
  
      // Now the rest of the p-s
      for(int i=0; i < n_shift; i++) { 
	if( i!=isz ) { 
	  moveToFastMemoryHint(p[i]);
	  moveToFastMemoryHint(psi[i]);
	}
      }
    
      // Blockers get freed here
    }

    // initialise psi-s
    for(int i=0; i < n_shift; i++) { 
      for(int n=0; n < N; ++n) {
	psi[i][n][sub] = zero;
      }
    }

    // initialise chi internal
    for(int i=0; i < N; i++) { 
      chi_internal[i][sub] = chi[i];
    }

    // -------- All memory setup and copies done. Timer starts here
    QDPIO::cout << "MinvCG starting" << endl;
    flopcount.reset();
    swatch.reset();
    swatch.start();

    // If chi has zero norm then the result is zero
    Double chi_norm_sq = norm2(chi_internal,sub);
    flopcount.addSiteFlops(4*Nc*Ns*N,sub);

    Double chi_norm = sqrt(chi_norm_sq);

    if( toBool( chi_norm < fuzz )) { 
      n_count = 0;
      swatch.stop();
      QDPIO::cout << "MinvCG: Finished. Iters taken = " << n_count << endl;
      flopcount.report("MinvCGArray", swatch.getTimeInSeconds());

      // Revert psi-s
      for(int i=0; i < n_shift; i++) { 
	revertFromFastMemoryHint(psi[i], true);
      }

      END_CODE();
      return;
    }

    multi1d<Double> rsd_sq(n_shift);
    multi1d<Double> rsdcg_sq(n_shift);

    Double cp = chi_norm_sq;
    for(int s = 0; s < n_shift; ++s)  {
      rsdcg_sq[s] = RsdCG[s] * RsdCG[s];  // RsdCG^2
      rsd_sq[s] = Real(cp) * rsdcg_sq[s]; // || chi ||^2 RsdCG^2
    }

  
    // r[0] := p[0] := Chi 
    for(int n=0; n < N; ++n) {
      r[n][sub] = chi_internal[n];

      for(int s=0; s < n_shift; s++) {
	p[s][n][sub] = chi_internal[n];
      }
    }



    //  b[0] := - | r[0] |**2 / < p[0], Ap[0] > ;/
    //  First compute  d  =  < p, A.p > 
    //  Ap = A . p  */
    A(Ap, p[isz], PLUS);
    flopcount.addFlops(A.nFlops());

    for(int n=0; n < N; ++n) {
      Ap[n][sub] += shifts[isz] * p[isz][n];
    }
    flopcount.addSiteFlops(4*Nc*Ns*N, sub);

    /*  d =  < p, A.p >  */
    Double d = innerProductReal(p[isz], Ap, sub);
    flopcount.addSiteFlops(4*Nc*Ns*N ,sub);

    Double b = -cp/d;

    /* Compute the shifted bs and z */
    multi1d<Double> bs(n_shift);
    multi2d<Double> z(2, n_shift);
    int iz;

    z[0][isz] = Double(1);
    z[1][isz] = Double(1);
    bs[isz] = b;
    iz = 1;

    for(int s = 0; s < n_shift; ++s)
    {
      if( s != isz ) {
	z[1-iz][s] = Double(1);
	z[iz][s] = Double(1) / (Double(1) - (Double(shifts[s])-Double(shifts[isz]))*b);
	bs[s] = b * z[iz][s];
      }
    }

    //  r[1] += b[0] A . p[0]; 
    for(int n=0; n < N; ++n) {
      r[n][sub] += Real(b)* Ap[n];
    }
    flopcount.addSiteFlops(4*Nc*Ns*N, sub);

    //  Psi[1] -= b[0] p[0] = - b[0] chi;
    // Psi[0] are all 0 so I can write this as a -= bs*chi
    for(int s = 0; s < n_shift; ++s) {
      for(int n=0; n < N; ++n) {
	psi[s][n][sub] -= Real(bs[s])*chi_internal[n];
      }
    }
    flopcount.addSiteFlops(4*Nc*Ns*N*n_shift, sub);

    //  c = |r[1]|^2   
    Double c = norm2(r,sub);   	       	        
    flopcount.addSiteFlops(4*Nc*Ns*N,sub);


    // Check convergence of first solution
    multi1d<bool> convsP(n_shift);
    for(int s = 0; s < n_shift; ++s) {
      convsP[s] = false;
    }

    bool convP = toBool( c < rsd_sq[isz] );


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
      for(int s = 0; s < n_shift; ++s) {

	// Always update p[isz] even if isz is converged
	// since the other p-s depend on it.
	if (s == isz) {

	  for(int n=0; n < N; ++n) {
	    p[s][n][sub] = r[n] + Real(a)*p[s][n];
	    //      p[s][n][sub] *= Real(a);	        
	    //	  p[s][n][sub] += r[n];	                
	  }
	  flopcount.addSiteFlops(4*Nc*Ns*N,sub);

	}
	else {
	  // Don't update other p-s if converged.
	  if( ! convsP[s] ) { 
	    as = a * z[iz][s]*bs[s] / (z[1-iz][s]*b);
	  
	    for(int n=0; n < N; ++n) {
	      p[s][n][sub] = Real(as)*p[s][n] + Real(z[iz][s])*r[n];
	      //p[s][n][sub] *= Real(as);	         
	      //p[s][n][sub] += Real(z[iz][s])*r[n];	 
	    }
	    flopcount.addSiteFlops(6*Nc*Ns*N, sub);
	  }
	}

      }

      //  cp  =  | r[k] |**2 
      cp = c;

      //  b[k] := | r[k] |**2 / < p[k], Ap[k] > ;
      //  First compute  d  =  < p, A.p >  
      //  Ap = A . p 
      A(Ap, p[isz], PLUS);
      flopcount.addFlops(A.nFlops());

      for(int n=0; n < N; ++n) {
	Ap[n][sub] += shifts[isz] *p[isz][n];
      }
      flopcount.addSiteFlops(4*Nc*Ns*N, sub);

      /*  d =  < p, A.p >  */
      d = innerProductReal(p[isz], Ap, sub);
      flopcount.addSiteFlops(4*Nc*Ns*N, sub);
    
      bp = b;
      b = -cp/d;

      // Compute the shifted bs and z 
      bs[isz] = b;
      iz = 1 - iz;
      for(int s = 0; s < n_shift; s++) {
      
	if (s != isz && !convsP[s] ) {
	  z0 = z[1-iz][s];
	  z1 = z[iz][s];
	  z[iz][s] = z0*z1*bp;
	  z[iz][s] /= b*a*(z1-z0) + z1*bp*(Double(1) - (shifts[s] - shifts[isz])*b);
	  bs[s] = b*z[iz][s]/z0;
	}
      }

      //  r[k+1] += b[k] A . p[k] ; 
      for(int n=0; n < N; ++n) {
	r[n][sub] += Real(b)*Ap[n];
      }
      flopcount.addSiteFlops(4*Nc*Ns*N, sub);

      //  Psi[k+1] -= b[k] p[k] ; 
      for(int s = 0; s < n_shift; ++s) {
	if (! convsP[s] ) {
	  for(int n=0; n < N; ++n) {
	    psi[s][n][sub] -= Real(bs[s])*p[s][n];
	  }
	  flopcount.addSiteFlops(4*Nc*Ns*N, sub);
	}
      }

      //  c  =  | r[k] |**2 
      c = norm2(r,sub);	                
      flopcount.addSiteFlops(4*Nc*Ns*N, sub);

      //    IF |psi[k+1] - psi[k]| <= RsdCG |psi[k+1]| THEN RETURN;
      // or IF |r[k+1]| <= RsdCG |chi| THEN RETURN;
      convP = true;
      for(int s = 0; s < n_shift; s++) {
	if (! convsP[s] ) {


	  // Convergence methods 
	  // Check norm of shifted residuals 
	  Double css = c * z[iz][s]* z[iz][s];
	  convsP[s] = toBool(  css < rsd_sq[s] );

	}
	convP &= convsP[s];
      }

      n_count = k;
    }

    swatch.stop();
    QDPIO::cout << "MinvCGArray finished: " << n_count << " iterations " << endl;
    flopcount.report("MinvCGArray", swatch.getTimeInSeconds());
    QDPIO::cout << "MinvCG Explicit solution check: " << endl;

    // Expicitly check the ALL solutions
    for(int s = 0; s < n_shift; ++s)
    {
      A(Ap, psi[s], PLUS);
      for(int n=0; n < N; ++n) {
	Ap[n][sub] += shifts[s] * psi[s][n];

	Ap[n][sub] -= chi_internal[n];
      }

      c = norm2(Ap,sub);	              

      QDPIO::cout << "MInvCG (conv): s = " << s 
		  << " shift = " << shifts[s]
		  << " r = " <<  Real(sqrt(c)/chi_norm) << endl;
		
    }

    for(int i=0; i < n_shift; i++) { 
      revertFromFastMemoryHint(psi[i],true);
    }

    if (n_count == MaxCG) {
      QDPIO::cout << "too many CG iterationns: " << n_count << endl;
      QDP_abort(1);

    }
    END_CODE();
    return;
  }


  /*! \ingroup invert */
  template<>
  void MInvCG(const LinearOperatorArray<LatticeFermion>& M,
	      const multi1d<LatticeFermion>& chi, 
	      multi1d< multi1d<LatticeFermion> >& psi, 
	      const multi1d<Real>& shifts,
	      const multi1d<Real>& RsdCG, 
	      int MaxCG,
	      int &n_count)
  {
    MInvCG_a(M, chi, psi, shifts, RsdCG, MaxCG, n_count);
  }


  /*! \ingroup invert */
  template<>
  void MInvCG(const DiffLinearOperatorArray<LatticeFermion, 
	                                    multi1d<LatticeColorMatrix>, 
	                                    multi1d<LatticeColorMatrix> >& M,
	      const multi1d<LatticeFermion>& chi, 
	      multi1d< multi1d<LatticeFermion> >& psi, 
	      const multi1d<Real>& shifts,
	      const multi1d<Real>& RsdCG, 
	      int MaxCG,
	      int &n_count)
  {
    MInvCG_a(M, chi, psi, shifts, RsdCG, MaxCG, n_count);
  }

}  // end namespace Chroma
