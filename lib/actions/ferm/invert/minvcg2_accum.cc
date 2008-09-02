// $Id: minvcg2_accum.cc,v 3.1 2008-09-02 20:10:18 bjoo Exp $

/*! \file
 *  \brief Multishift Conjugate-Gradient algorithm for a Linear Operator
 */

#include "linearop.h"
#include "actions/ferm/invert/minvcg2_accum.h"

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

  template<typename T>
  void MInvCG2Accum_a(const LinearOperator<T>& M, 
		 const T& chi, 
		 T& psi,
		 const Real& norm,
		 const multi1d<Real>& residues,
		 const multi1d<Real>& poles, 
		 const Real& RsdCG, 
		 int MaxCG,
		 int& n_count)
  {
    START_CODE();
    const Subset& sub = M.subset();


    int n_shift = poles.size();
    int s;
    Double cp;
    if (n_shift == 0) 
    {
      QDPIO::cerr << "MInvCGAccum: You must supply at least 1 mass: mass.size() = " 
		  << n_shift << endl;
      QDP_abort(1);
    }

    if (n_shift != residues.size()) { 
      QDPIO::cerr << "MInvCGAccum: Different number of poles and residues: n_oles=" << n_shift << " n_residues="<< residues.size() << endl;
      QDP_abort(1);
    }

    // Put in the A scale 
    multi1d<Real> scaled_res(n_shift);
    for(int i=0; i < n_shift; i++) { 
      scaled_res[i] = norm*residues(i);
    }


    /* Now find the smallest mass */
    int isz = 0;
    for(int findit=1; findit < n_shift; ++findit) {
      if ( toBool( poles[findit] < poles[isz])  ) { 
	isz = findit;
      }
    }


    // For this algorithm, all the psi have to be 0 to start
    // Only is that way the initial residuum r = chi
    psi[sub] = zero;
  
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

    // Check RHS is not too small 
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

    Double rsdcg_sq = RsdCG * RsdCG;
    Double rsd_sq = chi_norm_sq*rsdcg_sq;

  
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

    multi1d<T> X(n_shift); moveToFastMemoryHint(X);
    for(s = 0; s < n_shift; ++s) { 
      X[s][sub] = zero;
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
    r[sub] += Real(b)*MMp;                        flopcount.addSiteFlops(4*Nc*Ns,sub);

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
      // z[iz][s] = Double(1) / (Double(1) - (Double(poles[s])-Double(poles[isz]))*b);
      z[iz][s] = Double(1) / (Double(1) - Double(poles[s])*b);
      bs[s] = b * z[iz][s];
    }


    //  Psi[1] -= b[0] p[0] = - b[0] chi;
    for(s = 0; s < n_shift; ++s) {
      X[s][sub] = -Real(bs[s])*chi_internal;  flopcount.addSiteFlops(2*Nc*Ns,sub);
    }
  
    //  c = |r[1]|^2   
    Double c = norm2(r,sub);   	       	         flopcount.addSiteFlops(4*Nc*Ns,sub);

    // Check convergence of first solution
    multi1d<bool> convsP(n_shift);
    for(s = 0; s < n_shift; ++s) {
      convsP[s] = false;
    }

    bool convP = false;

    psi[sub] = scaled_res[0]*X[0];     flopcount.addSiteFlops(2*Nc*Ns,sub);
    for(s=1; s < n_shift; ++s) {
      psi[sub] += scaled_res[s]*X[s];  flopcount.addSiteFlops(4*Nc*Ns,sub);
    }


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
      p_0[sub] = r + Real(a)*p_0;                        flopcount.addSiteFlops(4*Nc*Ns,sub);
      //  p[k+1] := r[k+1] + a[k+1] p[k]; 
      //  Compute the shifted as */
      //  ps[k+1] := zs[k+1] r[k+1] + a[k+1] ps[k];
      for(s = 0; s < n_shift; ++s) {

	  // Don't update other p-s if converged.
	if( ! convsP[s] ) {

	  as = a * z[iz][s]*bs[s] / (z[1-iz][s]*b);
	  p[s][sub] = Real(z[iz][s])*r + Real(as)*p[s];  flopcount.addSiteFlops(6*Nc*Ns,sub);
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
      r[sub] += Real(b)*MMp;                                flopcount.addSiteFlops(4*Nc*Ns,sub);
      //  c  =  | r[k] |**2 
      c = norm2(r,sub);	                                   flopcount.addSiteFlops(4*Nc*Ns,sub);

      // Compute the shifted bs and z 
      iz = 1 - iz;

      for(s = 0; s < n_shift; s++) {
      	if ( !convsP[s] ) {
	  z0 = z[1-iz][s];
	  z1 = z[iz][s];
	  z[iz][s] = z0*z1*bp;
	  z[iz][s] /= b*a*(z1-z0) + z1*bp*(Double(1) - poles[s]*b);
	  bs[s] = b*z[iz][s]/z0;
	}
      }



      //  X[k+1] -= b[k] p[k] ; 
      for(s = 0; s < n_shift; ++s) {
	if (! convsP[s] ) {
	  X[s][sub] = Real(bs[s])*p[s];                 flopcount.addSiteFlops(2*Nc*Ns,sub);
	}
      }


      //    IF |psi[k+1] - psi[k]| <= RsdCG |psi[k]| THEN RETURN;
      // or IF |r[k+1]| <= RsdCG |chi| THEN RETURN;
      convP = true;
      
      for(s = 0; s < n_shift; s++) {
	if (! convsP[s] ) {
	  // Convergence methods 
	  // Check norm of shifted residuals 
	  Double css = c * z[iz][s]* z[iz][s];
	  convsP[s] = toBool( css < rsd_sq ); 
	}
	convP &= convsP[s];
      }

      T vectmp;
      vectmp[sub] = scaled_res[0]*X[0];    flopcount.addSiteFlops(4*Nc*Ns,sub);
      for(s=1; s < n_shift; s++) { 
	vectmp[sub] += scaled_res[s]*X[s]; flopcount.addSiteFlops(4*Nc*Ns,sub);
      }
      
      Double delta_s_norm = norm2(vectmp,sub); flopcount.addSiteFlops(4*Nc*Ns,sub);

      // Sum up the psi
      psi[sub] -= vectmp; flopcount.addSiteFlops(2*Nc*Ns,sub);
      Double s_norm = norm2(psi,sub);  flopcount.addSiteFlops(4*Nc*Ns,sub);

      bool convPSum = true;
      convP |= toBool( delta_s_norm <  rsdcg_sq*s_norm );
      
      n_count = k;
    } 

    swatch.stop();

    QDPIO::cout << "MInvCG2Accum: " << n_count << " iterations" << endl;
    flopcount.report("MInvCG2Accum", swatch.getTimeInSeconds());
    revertFromFastMemoryHint(psi,true);

    if (n_count == MaxCG) {
      QDP_error_exit("too many CG iterationns: %d\n", n_count);
    }

    END_CODE();
    return;
  }


  /*! \ingroup invert */
  template<>
  void MInvCG2Accum(const LinearOperator<LatticeFermion>& M,
	       const LatticeFermion& chi, 
	       LatticeFermion& psi,
	       const Real& norm,
	       const multi1d<Real>& residues,
	       const multi1d<Real>& poles, 
	       const Real& RsdCG, 
	       int MaxCG,
	       int& n_count)
  {
    MInvCG2Accum_a(M, chi, psi, norm, residues, poles, RsdCG, MaxCG, n_count);
  }


  /*! \ingroup invert */
  template<>
  void MInvCG2Accum(const DiffLinearOperator<LatticeFermion,
	                               multi1d<LatticeColorMatrix>,
	                               multi1d<LatticeColorMatrix> >& M,
	       const LatticeFermion& chi, 
	       LatticeFermion& psi,
	       const Real& norm,
	       const multi1d<Real>& residues,
	       const multi1d<Real>& poles, 
	       const Real& RsdCG, 
	       int MaxCG,
	       int& n_count)
  {
    MInvCG2Accum_a(M, chi, psi, norm, residues, poles, RsdCG, MaxCG, n_count);
  }

}  // end namespace Chroma
