// -*- C++ -*-
// $Id: stout_utils.cc,v 1.6 2009-02-04 21:16:03 bjoo Exp $
/*! \file
 *  \brief Stout utilities
 */

#include "chroma_config.h"
#include "chromabase.h"
#include "util/gauge/stout_utils.h"

//#if defined(BUILD_JIT_CLOVER_TERM)
//#include "util/gauge/stout_utils_ptx.h"
//#endif

namespace Chroma 
{ 

  //! Timings
  /*! \ingroup gauge */
  namespace StoutLinkTimings { 
    static double smearing_secs = 0;
    double getSmearingTime() { 
      return smearing_secs;
    }

    static double force_secs = 0;
    double getForceTime() { 
      return force_secs;
    }

    static double functions_secs = 0;
    double getFunctionsTime() { 
      return functions_secs;
    }
  }

  //! Utilities
  /*! \ingroup gauge */
  namespace Stouting 
  {

    /*! \ingroup gauge */
    void getQs(const multi1d<LatticeColorMatrix>& u, LatticeColorMatrix& Q, 
	       LatticeColorMatrix& QQ,
	       int mu,
	       const multi1d<bool>& smear_in_this_dirP,
	       const multi2d<Real>& rho)
    {
      START_CODE();
      
      LatticeColorMatrix C;
      getQsandCs(u, Q, QQ, C, mu, smear_in_this_dirP,rho);

      END_CODE();
    }


    /*! \ingroup gauge */
    void getQsandCs(const multi1d<LatticeColorMatrix>& u, LatticeColorMatrix& Q, 
		    LatticeColorMatrix& QQ,
		    LatticeColorMatrix& C, 
		    int mu,
		    const multi1d<bool>& smear_in_this_dirP,
		    const multi2d<Real>& rho)
    {
      START_CODE();
      
      C = zero;
      
      // If rho is nonzero in this direction then accumulate the staples
      for(int nu=0; nu < Nd; nu++) 
      { 
	// Accumulate mu-nu staple
	if( (mu != nu) && smear_in_this_dirP[nu] ) 
	{
	  LatticeColorMatrix U_nu_plus_mu = shift(u[nu], FORWARD, mu);
	  LatticeColorMatrix tmp_mat;
	  LatticeColorMatrix tmp_mat2; 
	  
	
	  // Forward staple
	  //             2
	  //       ^ ---------> 
	  //       |          |
	  //    1  |          |  3
	  //       |          |
	  //       |          V
	  //       x          x + mu
	  //
	  tmp_mat = shift(u[mu],FORWARD, nu);
	  tmp_mat2 = u[nu]*tmp_mat;
	  tmp_mat  = tmp_mat2*adj(U_nu_plus_mu);
	  
	  
	  // Backward staple
	  //             
	  //       |          ^ 
	  //       |          |
	  //    1  |          |  3
	  //       |     2    |
	  //       V--------->|          
	  //       x-nu        x - nu  + mu
	  //
	  //
	  //  we construct it on x-nu and shift it up to x, 
	  // (with a backward shift)
	  
	  // This is the staple on x-nu:
	  // tmp_1(x) = u_dag(x,nu)*u(x,mu)*u(x+mu,nu)
	  {
	    LatticeColorMatrix tmp_mat3;
	    tmp_mat3 = adj(u[nu])*u[mu];
	    tmp_mat2 = tmp_mat3*U_nu_plus_mu;
	    tmp_mat += shift(tmp_mat2, BACKWARD, nu);
	    tmp_mat *= rho(mu,nu);
	  }
	  C += tmp_mat;
	  
	}
      }
      
      // Now I can form the Q
      LatticeColorMatrix Omega;
      Omega = C*adj(u[mu]); // Q_mu is Omega mu here (eq 2 part 2)
      
      LatticeColorMatrix tmp2 = adj(Omega) - Omega;
      LatticeColorMatrix tmp3 = trace(tmp2);
      tmp3 *= Real(1)/Real(Nc);
      tmp2 -= tmp3;
      tmp2 *= Real(0.5);
      Q = timesI(tmp2);
      QQ = Q*Q;
      
      END_CODE();
    }
    
    /*! \ingroup gauge */
    // Do the force recursion from level i+1, to level i
    // The input fat_force F is modified.
    void deriv_recurse(multi1d<LatticeColorMatrix>& F,
		       const multi1d<bool>& smear_in_this_dirP,
		       const multi2d<Real>& rho,
		       const multi1d<LatticeColorMatrix>& u)
    {
      START_CODE();
      
      // Things I need
      // C_{\mu} = staple multiplied appropriately by the rho
      // Lambda matrices asper eq(73) 
      multi1d<LatticeColorMatrix> F_plus(Nd);
      
      // Save the fat force
      F_plus = F;
      
      
      multi1d<LatticeColorMatrix> Lambda(Nd);
      multi1d<LatticeColorMatrix> C(Nd);
      
      // The links at this level (unprimed in the paper).
      //const multi1d<LatticeColorMatrix>& u = smeared_links[level];
      
      for(int mu=0; mu < Nd; mu++) 
      {
	if( smear_in_this_dirP[mu] ) 
	{ 
	  LatticeColorMatrix Q,QQ;   // This is the C U^{dag}_mu suitably antisymmetrized
	  
	  // Get Q, Q^2, C, c0 and c1 -- this code is the same as used in stout_smear()
	  getQsandCs(u, Q, QQ, C[mu], mu, smear_in_this_dirP,rho);
	  
	  // Now work the f-s and b-s
	  multi1d<LatticeComplex> f;
	  multi1d<LatticeComplex> b_1;
	  multi1d<LatticeComplex> b_2;
	  
	  // Get the fs and bs  -- does internal resize to make them arrays of length 3
	  //	QDPIO::cout << __func__ << ": mu=" << mu << endl;
	  getFsAndBs(Q,QQ, f, b_1, b_2, true);
	  
	  
	  LatticeColorMatrix B_1 = b_1[0] + b_1[1]*Q + b_1[2]*QQ;
	  LatticeColorMatrix B_2 = b_2[0] + b_2[1]*Q + b_2[2]*QQ;
	  
	  
	  // Construct the Gamma ( eq 74 and 73 )
	  LatticeColorMatrix USigma = u[mu]*F_plus[mu];
	  LatticeColorMatrix Gamma = f[1]*USigma + f[2]*(USigma*Q + Q*USigma)
	    + trace(B_1*USigma)*Q
	    + trace(B_2*USigma)*QQ;
	  
	  // Take the traceless hermitian part to form Lambda_mu (eq 72)
	  Lambda[mu] = Gamma + adj(Gamma);    // Make it hermitian
	  LatticeColorMatrix tmp3 = (Double(1)/Double(Nc))*trace(Lambda[mu]); // Subtract off the trace
	  Lambda[mu] -= tmp3;
	  Lambda[mu] *= Double(0.5);         // overall factor of 1/2
	  
	  // The first 3 terms of eq 75
	  // Now the Fat force * the exp(iQ)
	  F[mu]  = F_plus[mu]*(f[0] + f[1]*Q + f[2]*QQ);
	  
#if 0
	  QDPIO::cout << __func__ << ": F[" << mu << "]= " << norm2(F[mu]) 
		      << "  F_plus[mu]=" << norm2(F_plus[mu])
		      << "  f[0]=" << norm2(f[0]) 
		      << "  f[1]=" << norm2(f[1]) 
		      << "  f[2]=" << norm2(f[2]) 
		      << "  B_1=" << norm2(B_1) 
		      << "  B_2=" << norm2(B_2) 
		      << "  Q=" << norm2(Q) 
		      << "  QQ=" << norm2(QQ) 
		      << endl;
#endif
	  
	} // End of if( smear_in_this_dirP[mu] )
	// else what is in F_mu is the right force
      }
      
      // At this point we should have 
      //
      //  F[mu] = F_plus[mu]*exp(iQ) 
      //
      //  We need the 8 staple terms left in dOmega/dU (last 6 terms in eq 75 + the iC{+}Lambda
      //  term in eq 75 which in reality just covers 2 staples.
      
      //  We have to make this a separate loop from the above, because we need to know the 
      //  Lambda[mu] and [nu] for all the avaliable mu-nu combinations
      for(int mu = 0; mu < Nd; mu++) 
      { 
	if( smear_in_this_dirP[mu] ) 
	{ 
	  LatticeColorMatrix staple_sum = zero;
	  // LatticeColorMatrix staple_sum_dag = adj(C[mu])*Lambda[mu];
	  for(int nu = 0; nu < Nd; nu++) { 
	    if((mu != nu) && smear_in_this_dirP[nu] ) { 
	      LatticeColorMatrix U_nu_plus_mu = shift(u[nu],FORWARD, mu);
	      LatticeColorMatrix U_mu_plus_nu = shift(u[mu],FORWARD, nu);
	      LatticeColorMatrix Lambda_nu_plus_mu = shift(Lambda[nu], FORWARD, mu);
	      LatticeColorMatrix tmp_mat;
	      LatticeColorMatrix tmp_mat2;
	      
	      
	      //  THe three upward staples
	      //  Staples 1 5 and 6 in the paper
	      // 
	      //  Staple 1
	      //      rho(nu,mu)  *                ( [ U_nu(x+mu) U^+_mu(x+nu) ] U^+_nu(x) ) Lambda_nu(x)
	      //  Staple 5 
	      //    - rho(nu,mu) * Lambda_nu(x+mu)*( [ U_nu(x+mu) U^+_mu(x+nu) ] U^+_nu(x) )
	      //  Staple 6
	      //      rho(mu,nu) *                   [ U_nu(x+mu) U^+_mu(x+nu) ] Lambda_mu(x + nu) U^{+}_nu(x)
	      //
	      //
	      //  Here the suggestive [] are common to all three terms.
	      //  Also Staple 1 and 5 share the additional U^+_nu(x) as indicated by suggestive ()
	      //  Staple 1 and 5 also share the same rho(nu,mu) but have different sign
	      
	      {
		LatticeColorMatrix tmp_mat3;
		LatticeColorMatrix tmp_mat4;
		
		tmp_mat = U_nu_plus_mu*adj(U_mu_plus_nu); // Term in square brackets common to all
		tmp_mat2 = tmp_mat*adj(u[nu]);            // Term in round brackets common to staples 1 and 5
		tmp_mat3 = tmp_mat2*Lambda[nu];           // Staple 1
		tmp_mat4 = Lambda_nu_plus_mu*tmp_mat2;
		tmp_mat3 -= tmp_mat4;                     // Staple 5 and minus sign
		tmp_mat3 *= rho(nu,mu);            // Common factor on staple 1 and 5
		
		tmp_mat4 = shift(Lambda[mu],FORWARD,nu);
		tmp_mat2 = tmp_mat*tmp_mat4;                     // Staple 6
		tmp_mat = tmp_mat2*adj(u[nu]);                   // and again
		tmp_mat *= rho(mu,nu);                    // rho(mu, nu) factor on staple 6
		tmp_mat += tmp_mat3;                          // collect staples 1 5 and 6 onto staple sum
		staple_sum += tmp_mat;                           // slap onto the staple sum
		
		// Tmp 3 disappears here
	      }
	      
	      // The three downward staples
	      //
	      // Paper Staples 2, 3, and 4;
	      //
	      // Staple 2:
	      //     rho_mu_nu * U^{+}_nu(x-nu+mu) [  U^+_mu(x-nu) Lambda_mu(x-nu)     ] U_nu(x-nu)
	      // Staple 3
	      //     rho_nu_mu * U^{+}_nu(x-nu+mu) [ Lambda_nu(x-nu+mu) U^{+}_mu(x-nu) ] U_nu(x-nu)
	      // Staple 4
	      //   - rho_nu_mu * U^{+}_nu(x-nu+mu) [ U^{+}_mu(x-nu) Lambda_nu(x-nu)    ] U_nu(x-nu)
	      //
	      // I have suggestively placed brackets to show that all these staples share a common
	      // first and last term. Secondly staples 3 and 4 share the same value of rho (but opposite sign)
	      //
	      // Finally this can all be communicated on site and then shifted altoghether to x-nu
	      {
		LatticeColorMatrix tmp_mat4;
		
		tmp_mat  = Lambda_nu_plus_mu*adj(u[mu]);      // Staple 3 term in brackets
		tmp_mat4 = adj(u[mu])*Lambda[nu];
		tmp_mat -= tmp_mat4;                          // Staple 4 term in brackets and -ve sign
		tmp_mat *= rho(nu,mu);                        // Staple 3 & 4 common rho value
		tmp_mat2 = adj(u[mu])*Lambda[mu];             // Staple 2 term in brackets
		tmp_mat2 *= rho(mu,nu);                       // Staple 2 rho factor
		tmp_mat += tmp_mat2;                          // Combine terms in brackets, signs and rho factors
		
		tmp_mat2 = adj(U_nu_plus_mu)*tmp_mat;         // Common first matrix
		tmp_mat = tmp_mat2*u[nu];                     // Common last matrix
		
		staple_sum += shift(tmp_mat, BACKWARD, nu);   // Shift it all back to x-nu
	      }
	      
	    } // end of if mu != nu
	  } // end nu loop
	  
	  // Add on this term - there is a relative minus sign which will be corrected by the sign on 
	  // on accumulation to F
	  staple_sum -= adj(C[mu])*Lambda[mu];
	  
	  F[mu] -= timesI(staple_sum);
	  
#if 0
	  QDPIO::cout << __func__ << ":b,  F[" << mu << "]= " << norm2(F[mu]) 
		      << "  staple=" << norm2(staple_sum)
		      << "  Lambda=" << norm2(Lambda[mu]) 
		      << "  C=" << norm2(C[mu]) 
		      << endl;
#endif
	} // End of if(smear_in_this_dirP[mu]
	// Else nothing needs done to the force
      } // end mu loop
      
      
      
      // Done
      END_CODE();
    }
     
    /*! \ingroup gauge */
    void getFs(const LatticeColorMatrix& Q,
	       const LatticeColorMatrix& QQ,
	       multi1d<LatticeComplex>& f)
    {
      START_CODE();

      // Now compute the f's -- use the same function as for computing the fs, bs etc in derivative
      // but don't compute the b-'s
      multi1d<LatticeComplex> b_1; // Dummy - not used      -- throwaway -- won't even get resized
      multi1d<LatticeComplex> b_2; // Dummy - not used here -- throwaway -- won't even get resized
	  
      getFsAndBs(Q,QQ,f,b_1,b_2,false);   // This routine computes the f-s

      END_CODE();
    }


    /* A namespace to hide the thread dispatcher in */
    namespace StoutUtils { 
      struct GetFsAndBsArgs { 
	const LatticeColorMatrix& Q;
	const LatticeColorMatrix& QQ;
	multi1d<LatticeComplex>& f;
	multi1d<LatticeComplex>& b1;
	multi1d<LatticeComplex>& b2;
	bool dobs;
      };



      
    inline 
    void getFsAndBsSiteLoop(int lo, int hi, int myId, 
			      GetFsAndBsArgs* arg)
    {
#ifndef QDP_IS_QDPJIT
      const LatticeColorMatrix& Q = arg->Q;
      const LatticeColorMatrix& QQ = arg->QQ;
      multi1d<LatticeComplex>& f = arg->f;
      multi1d<LatticeComplex>& b1 = arg->b1;
      multi1d<LatticeComplex>& b2 = arg->b2;
      bool dobs=arg->dobs;
      
      for(int site=lo; site < hi; site++)  
      { 
	// Get the traces
	PColorMatrix<QDP::RComplex<REAL>, Nc>  Q_site = Q.elem(site).elem();
	PColorMatrix<QDP::RComplex<REAL>, Nc>  QQ_site = QQ.elem(site).elem();
	PColorMatrix<QDP::RComplex<REAL>, Nc>  QQQ = QQ_site*Q_site;
	
	Real trQQQ; 
	trQQQ.elem()  = realTrace(QQQ);
	Real trQQ;
	trQQ.elem()   = realTrace(QQ_site);
	
	REAL c0    = ((REAL)1/(REAL)3) * trQQQ.elem().elem().elem().elem();  // eq 13
	REAL c1    = ((REAL)1/(REAL)2) * trQQ.elem().elem().elem().elem();	 // eq 15 
	
	
	if( c1 < 4.0e-3  ) 
	{ // RGE: set to 4.0e-3 (CM uses this value). I ran into nans with 1.0e-4
	  // ================================================================================
	  // 
	  // Corner Case 1: if c1 < 1.0e-4 this implies c0max ~ 3x10^-7
	  //    and in this case the division c0/c0max in arccos c0/c0max can be undefined
	  //    and produce NaN's
	  
	  // In this case what we can do is get the f-s a different way. We go back to basics:
	  //
	  // We solve (using maple) the matrix equations using the eigenvalues 
	  //
	  //  [ 1, q_1, q_1^2 ] [ f_0 ]       [ exp( iq_1 ) ]
	  //  [ 1, q_2, q_2^2 ] [ f_1 ]   =   [ exp( iq_2 ) ]
	  //  [ 1, q_3, q_3^2 ] [ f_2 ]       [ exp( iq_3 ) ]
	  //
	  // with q_1 = 2 u w, q_2 = -u + w, q_3 = - u - w
	  // 
	  // with u and w defined as  u = sqrt( c_1/ 3 ) cos (theta/3)
	  //                     and  w = sqrt( c_1 ) sin (theta/3)
	  //                          theta = arccos ( c0 / c0max )
	  // leaving c0max as a symbol.
	  //
	  //  we then expand the resulting f_i as a series around c0 = 0 and c1 = 0
	  //  and then substitute in c0max = 2 ( c_1/ 3)^(3/2)
	  //  
	  //  we then convert the results to polynomials and take the real and imaginary parts:
	  //  we get at the end of the day (to low order)
	  
	  //                  1    2 
	  //   f0[re] := 1 - --- c0  + h.o.t
	  //                 720     
	  //
	  //	         1       1           1        2 
	  //   f0[im] := - - c0 + --- c0 c1 - ---- c0 c1   + h.o.t
	  //               6      120         5040        
	  //
	  //
	  //             1        1            1        2 
	  //   f1[re] := -- c0 - --- c0 c1 + ----- c0 c1  +  h.o.t
	  //             24      360         13440        f
	  //
	  //                 1       1    2    1     3    1     2
	  //   f1[im] := 1 - - c1 + --- c1  - ---- c1  - ---- c0   + h.o.t
	  //                 6      120       5040       5040
	  //
	  //               1   1        1    2     1     3     1     2
	  //   f2[re] := - - + -- c1 - --- c1  + ----- c1  + ----- c0  + h.o.t
	  //               2   24      720       40320       40320    
	  //
	  //              1        1              1        2
	  //   f2[im] := --- c0 - ---- c0 c1 + ------ c0 c1  + h.o.t
	  //             120      2520         120960
	  
	  //  We then express these using Horner's rule for more stable evaluation.
	  // 
	  //  to get the b-s we use the fact that
	  //                                      b2_i = d f_i / d c0
	  //                                 and  b1_i = d f_i / d c1
	  //
	  //  where the derivatives are partial derivativs
	  //
	  //  And we just differentiate the polynomials above (keeping the same level
	  //  of truncation) and reexpress that as Horner's rule
	  // 
	  //  This clearly also handles the case of a unit gauge as no c1, u etc appears in the 
	  //  denominator and the arccos is never taken. In this case, we have the results in 
	  //  the raw c0, c1 form and we don't need to flip signs and take complex conjugates.
	  //
	  //  I checked the expressions below by taking the difference between the Horner forms
	  //  below from the expanded forms (and their derivatives) above and checking for the
	  //  differences to be zero. At this point in time maple seems happy.
	  //  ==================================================================================
	  
	  f[0].elem(site).elem().elem().real() = 1.0-c0*c0/720.0;
	  f[0].elem(site).elem().elem().imag() =  -(c0/6.0)*(1.0-(c1/20.0)*(1.0-(c1/42.0))) ;
	  
	  f[1].elem(site).elem().elem().real() =  c0/24.0*(1.0-c1/15.0*(1.0-3.0*c1/112.0)) ;
	  f[1].elem(site).elem().elem().imag() =  1.0-c1/6.0*(1.0-c1/20.0*(1.0-c1/42.0))-c0*c0/5040.0 ;
	  
	  f[2].elem(site).elem().elem().real() = 0.5*(-1.0+c1/12.0*(1.0-c1/30.0*(1.0-c1/56.0))+c0*c0/20160.0);
	  f[2].elem(site).elem().elem().imag() = 0.5*(c0/60.0*(1.0-c1/21.0*(1.0-c1/48.0)));
	  
	  if( dobs == true ) {
	    //  partial f0/ partial c0
	    b2[0].elem(site).elem().elem().real() = -c0/360.0;
	    b2[0].elem(site).elem().elem().imag() =  -(1.0/6.0)*(1.0-(c1/20.0)*(1.0-c1/42.0));
	    
	    // partial f0 / partial c1
	    //
	    b1[0].elem(site).elem().elem().real() = 0;
	    b1[0].elem(site).elem().elem().imag() = (c0/120.0)*(1.0-c1/21.0);
	    
	    // partial f1 / partial c0
	    //
	    b2[1].elem(site).elem().elem().real() = (1.0/24.0)*(1.0-c1/15.0*(1.0-3.0*c1/112.0));
	    b2[1].elem(site).elem().elem().imag() = -c0/2520.0;
	    
	    
	    // partial f1 / partial c1
	    b1[1].elem(site).elem().elem().real() = -c0/360.0*(1.0 - 3.0*c1/56.0 );
	    b1[1].elem(site).elem().elem().imag() = -1.0/6.0*(1.0-c1/10.0*(1.0-c1/28.0));
	    
	    // partial f2/ partial c0
	    b2[2].elem(site).elem().elem().real() = 0.5*c0/10080.0;
	    b2[2].elem(site).elem().elem().imag() = 0.5*(  1.0/60.0*(1.0-c1/21.0*(1.0-c1/48.0)) );
	    
	    // partial f2/ partial c1
	    b1[2].elem(site).elem().elem().real() = 0.5*(  1.0/12.0*(1.0-(2.0*c1/30.0)*(1.0-3.0*c1/112.0)) ); 
	    b1[2].elem(site).elem().elem().imag() = 0.5*( -c0/1260.0*(1.0-c1/24.0) );
	    
#if 0
	    {
	      multi1d<int> coord = Layout::siteCoords(Layout::nodeNumber(), site);
	      
	      QMP_fprintf(stdout, 
			  "%s: corner; site=%d coord=[%d,%d,%d,%d] f[0]=%g f[1]=%g f[2]=%g b1[0]=%g b1[1]=%g b1[2]=%g b2[0]=%g b2[1]=%g b2[2]=%g c0=%g c1=%g",
			  
			  __func__, site, coord[0], coord[1], coord[2], coord[3],
			  toDouble(localNorm2(f[0].elem(site))),
			  toDouble(localNorm2(f[1].elem(site))),
			  toDouble(localNorm2(f[2].elem(site))),
			  toDouble(localNorm2(b1[0].elem(site))),
			  toDouble(localNorm2(b1[1].elem(site))),
			  toDouble(localNorm2(b1[2].elem(site))),
			  toDouble(localNorm2(b2[0].elem(site))),
			  toDouble(localNorm2(b2[1].elem(site))),
			  toDouble(localNorm2(b2[2].elem(site))),
			  c0, c1);
	    }
#endif
	  } // Dobs==true
	}
	else 
	{ 
	  // ===================================================================================
	  // Normal case: Do as per paper
	  // ===================================================================================
	  bool c0_negativeP = c0 < 0;
	  REAL c0abs = fabs((double)c0);
	  REAL c0max = 2*pow( (double)(c1/(double)3), (double)1.5);
	  REAL theta;
	  
	  // ======================================================================================
	  // Now work out theta. In the paper the case where c0 -> c0max even when c1 is reasonable 
	  // Has never been considered, even though it can arise and can cause the arccos function
	  // to fail
	  // Here we handle it with series expansion
	  // =====================================================================================
	  REAL eps = (c0max - c0abs)/c0max;
	  
	  if( eps < 0 ) {
	    // ===============================================================================
	    // Corner Case 2: Handle case when c0abs is bigger than c0max. 
	    // This can happen only when there is a rounding error in the ratio, and that the 
	    // ratio is really 1. This implies theta = 0 which we'll just set.
	    // ===============================================================================
	    {
	      multi1d<int> coord = Layout::siteCoords(Layout::nodeNumber(), site);

//	      fprintf(stdout, 
//		      "%s: corner2; site=%d coord=[%d,%d,%d,%d] c0max=%g c0abs=%d eps=%g\n Setting theta=0",
//		      __func__, site, coord[0], coord[1], coord[2], coord[3],
//		      c0abs, c0max,eps);
	    }
	    theta = 0;
	  }
	  else if ( eps < 1.0e-3 ) {
	    // ===============================================================================
	    // Corner Case 3: c0->c0max even though c1 may be actually quite reasonable.
	    // The ratio |c0|/c0max -> 1 but is still less than one, so that a 
	    // series expansion is possible.
	    // SERIES of acos(1-epsilon): Good to O(eps^6) or with this cutoff to O(10^{-18}) Computed with Maple.
	    //  BTW: 1-epsilon = 1 - (c0max-c0abs)/c0max = 1-(1 - c0abs/c0max) = +c0abs/c0max
	    //
	    // ===============================================================================
	    REAL sqtwo = sqrt((REAL)2);
	    
	    theta = sqtwo*sqrt(eps)*( 1.0 + ( (1/(REAL)12) + ( (3/(REAL)160) + ( (5/(REAL)896) + ( (35/(REAL)18432) + (63/(REAL)90112)*eps ) *eps) *eps) *eps) *eps);
	    
	  } 
	  else {  
	    // 
	    theta = acos( c0abs/c0max );
	  }
	  
	  multi1d<REAL> f_site_re(3);
	  multi1d<REAL> f_site_im(3);
	  
	  multi1d<REAL> b1_site_re(3);
	  multi1d<REAL> b1_site_im(3);
	  
	  multi1d<REAL> b2_site_re(3);
	  multi1d<REAL> b2_site_im(3);
	  
	  
	  
	  REAL u = sqrt(c1/3)*cos(theta/3);
	  REAL w = sqrt(c1)*sin(theta/3);
	  
	  REAL u_sq = u*u;
	  REAL w_sq = w*w;
	  
	  REAL xi0,xi1;
	  {
	    bool w_smallP  = fabs(w) < 0.05;
	    if( w_smallP ) { 
	      xi0 = (REAL)1 - ((REAL)1/(REAL)6)*w_sq*( 1 - ((REAL)1/(REAL)20)*w_sq*( (REAL)1 - ((REAL)1/(REAL)42)*w_sq ) );
	    }
	    else {
	      xi0 = sin(w)/w;
	    }
	    
	    if( dobs==true) {
	      
	      if( w_smallP  ) { 
		xi1 = -1*( ((REAL)1/(REAL)3) - ((REAL)1/(REAL)30)*w_sq*( (REAL)1 - ((REAL)1/(REAL)28)*w_sq*( (REAL)1 - ((REAL)1/(REAL)54)*w_sq ) ) );
	      }
	      else { 
		xi1 = cos(w)/w_sq - sin(w)/(w_sq*w);
	      }
	    }
	  }
	  
	  REAL cosu = cos(u);
	  REAL sinu = sin(u);
	  REAL cosw = cos(w);
	  REAL sinw = sin(w);
	  REAL sin2u = sin(2*u);
	  REAL cos2u = cos(2*u);
	  REAL ucosu = u*cosu;
	  REAL usinu = u*sinu;
	  REAL ucos2u = u*cos2u;
	  REAL usin2u = u*sin2u;
	  
	  REAL denum = (REAL)9*u_sq - w_sq;
	  
	  {
	    REAL subexp1 = u_sq - w_sq;
	    REAL subexp2 = 8*u_sq*cosw;
	    REAL subexp3 = (3*u_sq + w_sq)*xi0;
	    
	    f_site_re[0] = ( (subexp1)*cos2u + cosu*subexp2 + 2*usinu*subexp3 ) / denum ;
	    f_site_im[0] = ( (subexp1)*sin2u - sinu*subexp2 + 2*ucosu*subexp3 ) / denum ;
	  }
	  {
	    REAL subexp = (3*u_sq -w_sq)*xi0;
	    
	    f_site_re[1] = (2*(ucos2u - ucosu*cosw)+subexp*sinu)/denum;
	    f_site_im[1] = (2*(usin2u + usinu*cosw)+subexp*cosu)/denum;
	  }
	  
	  
	  {
	    REAL subexp=3*xi0;
	    
	    f_site_re[2] = (cos2u - cosu*cosw -usinu*subexp) /denum ;
	    f_site_im[2] = (sin2u + sinu*cosw -ucosu*subexp) /denum ;
	  }
	  
	  if( dobs == true ) 
	    {
	      multi1d<REAL> r_1_re(3);
	      multi1d<REAL> r_1_im(3);
	      multi1d<REAL> r_2_re(3);
	      multi1d<REAL> r_2_im(3);
	      
	      //	  r_1[0]=Double(2)*cmplx(u, u_sq-w_sq)*exp2iu
	      //          + 2.0*expmiu*( cmplx(8.0*u*cosw, -4.0*u_sq*cosw)
	      //	      + cmplx(u*(3.0*u_sq+w_sq),9.0*u_sq+w_sq)*xi0 );
	      {
		REAL subexp1 = u_sq - w_sq;
		REAL subexp2 =  8*cosw + (3*u_sq + w_sq)*xi0 ;
		REAL subexp3 =  4*u_sq*cosw - (9*u_sq + w_sq)*xi0 ;
		
		r_1_re[0] = 2*(ucos2u - sin2u *(subexp1)+ucosu*( subexp2 )- sinu*( subexp3 ) );
		r_1_im[0] = 2*(usin2u + cos2u *(subexp1)-usinu*( subexp2 )- cosu*( subexp3 ) );
	      }
	      
	      // r_1[1]=cmplx(2.0, 4.0*u)*exp2iu + expmiu*cmplx(-2.0*cosw-(w_sq-3.0*u_sq)*xi0,2.0*u*cosw+6.0*u*xi0);
	      {
		REAL subexp1 = cosw+3*xi0;
		REAL subexp2 = 2*cosw + xi0*(w_sq - 3*u_sq);
		
		r_1_re[1] = 2*((cos2u - 2*usin2u) + usinu*( subexp1 )) - cosu*( subexp2 );
		r_1_im[1] = 2*((sin2u + 2*ucos2u) + ucosu*( subexp1 )) + sinu*( subexp2 );
	      }
	      
	      
	      // r_1[2]=2.0*timesI(exp2iu)  +expmiu*cmplx(-3.0*u*xi0, cosw-3*xi0);
	      {
		REAL subexp = cosw - 3*xi0;
		r_1_re[2] = -2*sin2u -3*ucosu*xi0 + sinu*( subexp );
		r_1_im[2] = 2*cos2u  +3*usinu*xi0 + cosu*( subexp );
	      }
	      
	      
	      //r_2[0]=-2.0*exp2iu + 2*cmplx(0,u)*expmiu*cmplx(cosw+xi0+3*u_sq*xi1,
	      //						 4*u*xi0);
	      {
		REAL subexp = cosw + xi0 + 3*u_sq*xi1;
		r_2_re[0] = -2*(cos2u + u*( 4*ucosu*xi0 - sinu*(subexp )) );
		r_2_im[0] = -2*(sin2u - u*( 4*usinu*xi0 + cosu*(subexp )) );
	      }
	      
	      
	      // r_2[1]= expmiu*cmplx(cosw+xi0-3.0*u_sq*xi1, 2.0*u*xi0);
	      // r_2[1] = timesMinusI(r_2[1]);
	      {
		REAL subexp =  cosw + xi0 - 3*u_sq*xi1;
		r_2_re[1] =  2*ucosu*xi0 - sinu*( subexp ) ;
		r_2_im[1] = -2*usinu*xi0 - cosu*( subexp ) ;
	      }
	      
	      //r_2[2]=expmiu*cmplx(xi0, -3.0*u*xi1);
	      {
		REAL subexp = 3*xi1;
		
		r_2_re[2] =    cosu*xi0 - usinu*subexp ;
		r_2_im[2] = -( sinu*xi0 + ucosu*subexp ) ;
	      }      
	      
	      REAL b_denum=2*denum*denum;
	      
	      
	      for(int j=0; j < 3; j++) { 
		
		{
		  REAL subexp1 = 2*u;
		  REAL subexp2 = 3*u_sq - w_sq;
		  REAL subexp3 = 2*(15*u_sq + w_sq);
		  
		  b1_site_re[j]=( subexp1*r_1_re[j] + subexp2*r_2_re[j] - subexp3*f_site_re[j] )/b_denum;
		  b1_site_im[j]=( subexp1*r_1_im[j] + subexp2*r_2_im[j] - subexp3*f_site_im[j] )/b_denum;
		}
		
		{ 
		  REAL subexp1 = 3*u;
		  REAL subexp2 = 24*u;
		  
		  b2_site_re[j]=( r_1_re[j]- subexp1*r_2_re[j] - subexp2 * f_site_re[j] )/b_denum;
		  b2_site_im[j]=( r_1_im[j] -subexp1*r_2_im[j] - subexp2 * f_site_im[j] )/b_denum;
		}
	      }

	      // Now flip the coefficients of the b-s
	      if( c0_negativeP ) 
	      {
		//b1_site[0] = conj(b1_site[0]);
		b1_site_im[0] *= -1;
		
		//b1_site[1] = -conj(b1_site[1]);
		b1_site_re[1] *= -1;
		
		//b1_site[2] = conj(b1_site[2]);
		b1_site_im[2] *= -1;
		
		//b2_site[0] = -conj(b2_site[0]);
		b2_site_re[0] *= -1;
		
		//b2_site[1] = conj(b2_site[1]);
		b2_site_im[1] *= -1;
		
		//b2_site[2] = -conj(b2_site[2]);
		b2_site_re[2] *= -1;
	      }
	      
	      // Load back into the lattice sized object
	      for(int j=0; j < 3; j++) {
		
		b1[j].elem(site).elem().elem().real() = b1_site_re[j];
		b1[j].elem(site).elem().elem().imag() = b1_site_im[j];
		
		b2[j].elem(site).elem().elem().real() = b2_site_re[j];
		b2[j].elem(site).elem().elem().imag() = b2_site_im[j];
	      }
#if 0
	      {
		multi1d<int> coord = Layout::siteCoords(Layout::nodeNumber(), site);
		REAL rat = c0abs/c0max;
		
		QMP_fprintf(stdout, 
			    "%s: site=%d coord=[%d,%d,%d,%d] f_site[0]=%g f_site[1]=%g f_site[2]=%g 1[0]=%g b1[1]=%g b1[2]=%g b2[0]=%g b2[1]=%g b2[2]=%g denum=%g c0=%g c1=%g c0max=%g rat=%g theta=%g",
			    
			    __func__, site, coord[0], coord[1], coord[2], coord[3],
			    toDouble(localNorm2(cmplx(Real(f_site_re[0]),Real(f_site_im[0])))),
			    toDouble(localNorm2(cmplx(Real(f_site_re[1]),Real(f_site_im[1])))),
			    toDouble(localNorm2(cmplx(Real(f_site_re[2]),Real(f_site_im[2])))),
			    toDouble(localNorm2(b1[0].elem(site))),
			    toDouble(localNorm2(b1[1].elem(site))),
			    toDouble(localNorm2(b1[2].elem(site))),
			    toDouble(localNorm2(b2[0].elem(site))),
			    toDouble(localNorm2(b2[1].elem(site))),
			    toDouble(localNorm2(b2[2].elem(site))),
			    denum, 
			    c0, c1, c0max,
			    rat, theta);
	      }
#endif
	      
	    } // end of if (dobs==true)
	  
	  // Now when everything is done flip signs of the b-s (can't do this before
	  // as the unflipped f-s are needed to find the b-s
	  
	  if( c0_negativeP ) {
	    
	    // f_site[0] = conj(f_site[0]);
	    f_site_im[0] *= -1;
	    
	    //f_site[1] = -conj(f_site[1]);
	    f_site_re[1] *= -1;
	    
	    //f_site[2] = conj(f_site[2]);
	    f_site_im[2] *= -1;
	    
	  }
	
	  // Load back into the lattice sized object
	  for(int j=0; j < 3; j++) { 
	    f[j].elem(site).elem().elem().real() = f_site_re[j];
	    f[j].elem(site).elem().elem().imag() = f_site_im[j];
	  }
	  
	} // End of if( corner_caseP ) else {}
      } // End site loop
#endif
    } // End Function

    } // End Namespace
	
    /*! \ingroup gauge */
    void getFsAndBs(const LatticeColorMatrix& Q,
		    const LatticeColorMatrix& QQ,
		    multi1d<LatticeComplex>& f,
		    multi1d<LatticeComplex>& b1,
		    multi1d<LatticeComplex>& b2,
		    bool dobs)
    {
      START_CODE();
      QDP::StopWatch swatch;
      swatch.reset();
      swatch.start();
      
      f.resize(3);

      b1.resize(3);
      b2.resize(3);

      
      int num_sites = Layout::sitesOnNode();
      StoutUtils::GetFsAndBsArgs args={Q,QQ,f,b1,b2,dobs};

#if defined(BUILD_JIT_CLOVER_TERM)
      //QDPIO::cout << "PTX getFsAndBs dobs = " << dobs << "\n";
      static CUfunction function;
      
      if (function == NULL)
	function = function_get_fs_bs_build( Q,QQ,f,b1,b2,dobs );
      
      // Execute the function
      function_get_fs_bs_exec(function, Q,QQ,f,b1,b2,dobs );
#else
      dispatch_to_threads(num_sites, args, StoutUtils::getFsAndBsSiteLoop);
#endif

      swatch.stop();
      StoutLinkTimings::functions_secs += swatch.getTimeInSeconds();
      END_CODE();
    }

    /*! \ingroup gauge */
    void smear_links(const multi1d<LatticeColorMatrix>& current, 
		     multi1d<LatticeColorMatrix>& next,
		     const multi1d<bool>& smear_in_this_dirP,
		     const multi2d<Real>& rho)
    {
      START_CODE();
      
      for(int mu = 0; mu < Nd; mu++) 
      {
	if( smear_in_this_dirP[mu] ) 
	{
	  LatticeColorMatrix Q, QQ;
	  
	  // Q contains the staple term. C is a throwaway
	  getQs(current, Q, QQ, mu, smear_in_this_dirP, rho);
	  
	  // Now compute the f's
	  multi1d<LatticeComplex> f;   // routine will resize these
	  getFs(Q,QQ,f);   // This routine computes the f-s
	  
	  // Assemble the stout links exp(iQ)U_{mu} 
	  next[mu]=(f[0] + f[1]*Q + f[2]*QQ)*current[mu];      
	}
	else { 
	  next[mu]=current[mu];  // Unsmeared
	}
	
      }
      
      END_CODE();
    }
    

    /*! \ingroup gauge */
    void stout_smear(LatticeColorMatrix& next,
		     const multi1d<LatticeColorMatrix>& current, 
		     int mu,
		     const multi1d<bool>& smear_in_this_dirP,
		     const multi2d<Real>& rho)
    {
      START_CODE();
      
      LatticeColorMatrix Q, QQ;
	  
      // Q contains the staple term. C is a throwaway
      getQs(current, Q, QQ, mu, smear_in_this_dirP, rho);
	  
      // Now compute the f's
      multi1d<LatticeComplex> f;   // routine will resize these
      getFs(Q,QQ,f);   // This routine computes the f-s
	  
      // Assemble the stout links exp(iQ)U_{mu} 
      next = (f[0] + f[1]*Q + f[2]*QQ)*current[mu];      
      
      END_CODE();
    }
    
  }

}
