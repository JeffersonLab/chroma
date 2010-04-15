// -*- C++ -*-
// $Id: reunit.cc,v 3.3 2007-02-22 21:11:49 bjoo Exp $

/*! \file
 *  \ingroup gauge
 *  \author Subsetting added by A. Hart
 *  \param[in,out] xa  The descriptor of matrices to be reunitarized.
 *            Must be of type LatticeColorMatrix
 *  \param[in] bad Descriptor of flags indicating sites violating unitarity.
 *            Only used if ruflag = REUNITARIZE_LABEL or
 *            REUNITARIZE_ERROR.
 *  \param[in] ruflag Can also be REUNITARIZE in which case the
 *            matrices are reunitarized but no complaints are made.
 *  \param[out] numbad Total number of matrices violating unitarity.
 *            ONLY USED IF ruflag is testing for ERROR or LABEL. 
 *  \param[in] mstag  An (un)ordered subset of sites
 *  \brief Reunitarize in place a color matrix to SU(N)
 *
 *  Reunitarize (to a SU(N)) inplace the matrix XA under some option
 */

#include "chromabase.h"
#include "util/gauge/reunit.h"


namespace Chroma {


  namespace ReunitEnv {
    static double time_spent =0;
    double getTime() { return time_spent; }
  };

  template<typename Q, typename C, typename R, typename S>
  inline
  void reunit_t(Q& xa, 
		LatticeBoolean& bad, 
		int& numbad, 
		enum Reunitarize ruflag,
		const S& mstag)
  {
    START_CODE();
    
    QDP::StopWatch swatch;
    swatch.reset();
    swatch.start();

    multi2d<C> a(Nc, Nc);
    multi2d<C> b(Nc, Nc);
    R t1;
    C t2;
    R t3;
    R t4;
    R sigmasq = 0;
    multi1d<C> row(Nc);
    
    // The initial number of matrices violating unitarity.
    numbad = 0;
    
    // some kind of small floating point number, should be prec. dep.
    R fuzz = 1.0e-5; 
    
    // Extract initial components 
    for(int i=0; i < Nc; ++i)
      for(int j=0; j < Nc; ++j)
	(a[i][j])[mstag] = peekColor(xa, i, j);
    
    // Use the Nc-dependent reunitarizers
    switch (Nc)
      {
	//---------------------------------------------------------
	// U(1)
      case 1:
	/* normalise the complex number */
	/* sigmasq = sqrt(u^t . u) */
	sigmasq[mstag] = sqrt(localNorm2(a[0][0]));
	
	/* rescale the field */
	/* u <- u/sigmasq */
	(a[0][0])[mstag] /= sigmasq;
	
	/* Now, do various things depending on the input flag. */
	/* For use later, calculate the mean squared deviation */
	if ( ruflag == REUNITARIZE_ERROR ||
	     ruflag == REUNITARIZE_LABEL )
	  {
	    sigmasq[mstag] = fabs(1-sigmasq);
	  }
	
	/* Do things depending on the mean squared deviation */
	switch (ruflag)
	  {
	  case REUNITARIZE_ERROR:
	    /* Gripe and stop if unitarity is violated. */
	    numbad = toInt(sum(where(sigmasq > fuzz, LatticeInteger(1), LatticeInteger(0))));
	    if ( numbad > 0 )
	      QDP_error_exit("Unitarity violated", numbad);
	    break;
	    
	  case REUNITARIZE_LABEL:
	    /* Label the bad guys if unitarity is violated. */
	    bad[mstag] = sigmasq > fuzz;
	    numbad = toInt( sum(where(bad, LatticeInteger(1), LatticeInteger(0)), 
				mstag) );
	    break;
	  default:
	    break;
	  }
	break;
	
	//---------------------------------------------------------
	// SU(2)
      case 2:
	/* If you want to check unitarity, I have to save the second row somewhere */
	if ( ruflag == REUNITARIZE_ERROR  ||
	     ruflag == REUNITARIZE_LABEL )
	  {
	    for(int c = 0; c < Nc; ++c)
	      (row[c])[mstag] = a[c][Nc-1];
	  }
	
	/* normalise the first row */
	/* t1 = sqrt(u^t . u) */
	t1[mstag] = localNorm2(a[0][0]);
	for(int c = 1; c < Nc; ++c)
	  t1[mstag] += localNorm2(a[c][0]);
	t1[mstag] = sqrt(t1);
	
	
	/* overwrite the first row with the rescaled value */
	/* u <- u/t1 */
	t4[mstag] = 1 / t1;
	for(int c = 0; c < Nc; ++c)
	  (a[c][0])[mstag] *= t4;
	
	/* construct the second row from the first row */
	(a[1][1])[mstag] = adj(a[0][0]);
	(a[0][1])[mstag] = -adj(a[1][0]);
	
	/* Now, do various things depending on the input flag. */
	/* For use later, calculate the mean squared deviation */
	if ( ruflag == REUNITARIZE_ERROR ||
	     ruflag == REUNITARIZE_LABEL )
	  {
	    /* Calculate the mean square deviation.  */
	    /* sigmasq = (1 - t1)**2 */
	    sigmasq[mstag] = pow(1-t1,2);
	    
	    /* sigmasq <- sqrt(sigmasq + |crow(.)-a(1,.)|**2) */
	    /* overwrite row */
	    for(int c = 0; c < Nc; ++c)
	      sigmasq[mstag] += localNorm2(row[c] - a[c][Nc-1]);
	    
	    sigmasq[mstag] = sqrt(sigmasq);
	  }
	
	
	/* Do things depending on the mean squared deviation */
	switch (ruflag)
	  {
	  case REUNITARIZE_ERROR:
	    /* Gripe and stop if unitarity is violated. */
	    numbad = toInt(sum(where(sigmasq > fuzz, LatticeInteger(1), LatticeInteger(0))));
	    if ( numbad > 0 )
	      QDP_error_exit("Unitarity violated", numbad);
	    break;
	    
	  case REUNITARIZE_LABEL:
	    /* Label the bad guys if unitarity is violated. */
	    bad[mstag] = sigmasq > fuzz;
	    numbad = toInt(sum(where(bad, LatticeInteger(1), LatticeInteger(0)),
			       mstag));
	    break;
	    
	  default:
	    break;
	  }
	break;
	
	//---------------------------------------------------------
	// SU(3)
      case 3:
	/* If you want to check unitarity, I have to save the third row somewhere */
	if ( ruflag == REUNITARIZE_ERROR  ||
	     ruflag == REUNITARIZE_LABEL )
	  {
	    for(int c = 0; c < Nc; ++c)
	      (row[c])[mstag] = a[c][Nc-1];
	  }
	
	
	/* normalise the first row */
	/* t1 = sqrt(u^t . u) */
	t1[mstag] = localNorm2(a[0][0]);
	for(int c = 1; c < Nc; ++c)
	  t1[mstag] += localNorm2(a[c][0]);
	t1[mstag] = sqrt(t1);
	
	
	/* overwrite the first row with the rescaled value */
	/* u <- u/t1 */
	t4[mstag] = 1 / t1;
	for(int c = 0; c < Nc; ++c)
	  (a[c][0])[mstag] *= t4;
	
	/* calculate the orthogonal component to the second row */
	/* t2 <- u^t.v */
	t2[mstag] = adj(a[0][0]) * a[0][1];
	for(int c = 1; c < Nc; ++c)
	  t2[mstag] += adj(a[c][0]) * a[c][1];
	
	/* orthogonalize the second row relative to the first row */
	/* v <- v - t2*u */
	for(int c = 0; c < Nc; ++c)
	  (a[c][1])[mstag] -= t2 * a[c][0]; 
	
	/* normalise the second row */
	/* t3 = sqrt(u^t . u) */
	t3[mstag] = localNorm2(a[0][1]);
	for(int c = 1; c < Nc; ++c)
	  t3[mstag] += localNorm2(a[c][1]);
	t3[mstag] = sqrt(t3);
	
	/* overwrite the second row with the rescaled value */
	/* v <- v/t3 */
	t4[mstag] = 1 / t3;
	for(int c = 0; c < Nc; ++c)
	  (a[c][1])[mstag] *= t4;
	
	
	/* the third row is the cross product of the new first and second rows */
	/* column 1: w(0) = u(1)*v(2) - u(2)*v(1) */
	(a[0][2])[mstag] = adj(a[1][0]) * adj(a[2][1]) - 
	  adj(a[2][0]) * adj(a[1][1]);
	
	/* column 2: w(1) = u(2)*v(0) - u(0)*v(2) */
	(a[1][2])[mstag] = adj(a[2][0]) * adj(a[0][1]) - 
	  adj(a[0][0]) * adj(a[2][1]);
	
	/* column 3: w(3) = u(1)*v(2) - u(2)*v(1) */
	(a[2][2])[mstag] = adj(a[0][0]) * adj(a[1][1]) - 
	  adj(a[1][0]) * adj(a[0][1]);
	
	/* Now, do various things depending on the input flag. */
	/* For use later, calculate the mean squared deviation */
	if ( ruflag == REUNITARIZE_ERROR ||
	     ruflag == REUNITARIZE_LABEL )
	  {
	    /* Calculate the mean square deviation.  */
	    /* sigmasq = (1 - t1)**2 + |t2|**2 + (1 - t3)**2 */
	    sigmasq[mstag] = pow(1-t1,2) + localNorm2(t2) + pow(1-t3,2);
	    
	    /* sigmasq <- sqrt(sigmasq + |crow(.)-a(2,.)|**2) */
	    /* overwrite row */
	    for(int c = 0; c < Nc; ++c)
	      sigmasq[mstag] += localNorm2(row[c] - a[c][Nc-1]);
	    
	    sigmasq[mstag] = sqrt(sigmasq);
	  }
	
	
	/* Do things depending on the mean squared deviation */
	switch (ruflag)
	  {
	  case REUNITARIZE_ERROR:
	    /* Gripe and stop if unitarity is violated. */
	    numbad = toInt(sum(where(sigmasq > fuzz, LatticeInteger(1), 
				     LatticeInteger(0)),
			       mstag));
	    if ( numbad > 0 )
	      QDP_error_exit("Unitarity violated", numbad);
	    break;
	  case REUNITARIZE_LABEL:
	    /* Label the bad guys if unitarity is violated. */
	    bad[mstag] = sigmasq > fuzz;
	    numbad = toInt(sum(where(bad,LatticeInteger(1), LatticeInteger(0)),
			       mstag));
	    
	  default:
	    break;
	  }
	break;
	
	//---------------------------------------------------------
	// SU(N > 3)
      default:
	if ( Nc > 3 )
	  {
	    /* If you want to check unitarity, I have to save the third row somewhere */
	    if ( ruflag == REUNITARIZE_ERROR  ||
		 ruflag == REUNITARIZE_LABEL )
	      {
		for(int c = 0; c < Nc; ++c)
		  (row[c])[mstag] = a[c][Nc-1];
	      }
	    
	    /* normalise the first row */
	    /* t1 = sqrt(u^t . u) */
	    t1[mstag] = localNorm2(a[0][0]);
	    for(int c = 1; c < Nc; ++c)
	      t1[mstag] += localNorm2(a[c][0]);
	    t1[mstag] = sqrt(t1);
	    
	    if ( ruflag == REUNITARIZE_ERROR ||
		 ruflag == REUNITARIZE_LABEL )
	      {
		/* Calculate the mean square deviation.  */
		/* sigmasq = (1 - t1)**2 */
		sigmasq[mstag] = pow(1-t1,2);
	      }
	    
	    /* overwrite the first row with the rescaled value */
	    /* u <- u/t1 */
	    t3[mstag] = 1 / t1;
	    for(int c = 0; c < Nc; ++c)
	      (a[c][0])[mstag] *= t3;
	    
	    /* Do Gramm-Schmidt on the remaining rows */
	    for(int j = 1; j < Nc; j++ )
	      {
		for(int i = 0; i < j; i++ )
		  {
		    /* t2 <- u^t.v */
		    t2[mstag] = adj(a[0][i]) * a[0][j];
		    for(int c = 1; c < Nc; ++c)
		      {
			t2[mstag] += adj(a[c][i]) * a[c][j];
		      }
		    
		    if ( (ruflag == REUNITARIZE_ERROR ||
			  ruflag == REUNITARIZE_LABEL) && j < (Nc-1) )
		      sigmasq[mstag] += localNorm2(t2);
		    
		    /* orthogonalize the j-th row relative to the i-th row */
		    /* v <- v - t2*u */
		    for(int c = 0; c < Nc; ++c)
		      {
			(a[c][j])[mstag] -= t2 * a[c][i];
		      }
		  }
		
		/* normalise the j-th row */
		/* t1 = sqrt(v^t . v) */
		t1[mstag] = localNorm2(a[0][j]);
		for(int c = 1; c < Nc; ++c)
		  t1[mstag] += localNorm2(a[c][j]);
		t1[mstag] = sqrt(t1);
		
		/* overwrite the j-th row with the rescaled value */
		/* v <- v/t1 */
		t3[mstag] = 1 / t1;
		for(int c = 0; c < Nc; ++c)
		  (a[c][j])[mstag] *= t3;
		
		if ( (ruflag == REUNITARIZE_ERROR ||
		      ruflag == REUNITARIZE_LABEL) && j < (Nc-1) )
		  {
		    /* Calculate the mean square deviation.  */
		    /* sigmasq = (1 - t1)**2 */
		    sigmasq[mstag] += pow(1-t1,2);
		  }
	      }
	    
	    /* Now we have a unitary matrix. We need to multiply the last
	       row with a phase to make the determinant 1. */
	    /* We compute the determinant by LU decomposition */
	    for(int j = 0; j < Nc; j++)
	      for(int i = 0; i < Nc; i++)
		(b[j][i])[mstag] = a[j][i];
	    
	    for(int j = 0; j < Nc; j++)
	      {
		for(int i = 0; i <= j; i++)
		  {
		    t2[mstag] = b[j][i];
		    for(int c = 0; c < i; c++)
		      t2[mstag] -= b[c][i] * b[j][c];
		    
		    (b[j][i])[mstag] = t2;
		  }
		
		for(int i = (j+1); i < Nc; i++)
		  {
		    t2[mstag] = b[j][i];
		    for(int c = 0; c < j; c++)
		      t2[mstag] -= b[c][i] * b[j][c];
		    
		    (b[j][i])[mstag] = adj(b[j][j]) * t2 / localNorm2(b[j][j]);
		  }
	      }
	    
	    /* The determinant */
	    t2[mstag] = b[0][0] * b[1][1];
	    for(int c = 2; c < Nc; c++)
	      t2[mstag] *= b[c][c];
	    
	    /* The phase of the determinant */
	    t4[mstag] = atan2(imag(t2), real(t1));
	    t2[mstag] = cmplx(cos(t4), -sin(t4));
	    for(int c = 0; c < Nc; ++c)
	      (a[c][Nc-1])[mstag] *= t2;
	    
            
	    /* Now, do various things depending on the input flag. */
	    /* For use later, finish calculating the mean squared deviation */
	    if ( ruflag == REUNITARIZE_ERROR ||
		 ruflag == REUNITARIZE_LABEL )
	      {
		for(int c = 0; c < Nc; ++c)
		  sigmasq[mstag] += localNorm2(row[c] - a[c][Nc-1]);
		
		sigmasq[mstag] = sqrt(sigmasq);
	      }
	    
	    numbad = toInt(sum(where(sigmasq > fuzz, LatticeInteger(1), 
				     LatticeInteger(0))));
	    
	    /* Do things depending on the mean squared deviation */
	    switch (ruflag)
	      {
	      case REUNITARIZE_ERROR:
		/* Gripe and stop if unitarity is violated. */
		numbad = toInt(sum(where(sigmasq > fuzz, LatticeInteger(1), 
					 LatticeInteger(0))));
		if ( numbad > 0 )
		  QDP_error_exit("Unitarity violated", numbad);
		break;
	      case REUNITARIZE_LABEL:
		/* Label the bad guys if unitarity is violated. */
		bad[mstag] = sigmasq > fuzz;
		numbad = toInt(sum(where(bad, LatticeInteger(1), 
					 LatticeInteger(0)),
				   mstag));
	      default:
		break;
	      }
	  }
	else
	  QDP_error_exit("Invalid Nc for reunit, Nc=%d", Nc);
      }
    
    // Insert final reunitarized components 
    for(int i=0; i < Nc; ++i)
      for(int j=0; j < Nc; ++j)
	pokeColor(xa[mstag], a[i][j], i, j);
    
    swatch.stop();
    ReunitEnv::time_spent += swatch.getTimeInSeconds();
    END_CODE();
  }
  
  // Overloaded definitions
  // SINGLE
  void reunit(LatticeColorMatrixF3& xa)
  {
    START_CODE();

    LatticeBoolean bad;
    int numbad;
    
    reunit_t<LatticeColorMatrixF3, LatticeComplexF, LatticeRealF, Subset>(xa, bad, numbad, REUNITARIZE, all);
    
    END_CODE();
  }

   // Overloaded definitions
  // DOUBLE
  void reunit(LatticeColorMatrixD3& xa)
  {
    START_CODE();

    LatticeBoolean bad;
    int numbad;
    
    reunit_t<LatticeColorMatrixD3, LatticeComplexD, LatticeRealD, Subset>(xa, bad, numbad, REUNITARIZE, all);
    
    END_CODE();
  }


  // SINGLE
  void reunit(LatticeColorMatrixF3& xa,
	      const Subset& mstag)
  {
    START_CODE();

    LatticeBoolean bad;
    int numbad;
    
    reunit_t<LatticeColorMatrixF3, LatticeComplexF, LatticeRealF, Subset>(xa, bad, numbad, REUNITARIZE, mstag);
    
    END_CODE();
  }

  // DOUBLE
  void reunit(LatticeColorMatrixD3& xa,
	      const Subset& mstag)
  {
    START_CODE();

    LatticeBoolean bad;
    int numbad;
    
    reunit_t<LatticeColorMatrixD3, LatticeComplexD, LatticeRealD, Subset>(xa, bad, numbad, REUNITARIZE, mstag);
    
    END_CODE();
  }

  // Overloaded definitions, with numbad and ruflag
  // Single

  void reunit(LatticeColorMatrixF3& xa,
	      int& numbad, 
	      enum Reunitarize ruflag)
  {
    START_CODE();

    LatticeBoolean bad;
    
    reunit_t<LatticeColorMatrixF3, LatticeComplexF, LatticeRealF, Subset>(xa, bad, numbad, REUNITARIZE, all);
    
    END_CODE();
  }

  // DOUBLE
  void reunit(LatticeColorMatrixD3& xa,
	      int& numbad, 
	      enum Reunitarize ruflag)
  {
    START_CODE();

    LatticeBoolean bad;
    
    reunit_t<LatticeColorMatrixD3, LatticeComplexD, LatticeRealD, Subset>(xa, bad, numbad, REUNITARIZE, all);
    
    END_CODE();
  }
  

  // SINGLE
  void reunit(LatticeColorMatrixF3& xa,
	      int& numbad, 
	      enum Reunitarize ruflag,
	      const Subset& mstag)
  {
    START_CODE();

    LatticeBoolean bad;
    
    reunit_t<LatticeColorMatrixF3, LatticeComplexF, LatticeRealF, Subset>(xa, bad, numbad, REUNITARIZE, mstag);
    
    END_CODE();
  }
  
   // DOUBLE
   void reunit(LatticeColorMatrixD3& xa,
	      int& numbad, 
	      enum Reunitarize ruflag,
	      const Subset& mstag)
  {
    START_CODE();

    LatticeBoolean bad;
    
    reunit_t<LatticeColorMatrixD3, LatticeComplexD, LatticeRealD, Subset>(xa, bad, numbad, REUNITARIZE, mstag);
    
    END_CODE();
  }
  // Overloaded definitions, with bad, numbad and ruflag
  // SINGLE
  void reunit(LatticeColorMatrixF3& xa, 
	      LatticeBoolean& bad, 
	      int& numbad, 
	      enum Reunitarize ruflag)
  {
    reunit_t<LatticeColorMatrixF3, LatticeComplexF, LatticeRealF, Subset>(xa, bad, numbad, ruflag, all);
  }

  // DOUBLE
  void reunit(LatticeColorMatrixD3& xa, 
	      LatticeBoolean& bad, 
	      int& numbad, 
	      enum Reunitarize ruflag)
  {
    reunit_t<LatticeColorMatrixD3, LatticeComplexD, LatticeRealD, Subset>(xa, bad, numbad, ruflag, all);
  }

  // Single
  void reunit(LatticeColorMatrixF3& xa, 
	      LatticeBoolean& bad, 
	      int& numbad, 
	      enum Reunitarize ruflag,
	      const Subset& mstag)
  {
    reunit_t<LatticeColorMatrixF3, LatticeComplexF, LatticeRealF, Subset>(xa, bad, numbad, ruflag, mstag);
  }
  
  
 // Double
  void reunit(LatticeColorMatrixD3& xa, 
	      LatticeBoolean& bad, 
	      int& numbad, 
	      enum Reunitarize ruflag,
	      const Subset& mstag)
  {
    reunit_t<LatticeColorMatrixD3, LatticeComplexD, LatticeRealD, Subset>(xa, bad, numbad, ruflag, mstag);
  }

} // End namespace
