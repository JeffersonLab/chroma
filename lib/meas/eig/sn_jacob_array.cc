// $Id: sn_jacob_array.cc,v 3.1 2007-02-22 21:11:49 bjoo Exp $
/*! \file
 *  \brief Single-node Jacobi routine
 */

#include "chromabase.h"
#include "meas/eig/sn_jacob_array.h"

namespace Chroma {

//! Single-node Jacobi rotation
/*!
 * \ingroup eig
 *
 * This subroutine contains a "single node" Jacobi routine
 * to be used with the Ritz functional eigenvialue/vector finder.
 *
 *
 *  Psi		Eigenvectors			(Modify)
 *  N_eig	Eigenvalue number 		(Read)
 *  lambda	Diagonals / Eigenvalues		(Modify)
 *  off_diag	Upper triang off-diag matrix elems	(Modify)
 *  Toler	Tolerance for off-diag elems	(Read)
 *  N_max	Maximal number of Jacobi iters	(Read)
 *  Ncb		Number of sublattices		(Read)
 *  \param sub         Subset to use                   (Read) 
 *
 * \return 	Number of Jacobi iters		(Write) 
 */
template <typename T>
int SN_Jacob_Array_t(multi2d<T>& psi, 
		     const int N_eig, 
		     multi1d<Real>& lambda, 
		     multi1d<Complex>& off_diag, 
		     Real tolerance, 
		     int N_max,
		     const Subset& sub)
{
  START_CODE();
  
  int N5 = psi.size1();

  int n_count;
  multi1d<T> psi_t1(N5);
  multi1d<T> psi_t2(N5);
  Complex ctmp1;
  Complex ctmp2;
  Complex v12;
  Complex v21;
  Real v11;
  Real dd;
  Real ftmp;
  Real diff_l;
  Real theta;
  Real t;
  Real acc;
  Real c;
  Real s;
  Real al1;
  Real al2;
  int cb;
  int i;
  int j;
  int ij;
  int m;
  int mi;
  int mj;
  int i_rot;
  int k;

  Real tol_sq = tolerance * tolerance;
            
  for(k = 1; k <= N_max; k++) 
  {
    i_rot = 0;
    ij = 0;
    
    for(j = 1; j < N_eig; j++) 
    {
      for(i = 0; i < j; i++) 
      {
	dd = real(conj(off_diag[ij]) * off_diag[ij]);
	ftmp = fabs(tol_sq * lambda[i] * lambda[j]);

	if( toBool( dd > ftmp) ) 
	{
	  // Make a rotation to set off-diagonal part to zero 
	  i_rot++;
	  dd = sqrt(dd);

	  acc = Real(100) * dd;

	  diff_l = lambda[j] -  lambda[i];
	  ftmp = fabs(diff_l);

	  if( toBool( (ftmp+acc) == ftmp ) ) 
	  {
	    t = dd / diff_l;
	  }
	  else 
	  {
	    theta = Real(0.5) * diff_l / dd;
	    t = sqrt(Real(1) + theta*theta);
	    ftmp = fabs(theta);
	    t = Real(1) / (ftmp+t);
	    if( toBool( theta < Real(0) ) ) { 
	      t = -t;
	    }
	  }

	  if( toBool( diff_l >= Real(0) ) ) {
	    c = Real(1) / sqrt( Real(1) + t*t);
	    s = - t * c;
	  }
	  else
	  {
	    s = Real(1) / sqrt( Real(1) + t*t);
	    c = t * s;
	  }

	  ftmp = c * c;
	  al1 = ftmp * lambda[i];
	  al2 = ftmp * lambda[j];
	  ftmp = s * s;
	  al1 += ftmp * lambda[j];
	  al2 += ftmp * lambda[i];
	  ftmp = Real(2) * dd * s * c;
	  al1 += ftmp;
	  al2 -= ftmp;
	  lambda[i] = al1;
	  lambda[j] = al2;
	  v11 = c;
	  v12 = (s / dd) * off_diag[ij];
	  v21 = -conj(v12);
	  off_diag[ij] = Real(0);

	  // Now rotate the eigenvectors */
	  for(int n=0; n < N5; n++) {
	    psi_t1[n][sub] = psi[i][n] * v11;

	    /* psi_t1 += psi[j][cb] * adj(v12); Wrong?? */
	    psi_t1[n][sub] -= psi[j][n] * v21;
	    psi_t2[n][sub] = psi[j][n] * v11;
	    psi_t2[n][sub] -= psi[i][n] * v12;
	    psi[i][n][sub] = psi_t1[n];
	    psi[j][n][sub] = psi_t2[n];
	  }

	  // Rotate the other matrix elements 
	  for(m = 0; m < N_eig; m++) {
	    if( m != i && m != j ) {
	    
	      if( m < i ) {
	     
		mi = i * (i-1) / 2 + m;
		mj = j * (j-1) / 2 + m;
		ctmp1 = off_diag[mi] * v11;
		/* ctmp1 += off_diag[mj] * adj(v12); ok */
		ctmp1 -= off_diag[mj] * v21;
		ctmp2 = off_diag[mj] * v11;
		ctmp2 -= off_diag[mi] * v12;
		off_diag[mi] = ctmp1;
		off_diag[mj] = ctmp2;
	      }
	      else if( m < j) {
	      
		mi = m * (m-1) / 2 + i;
		mj = j * (j-1) / 2 + m;
		ctmp1 = conj(off_diag[mi]) * v11;
		ctmp1 -= off_diag[mj] * v21;
		ctmp2 = off_diag[mj] * v11;
		ctmp2 -= conj(off_diag[mi]) * v12;
		off_diag[mi] = conj(ctmp1);
		off_diag[mj] = ctmp2;
	      }
	      else {
	      
		mi = m * (m-1) / 2 + i;
		mj = m * (m-1) / 2 + j;
		ctmp1 = conj(off_diag[mi]) * v11;
		ctmp1 -= conj(off_diag[mj]) * v21;
		ctmp2 = conj(off_diag[mj]) * v11;
		ctmp2 -= conj(off_diag[mi]) * v12;
		off_diag[mi] = conj(ctmp1);
		off_diag[mj] = conj(ctmp2);
	      }
	    }
	  }
	}
	
	ij++;
      }
    }

    if( i_rot == 0 ) 
    {
      n_count = k;
      QDPIO::cout << "Jacobi converged after " << k << " iters" << endl;
      
      // Sort the eigenvalues
      // In order of increasing modulus
      for(j = 1; j < N_eig; j++) 
      {
	for(i = 0; i < j; i++) 
	{
	  ftmp = fabs(lambda[j]);
	  dd = fabs(lambda[i]);

	  // if( | lambda[j] | < | lambda[i] | )
	  if( toBool( ftmp < dd ) ) 
	  {
	    ftmp = lambda[i];
	    lambda[i] = lambda[j];
	    lambda[j] = ftmp;
	    
	    for(int n=0; n < N5; n++) 
	    {
	      psi_t1[n][sub] = psi[i][n];
	      psi[i][n][sub] = psi[j][n];
	      psi[j][n][sub] = psi_t1[n];
	    }
	  }
	}
      }
      END_CODE();
      return n_count;
    }
  }

  n_count = k;
  QDP_error_exit("too many Jacobi iterations: %d\n" ,k);
  END_CODE();

  return n_count;
}


int SN_Jacob_Array(multi2d<LatticeFermion>& psi, 
		   const int N_eig, 
		   multi1d<Real>& lambda, 
		   multi1d<Complex>& off_diag, 
		   Real tolerance, 
		   int N_max,
		   const Subset& sub)
{
  return SN_Jacob_Array_t(psi, N_eig, lambda, off_diag, tolerance, N_max, sub);
}

}  // end namespace Chroma
