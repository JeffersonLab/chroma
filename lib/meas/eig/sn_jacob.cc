// $Id: sn_jacob.cc,v 3.1 2007-02-22 21:11:49 bjoo Exp $
/*! \file
 *  \brief Single-node Jacobi routine
 */

#include "chromabase.h"
#include "meas/eig/sn_jacob.h"

namespace Chroma {

//! Single-node Jacobi rotation
/*!
 * \ingroup eig
 *
 * This subroutine contains a "single node" Jacobi routine
 * to be used with the Ritz functional eigenvialue/vector finder.
 *
 *  \param psi		Eigenvectors			(Modify)
 *  \param N_eig	Eigenvalue number 		(Read)
 *  \param lambda	Diagonals / Eigenvalues		(Modify)
 *  \param off_diag	Upper triang off-diag matrix elems	(Modify)
 *  \param tolererance	Tolerance for off-diag elems	(Read)
 *  \param N_max	Maximal number of Jacobi iters	(Read)
 *  \param sub          Subset to use                   (Read) 
 *
 *  \return 	        Number of Jacobi iters		(Write) 
 */
template <typename T>
int SN_Jacob_t(multi1d<T>& psi, 
	       const int N_eig, 
	       multi1d<Real>& lambda, 
	       multi1d<Complex>& off_diag, 
	       Real tolerance, 
	       int N_max,
	       const Subset& sub) 
{
  START_CODE();
  
  int n_count;
  T psi_t1;
  T psi_t2;
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


	  if( toBool( (ftmp+acc) == ftmp ) ) {
	    
	    t = dd / diff_l;
	  }
	  else {
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
	  psi_t1[sub] = psi[i] * v11;
	  /* psi_t1 += psi[j][cb] * adj(v12); Wrong?? */
	  psi_t1[sub] -= psi[j] * v21;
	  psi_t2[sub] = psi[j] * v11;
	  psi_t2[sub] -= psi[i] * v12;
	  psi[i][sub] = psi_t1;
	  psi[j][sub] = psi_t2;
	  

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

    if( i_rot == 0 ) {
    
      n_count = k;
      QDPIO::cout << "Jacobi converged after " << k << " iters" << endl;
      
      // Sort the eigenvalues
      // In order of increasing modulus
      for(j = 1; j < N_eig; j++) {
	for(i = 0; i < j; i++) {
	  
	  ftmp = fabs(lambda[j]);
	  dd = fabs(lambda[i]);

	  // if( | lambda[j] | < | lambda[i] | )
	  if( toBool( ftmp < dd ) ) {
	  
	    ftmp = lambda[i];
	    lambda[i] = lambda[j];
	    lambda[j] = ftmp;
	    
	    
	    psi_t1[sub] = psi[i];
	    psi[i][sub] = psi[j];
	    psi[j][sub] = psi_t1;
	  }
	}
      }
      END_CODE();
      return n_count;
    }
  }

  n_count = k;
  QDP_error_exit("too many Jacobi iterations: %d\n" , k);
  END_CODE();
  return n_count;
}


int SN_Jacob(multi1d<LatticeFermion>& psi, 
	     const int N_eig, 
	     multi1d<Real>& lambda, 
	     multi1d<Complex>& off_diag, 
	     Real tolerance, 
	     int N_max,
	     const Subset& sub)
{
  return SN_Jacob_t(psi, N_eig, lambda, off_diag, tolerance, N_max, sub);
}

}  // end namespace Chroma
