// $Id: sn_jacob.cc,v 1.1 2004-01-04 21:56:04 edwards Exp $
/*! \file
 *  \brief Single-node Jacobi routine
 */

#error "CONVERSION NOT COMPLETE"

#include "chromabase.h"
#include "meas/eig/ritz.h"

using namespace QDP;

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
 *  N_Count	Number of Jacobi iters		(Write) 
 */

void SN_Jacob(psi, N_eig, lambda, off_diag, Toler, N_max,
	      Ncb, n_count)

multi2d<LatticeFermion> psi(Ncb, N_eig);
int N_eig;
multi1d<Real> lambda(N_eig);
multi1d<Complex> off_diag(N_eig*(N_eig-1)/2);
Real Toler;
int N_max;
int Ncb;
int n_count;
{
  START_CODE("SN_Jacob");
  
  LatticeFermion psi_t1;
  LatticeFermion psi_t2;
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
  
  Real tol_sq = Toler * Toler;
            
  for(int k = 0; k <= N_max; k++)
  {
    i_rot = 0;
    ij = 0;
    for(j = 1; j < N_eig; j++)
      for(i = 0; i < j; i++)
      {
	dd = real(adj(off_diag[ij]) * off_diag[ij]);
	ftmp = fabs(tol_sq * lambda[i] * lambda[j]);

	if( dd > ftmp)
	{
	  /* Make a rotation to set off-diagonal part to zero */
	  i_rot++;
	  dd = sqrt(dd);
	  acc = Real(100) * dd;
	  diff_l = lambda[j] -  lambda[i];
	  ftmp = fabs(diff_l);

	  if( (ftmp+acc) == ftmp )
	  {
	    t = dd / diff_l;
	  }
	  else
	  {
	    theta = 0.5 * diff_l / dd;
	    t = sqrt(1.0 + theta*theta);
	    ftmp = fabs(theta);
	    t = 1.0 / (ftmp+t);
	    if( theta < 0.0 ) 
	      t = -t;
	  }

	  if( diff_l >= 0 )
	  {
	    c = 1.0 / sqrt(1.0 + t*t);
	    s = - t * c;
	  }
	  else
	  {
	    s = 1.0 / sqrt(1.0 + t*t);
	    c = t * s;
	  }

	  ftmp = c * c;
	  al1 = ftmp * lambda[i];
	  al2 = ftmp * lambda[j];
	  ftmp = s * s;
	  al1 += ftmp * lambda[j];
	  al2 += ftmp * lambda[i];
	  ftmp = 2.0 * dd * s * c;
	  al1 += ftmp;
	  al2 -= ftmp;
	  lambda[i] = al1;
	  lambda[j] = al2;
	  v11 = c;
	  v12 = (s / dd) * off_diag[ij];
	  v21 = -adj(v12);
	  off_diag[ij] = 0;

	  /* Now rotate the eigenvectors */
	  for(cb = 0; cb < Ncb; cb++)
	  {
	    psi_t1 = psi[i][cb] * v11;
	    /* psi_t1 += psi[j][cb] * adj(v12); Wrong?? */
	    psi_t1 -= psi[j][cb] * v21;
	    psi_t2 = psi[j][cb] * v11;
	    psi_t2 -= psi[i][cb] * v12;
	    psi[i][cb] = psi_t1;
	    psi[j][cb] = psi_t2;
	  }

	  /* Rotate the other matrix elements */
	  for(m = 0; m < N_eig; m++)
	    if( m != i && m != j )
	    {
	      if( m < i )
	      {
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
	      else if( m < j)
	      {
		mi = m * (m-1) / 2 + i;
		mj = j * (j-1) / 2 + m;
		ctmp1 = adj(off_diag[mi]) * v11;
		ctmp1 -= off_diag[mj] * v21;
		ctmp2 = off_diag[mj] * v11;
		ctmp2 -= adj(off_diag[mi]) * v12;
		off_diag[mi] = adj(ctmp1);
		off_diag[mj] = ctmp2;
	      }
	      else
	      {
		mi = m * (m-1) / 2 + i;
		mj = m * (m-1) / 2 + j;
		ctmp1 = adj(off_diag[mi]) * v11;
		ctmp1 -= adj(off_diag[mj]) * v21;
		ctmp2 = adj(off_diag[mj]) * v11;
		ctmp2 -= adj(off_diag[mi]) * v12;
		off_diag[mi] = adj(ctmp1);
		off_diag[mj] = adj(ctmp2);
	      }
	    }
	}

	ij++;
      }

    if( i_rot == 0 )
    {
      n_count = k;
      QDPIO::cout << "Jacobi converged after " << k << " iters" << endl;

      /* Sort the eigenvalues */
      for(j = 1; j < N_eig; j++)
	for(i = 0; i < j; i++)
	{
	  ftmp = fabs(lambda[j]);
	  dd = fabs(lambda[i]);
	  /* if( lambda[j] < lambda[i] ) */
	  if( ftmp < dd )
	  {
	    ftmp = lambda[i];
	    lambda[i] = lambda[j];
	    lambda[j] = ftmp;
	    for(cb = 0; cb < Ncb; cb++)
	    {
	      psi_t1 = psi[i][cb];
	      psi_t2 = psi[j][cb];
	      psi[j][cb] = psi_t1;
	      psi[i][cb] = psi_t2;
	    }
	  }
	}

      END_CODE("SN_Jacob");
      return;
    }
  }

  n_count = k;
  QDPIO::cerr << "too many Jacobi iterations " << k << endl;
  QDP_abort(1);
  END_CODE("SN_Jacob");
}
