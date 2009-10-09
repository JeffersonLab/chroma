// $Id: polar_dec.cc,v 3.3 2009-10-09 15:33:43 bjoo Exp $
/*! \file
 *  \brief Decompose a complex matrix as C = exp(i\alpha) V P
 */

#include "chromabase.h"
#include "meas/gfix/polar_dec.h"
#include "util/gauge/reunit.h"

namespace Chroma { 
//! Decompose a complex matrix as C = exp(i\alpha) V P
/*!
 * \ingroup gfix
 *
 * Decompose a complex matrix as C = exp(i\alpha) V P
 * with V SU(Nc) and P = (C^\dagger C)^{1/2} positive hermitian
 *
 * \param c        complex Nc x Nc matrix ( Modify )
 *                 on exit it contains the hermitian matrix P
 * \param v        the projected SU(Nc) Matrix ( Write )
 * \param alpha    the phase ( Write )
 * \param JacAccu  accuracy in the Jacobi iteration ( Read )
 * \param JacMax   maximum number of Jacobi iterations ( Read ) 
 */

void polar_dec(LatticeColorMatrix& c, LatticeColorMatrix& v,
	       LatticeReal& alpha, const Real& JacAccu, int JacMax)
{
  LatticeColorMatrix u_tmp;
  LatticeColorMatrix w_tmp;
  LatticeColorMatrix v_tmp;

  multi2d<LatticeComplex> mat1(Nc, Nc);
  multi2d<LatticeComplex> mat2(Nc, Nc);

  LatticeComplex lc1;
  LatticeComplex lc2;
  LatticeComplex v12;
  LatticeComplex v21;
  multi1d<LatticeReal> diag(Nc);
  LatticeReal rescale;
  LatticeReal hundred;
  LatticeReal v11;
  LatticeReal dd;
  LatticeReal diff_diag;
  LatticeReal tt;
  LatticeReal cc;
  LatticeReal ss;
  LatticeReal al1;
  LatticeReal al2;
  LatticeReal lr1;
  LatticeReal lr2;
  LatticeBoolean l_rot;
  LatticeBoolean lb1;
  LatticeBoolean lb2;
  LatticeBoolean lb3;

  Double off_d;
  Double diff_sq;
  Double det_diff;
  Double unit_test;
  int iter;
  int i_rot;
  int n_rot;
  int numbad;

  START_CODE();
  
  Real acc_sq = JacAccu * JacAccu / Real(100);
              
  u_tmp = adj(c) * c;

  /* Bring diagonal elements to be O(1) */
  rescale = real(trace(u_tmp));
  lr1 = 1;
  lr1 = lr1 / rescale;
  u_tmp = u_tmp * lr1;

  // Extract initial components 
  for(int i=0; i < Nc; ++i)
    for(int j=0; j < Nc; ++j)
      mat1[i][j] = peekColor(u_tmp, i, j);

  for(int i=0; i<Nc; i++)
  {
    for(int j=0; j<i; j++)
    {
      mat2[i][j] = 0;
      mat2[j][i] = 0;
    }
    mat2[i][i] = 1;
    diag[i] = real(mat1[i][i]);
  }

  /* Diagonalize mat1 by Jacobi iteration and keep "rotation" in mat2 */
                                    
  hundred = 100;
  iter = 0;
  i_rot = Layout::vol();

  while ( iter < JacMax && i_rot > 0 )
  {
    iter++;
    i_rot = 0;

    for(int j = 1; j < Nc; j++)
      for(int i = 0; i < j; i++)
      {
	dd = real(mat1[i][j] * mat1[j][i]);
	lr1 = fabs(diag[i] * diag[j] * acc_sq);

	l_rot = dd > lr1;
	n_rot = toInt(sum(where(l_rot, LatticeInteger(1), LatticeInteger(0))));

	if( n_rot > 0 )
	{
	  i_rot += n_rot;
	  dd = sqrt(dd);
	  lr1 = dd * hundred;
	  diff_diag = diag[j] - diag[i];

	  lr2 = fabs(diff_diag);
	  lr1 += lr2;
	  lb1 = lr2 == lr1;
	  lb2 = l_rot & lb1;
	  lr1 = dd / diff_diag;
	  copymask(tt, lb2, lr1);

	  lb1 = !lb1;
	  lb2 = l_rot & lb1;
	  lr1 = diff_diag / dd;
	  LatticeReal theta = 0.5 * lr1;
	  lr1 = fabs(theta) + sqrt(1 + theta * theta);
	  lr2 = 1 / lr1;
	  copymask(tt, lb2, lr2);

	  lr1 = 0;
	  lb1 = theta < lr1;
	  lb3 = lb2 & lb1;
	  copymask(tt, lb3, LatticeReal(-tt));

	  lr1 = 1 / sqrt(1 + tt*tt);
	  lr2 = tt * lr1;
	  copymask(cc, l_rot, lr1);
	  copymask(ss, l_rot, lr2);

	  lr1 = cc * cc;
	  al1 = lr1 * diag[i];
	  al2 = lr1 * diag[j];
	  lr1 = ss * ss;
	  al1 += lr1 * diag[j];
	  al2 += lr1 * diag[i];
	  lr1 = cc * ss * dd;

	  lr1 *= 2;
	  al1 -= lr1;
	  al2 += lr1;
	  copymask(diag[i], l_rot, al1);
	  copymask(diag[j], l_rot, al2);

	  v11 = cc;
	  lr1 = ss / dd;
	  v12 = mat1[j][i] * lr1;
	  v21 = -adj(v12);
	  lc1 = 0;
	  copymask(mat1[j][i], l_rot, lc1);
	  copymask(mat1[i][j], l_rot, lc1);

	  /* Rotate the remaining matrix elements */
	  for(int k = 0; k < Nc; k++)
	    if( k != i && k != j )
	    {
	      lc1 = mat1[k][i] * v11 - mat1[k][j] * v12;
	      lc2 = mat1[k][j] * v11 - mat1[k][i] * v21;
	      copymask(mat1[k][i], l_rot, lc1);
	      copymask(mat1[k][j], l_rot, lc2);
	      copymask(mat1[i][k], l_rot, LatticeComplex(adj(lc1)));
	      copymask(mat1[j][k], l_rot, LatticeComplex(adj(lc2)));
	    }

	  /* Rotate the eigenvector matrix mat2 */
	  for(int k = 0; k < Nc; k++)
	  {
	    lc1 = mat2[i][k] * v11 + mat2[j][k] * v21;
	    lc2 = mat2[j][k] * v11 + mat2[i][k] * v12;
	    copymask(mat2[i][k], l_rot, lc1);
	    copymask(mat2[j][k], l_rot, lc2);
	  }
	}
      }
  }

                                    
  for(int i=0; i < Nc; ++i)
    for(int j=0; j < Nc; ++j)
      pokeColor(w_tmp, mat2[i][j], i, j);

  v_tmp = w_tmp;

  /* Check unitarity */
  reunit(v_tmp, lb1, numbad, REUNITARIZE_LABEL);
  if (numbad > 0)
    printf("Warning: %d rotation matrices not unitary\n", numbad);
/*  */

  /* Check diagonalization */
  lr1 = 0;
  off_d = 0;
  numbad = 0;
  for(int i=0; i<Nc; i++)
  {
    for(int j=0; j<i; j++)
    {
      off_d += norm2(mat1[j][i]);
      mat1[j][i] = 0;
      mat1[i][j] = 0;
    }
    numbad += toInt(sum(where(diag[i] <= lr1, LatticeInteger(1), LatticeInteger(0))));
    mat1[i][i] = cmplx(diag[i],lr1);
  }

  for(int i=0; i < Nc; ++i)
    for(int j=0; j < Nc; ++j)
      pokeColor(v, mat1[i][j], i, j);

  v_tmp = v * w_tmp ;
  v = adj(w_tmp) * v_tmp - u_tmp;
  diff_sq = norm2(v) / Double(Layout::vol());
  off_d /= Double(Layout::vol());
  if (numbad > 0)
  {
    QDPIO::cout << "Warning: " << numbad << " C matrices have zero determiant" << endl;
#if 0
    push(nml_out,"Bad_C_matrices");
    write(nml_out, "numbad", numbad);
    pop(nml_out);
#endif
  }
  
  /* Rescale eigenvalues and construct P^{-1} in v and c*P^{-1} in u_tmp */
  for(int i=0; i<Nc; i++)
  {
    diag[i] = sqrt(diag[i] * rescale);
    lr2 = 1 / diag[i];
    mat1[i][i] = cmplx(lr2,lr1);
  }

  for(int i=0; i < Nc; ++i)
    for(int j=0; j < Nc; ++j)
      pokeColor(v, mat1[i][j], i, j);

  u_tmp = v * w_tmp;
  v = adj(w_tmp) * u_tmp;
  u_tmp = c * v;

  /* Constuct P in c */
  for(int i=0; i<Nc; i++)
  {
    mat1[i][i] = cmplx(diag[i],lr1);
  }

  for(int i=0; i < Nc; ++i)
    for(int j=0; j < Nc; ++j)
      pokeColor(c, mat1[i][j], i, j);

  v_tmp = c * w_tmp;
  c = adj(w_tmp) * v_tmp;

  /* Check that that u_tmp^dagger * u_tmp = 1 */
  v_tmp = 1;
  v_tmp -= adj(u_tmp) * u_tmp;
  unit_test = norm2(v_tmp);
  unit_test /= Double(Layout::vol());

    
  /* Now we just need to remove the phase alpha from u_tmp to make it v */
  /* We compute the determinant by LU decomposition */
  for(int i=0; i < Nc; ++i)
    for(int j=0; j < Nc; ++j)
      mat1[i][j] = peekColor(u_tmp, i, j);

  lr2 = 1;
  for(int j = 0; j < Nc; j++)
  {
    for(int i = 0; i <= j; i++)
    {
      lc1 = mat1[j][i];
      for(int k = 0; k < i; k++)
	lc1 -= mat1[k][i] * mat1[j][k];

      mat1[j][i] = lc1;
    }

    for(int i = (j+1); i < Nc; i++)
    {
      lc1 = mat1[j][i];
      for(int k = 0; k < j; k++)
	lc1 -= mat1[k][i] * mat1[j][k];

      lc2 = adj(mat1[j][j]) * mat1[j][j];
      lr1 = real(lc2);
      lr1 = lr2 / lr1;
      lc2 = lc1 * lr1;
      lc1 = adj(mat1[j][j]) * lc2;
      mat1[j][i] = lc1;
    }
  }

  /* The determinant */
  lc1 = mat1[0][0] * mat1[1][1];
  for(int k = 2; k < Nc; k++)
  {
    lc2 = mat1[k][k] * lc1;
    lc1 = lc2;
  }
  lr1 = real(lc1);
  lr2 = imag(lc1);
  alpha = atan2(lr2, lr1);

  lr1 = 1;
  lr2 = real(trace(adj(lc1) * lc1));
  lr1 -= lr2;
  det_diff = norm2(lr1) / Double(Layout::vol());

  /* Make u_tmp unitary in v */
  lr1 = Real(Nc);
  alpha = alpha / lr1;
  lr2 = log(lr2);
  lr2 = lr2 / lr1;
  lr2 = -lr2;
  lr2 = exp(lr2);
  cc = cos(alpha);
  ss = sin(alpha);
  lr1 = cc * lr2;
  lr2 = -(ss * lr2);
  lc1 = cmplx(lr1,lr2);
  v = u_tmp * lc1;
        
  /* Make v unitary */
  reunit(v, lb1, numbad, REUNITARIZE_LABEL);

#if 0
  XMLBufferWriter xml_out;
  push(xml_out,"Diagonalization_test");
  write(xml_out, "iter",iter);
  write(xml_out, "off_d",off_d);
  write(xml_out, "diff_sq", diff_sq);
  write(xml_out, "unit_test", unit_test);
  write(xml_out, "det_diff", det_diff);
  write(xml_out, "numbad", numbad);
  pop(xml_out);
#endif

#if 0
  // DEBUG DEBUG DEBUG
  QDPIO::cout << "polar_dec::iter " <<  iter << "\n" ;
  QDPIO::cout << "polar_dec::off_d " <<  off_d << "\n" ;
  QDPIO::cout << "polar_dec::diff_sq " <<  diff_sq << "\n" ;
  QDPIO::cout << "polar_dec::unit_test " <<  unit_test << "\n" ;
  QDPIO::cout << "polar_dec::det_diff " <<  det_diff<< "\n" ;
  QDPIO::cout << "polar_dec::numbad " <<   numbad<< "\n" ;
#endif  

  END_CODE();
}

}; // End namespace Chroma
