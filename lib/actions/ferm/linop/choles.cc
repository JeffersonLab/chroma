// $Id: choles.cc,v 1.1 2005-02-17 02:52:36 edwards Exp $

#error "NOT FULLY CONVERTED"

/*# CHOLES   - Cholesky decompose a matrix and use it to */
/*#            compute  lower(A^-1) = lower((L.L^dag)^-1) */
/*#            Adapted from Golub and Van Loan, Matrix Computations, 2nd, Sec 4.2.4 */

/*# Arguments: */

/*#  Mat         Full hermitian matrix              (Modify) */

include(types.mh)

SUBROUTINE(choles, mat, n)

multi2d<Complex> mat(n, n);
int n;

{ /* Local variables */
  include(COMMON_DECLARATIONS)

  multi1d<Real> diag_g(n);
  multi1d<Complex> v1(n);
  Complex sum;
  Real one;
  Real zero;
  int i;
  int j;
  int k;
  
  START_CODE();
  
        one = 1;
  zero = 0;
  
  /*# Cholesky decompose  A = L.L^dag */
  
  /*# NOTE!!: I can store this matrix in  mat, but will need a */
  /*#   temporary  diag */
  for(j = 0; j < n; ++j)
  {
    /*# Multiply mater mass term against basis vector.  */
    /*# Actually, I need a column of the lower triang matrix mat. */
    for(i = j; i < n; ++i)
    {
      v1[i] = adj(mat[i][j]);
    }

    /*# Back to cholesky */
    /*# forward substitute */
    for(k = 0; k < j; ++k)
    {
      for(i = j; i < n; ++i)
      {
	v1[i] -= adj(mat[k][j]) * mat[k][i];
      }
    }

    /*# The diagonal is (should be!!) real and positive */
    diag_g[j] = real(v1[j]);
    diag_g[j] = sqrt(diag_g[j]);
    diag_g[j] = one / diag_g[j];

    /*# backward substitute */
    for(i = j+1; i < n; ++i)
    {
      mat[j][i] = v1[i] * diag_g[j];
    }
  }


  /*# Use forward and back substitution to construct  mat = lower(A^-1) */
  for(k = 0; k < n; ++k)
  {
    for(i = 0; i < k; ++i)
      v1[i] = 0;

    /*# Forward substitution */
    v1[k] = cmplx(diag_g[k],zero);

    for(i = k+1; i < n; ++i)
    {
      sum = 0;

      for(j = k; j < i; ++j)
      {
	sum -= mat[j][i] * v1[j];
      }

      v1[i] = sum * diag_g[i];
    }

    /*# Backward substitution */
    v1[n-1] = v1[n-1] * diag_g[n-1];

    for(i = n-2; i >= k; --i)
    {
      sum = v1[i];

      for(j = i+1; j < n; ++j)
      {
	sum -= adj(mat[i][j]) * v1[j];
      }

      v1[i] = sum * diag_g[i];
    }

    /*# Overwrite column k of mat */
    for(i = k; i < n; ++i)
    {
      mat[k][i] = v1[i];
    }
  }
  
  /* Make it hermitian */
  for(j = 0; j < n; ++j)
    for(i = j+1; i < n; ++i)
    {
      mat[i][j] = adj(mat[j][i]);
    }

        
  END_CODE();
}
