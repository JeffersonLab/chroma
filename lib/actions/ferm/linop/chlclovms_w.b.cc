// $Id: chlclovms_w.b.cc,v 2.0 2005-09-25 21:04:28 edwards Exp $

#error "NOT FULLY CONVERTED"

/*# This routine is specific to Wilson fermions! */

/*# CHLCLOVMS - Cholesky decompose the clover mass term and uses it to */
/*#             compute  lower(A^-1) = lower((L.L^dag)^-1) */
/*#             Adapted from Golub and Van Loan, Matrix Computations, 2nd, Sec 4.2.4 */

/*# construct   */

/*# Arguments: */

/*#  Clov         Lower tri. mat. (no diag)           (Read) */
/*#  Diag_clov    diag(clov)   	                    (Read) */
/*#  Invclov      lower(clov^-1) (no diag)            (Write) */
/*#  Diag_invclov diag(invclov)                       (Write) */
/*#  DetP         flag whether to compute determinant (Read) */
/*#  logdet       logarithm of the determinant        (Write) */

include(types.mh)

SUBROUTINE(chlclovms, clov, diag_clov, invclov, diag_invclov, DetP, logdet)


LATTICE_TRIANGULAR(clov);
LATTICE_DIAG_TRIANGULAR(diag_clov);
LATTICE_TRIANGULAR(invclov);
LATTICE_DIAG_TRIANGULAR(diag_invclov);
int DetP;
Double logdet;

{ /* Local variables */
  include(COMMON_DECLARATIONS)

  multi1d<LatticeReal> diag_g(`(2*Nc)');
  multi1d<LatticeComplex> v1(`(2*Nc)');
  LatticeComplex sum;
  LatticeReal one;
  LatticeReal zero;
  LatticeReal log_diag;
  LatticeReal lrtmp;
  int i;
  int j;
  int k;
  int elem_ij;
  int elem_ji;
  int elem_ik;
  int elem_jk;
  int n;
  int s;
  
  START_CODE();
  
  n = 2*Nc;
  
  if ( n < 3 )
    QDP_error_exit("Matrix is too small", Nc, Ns);
  
            one = 1;
  zero = 0;
  
  if ( DetP == YES )
  {
            log_diag = 0;
  }
  
  /*# Cholesky decompose  A = L.L^dag */
  
  /*# NOTE!!: I can store this matrix in  invclov, but will need a */
  /*#   temporary  diag */
  for(s = 0;s  <= ( 1); ++s )
  {
    for(j = 0;j  < ( n); ++j )
    {
      /*# Multiply clover mass term against basis vector.  */
      /*# Actually, I need a column of the lower triang matrix clov. */
      v1[j] = cmplx(diag_clov[s][j],zero);
      for(i = j+1;i  < ( n); ++i )
      {
	elem_ij = i*(i-1)/2 + j;
	v1[i] = clov[s][elem_ij];
      }

      /*# Back to cholesky */
      /*# forward substitute */
      for(k = 0;k  < ( j); ++k )
      {
	for(i = j;i  < ( n); ++i )
	{
	  elem_jk = j*(j-1)/2 + k;
	  elem_ik = i*(i-1)/2 + k;
	  v1[i] -= adj(invclov[s][elem_jk]) * invclov[s][elem_ik];
	}
      }

      /*# The diagonal is (should be!!) real and positive */
      diag_g[j] = real(v1[j]);

      /*#+ */
      /*# Squeeze in computation of log(Det), if desired */
      /*#- */
      if ( DetP == YES )
      {
	lrtmp = log(diag_g[j]);
	log_diag += lrtmp;
      }

      diag_g[j] = sqrt(diag_g[j]);
      diag_g[j] = one / diag_g[j];

      /*# backward substitute */
      for(i = j+1;i  < ( n); ++i )
      {
	elem_ij = i*(i-1)/2 + j;
	invclov[s][elem_ij] = v1[i] * diag_g[j];
      }
    }


    /*# Use forward and back substitution to construct  invclov = lower(A^-1) */
    for(k = 0;k  < ( n); ++k )
    {
      for(i = 0;i  < ( k); ++i )
	v1[i] = 0;

      /*# Forward substitution */
      v1[k] = cmplx(diag_g[k],zero);

      for(i = k+1;i  < ( n); ++i )
      {
	sum = 0;

	for(j = k;j  < ( i); ++j )
	{
	  elem_ij = i*(i-1)/2 + j;
	  sum -= invclov[s][elem_ij] * v1[j];
	}

	v1[i] = sum * diag_g[i];
      }

      /*# Backward substitution */
      v1[n-1] = v1[n-1] * diag_g[n-1];

      for(i = n-2;i  >= ( k); --i )
      {
	sum = v1[i];

	for(j = i+1;j  < ( n); ++j )
	{
	  elem_ji = j*(j-1)/2 + i;
	  sum -= adj(invclov[s][elem_ji]) * v1[j];
	}

	v1[i] = sum * diag_g[i];
      }

      /*# Overwrite column k of invclov */
      diag_invclov[s][k] = real(v1[k]);

      for(i = k+1;i  < ( n); ++i )
      {
	elem_ik = i*(i-1)/2 + k;
	invclov[s][elem_ik] = v1[i];
      }
    }
  }
  
  if ( DetP == YES )
  {
    logdet = sum(log_diag);
          }
  else
  {
    logdet = 0;
  }
  
            
  END_CODE();
})
