// $Id: ldumul_w.cc,v 3.0 2006-04-03 04:58:50 edwards Exp $

#error "NOT FULLY CONVERTED"

/*# This routine is specific to Wilson fermions! */

/*# LDUMUL - Performs the operation */

/*#    chi <-   (L + D + L^dag) . psi */

/*# where   */
/*#   L       is stored as a lower triangular matrix (no diagonal) */
/*#   diag(L) is stored separately as a real. */

/*# Arguments: */

/*#  L	       Lower tri. mat. (no diag)       (Read) */
/*#  Diag_L    diag(L)   	               (Read) */
/*#  Psi       Pseudofermion Source   	       (Read) */
/*#  Chi       Pseudofermion output    	       (Write) */
include(types.mh)

SUBROUTINE(ldumul, L, diag_L, psix, chix)

LATTICE_TRIANGULAR(L);
LATTICE_DIAG_TRIANGULAR(diag_L);
LatticeFermion psix;
LatticeFermion chix;
{ /* Local variables */
  include(COMMON_DECLARATIONS)

  multi2d<LatticeComplex> psi(n, 2);
  multi2d<LatticeComplex> chi(n, 2);

  unsigned i;
  unsigned j;
  unsigned n;
  unsigned s;
  unsigned elem_ij;
  unsigned elem_ji;
  
  START_CODE();
  
  n = 2*Nc;
  
  if ( Ns != 4 )
    QDP_error_exit("code requires Ns == 4", Ns);
  
    
  psi = CAST(psix);

  for(s = 0; s < 2; ++s)
  {
    for(i = 0; i < n; ++i)
    {
      chi[s][i] = diag_L[s][i] * psi[s][i];

      for(j = 0; j < i; ++j)
      {
	elem_ij = i*(i-1)/2 + j;

	chi[s][i] += L[s][elem_ij] * psi[s][j];
      }

      for(j = i+1; j < n; ++j)
      {
	elem_ji = j*(j-1)/2 + i;

	chi[s][i] += adj(L[s][elem_ji]) * psi[s][j];
      }
    }
  }
  
  chix = CAST(chi);

    
  END_CODE();
}
