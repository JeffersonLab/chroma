/* This routine is specific to Wilson fermions! */

/* Dummy version of subroutine for pseudofermionic force calculation       */
/* in the polynomial hybrid Monte Carlo algorithm.                         */

/* (S. Sint, 9/97)  						           */

/* polydsduf -- computes the derivative of the pseudofermionic action with */ 
/* respect to the link field                                               */

/* u -- gauge field           ( Read )   */
/*         |  dS      dS_f               */
/* ds_u -- | ----   + -----   ( Modify ) */
/*         |  dU       dU                */
/* psi -- pseudofermion       ( Read )   */

void polydsduf(multi1d<LatticeColorMatrix>& ds_u,
	       const multi1d<LatticeColorMatrix>& u,
	       const LatticeFermion& psi)
{
  START_CODE("subroutine");;
 
  QDP_error_exit("Algorithm not yet implemented for Wilson fermions", Algorithm);

  END_CODE("subroutine");;
}
