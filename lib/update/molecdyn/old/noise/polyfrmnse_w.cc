/* $Id: polyfrmnse_w.cc,v 3.0 2006-04-03 04:59:10 edwards Exp $ ($Date: 2006-04-03 04:59:10 $) */

/* This routine is specific to wilson fermions!                           */

/* Arguments:                                                             */

/* u           --  gauge field ( Read )                                   */
/* chi         --  pseudofermion field - result of inversion ( Write )    */
/* chi_norm    --  | chi | ( Write )                                      */
/* Ncb         --  Number of checkerboards. Ignored for now     (Read)    */
/* n_congrd    --  number of CG iterations ( Write )                      */

/* Local variables:                                                       */

/* Functional description:                                                */

/* This routine WILL calculate the pseudofermion field (chi) for the      */
/* case of polynomial evolution once it is implemented for Wilson         */
/* fermions.                                                              */

include(types.mh)

SUBROUTINE(PolyFrmNse, u, chi, Ncb, n_congrd)

multi1d<LatticeColorMatrix> u(Nd);
multi1d<LatticeFermion> chi(Ncb);
int Ncb;
int n_congrd;
{
  include(COMMON_DECLARATIONS)
  
  START_CODE();

  QDP_error_exit("Polynomial evolution not implemented for Wilson fermions");

  END_CODE();
}

