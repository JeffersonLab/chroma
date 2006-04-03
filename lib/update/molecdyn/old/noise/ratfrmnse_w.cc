/* $Id: ratfrmnse_w.cc,v 3.0 2006-04-03 04:59:10 edwards Exp $ ($Date: 2006-04-03 04:59:10 $) */
/* RatFrmNse                                                              */
/* This is a dummy routine for Wilson fermions!                           */

/* Arguments:                                                             */

/* u           --  gauge field ( Read )                                   */
/* chi         --  pseudofermion field - result of inversion ( Write )    */
/* Ncb         --  Number of checkerboards. Ignored for now     (Read)    */
/* n_congrd    --  number of CG iterations ( Write )                      */

/* Local variables:                                                       */

/* eta         --  Gaussian random noise                                  */


/* Functional description:                                                */

/* This routine calculates the pseudofermion field (chi) for the case     */
/* of rational evolution                                                  */
/*                                                                        */
/*           chi =  d(Q)*[n(Q)]^(-1) * eta                                */
/*                                                                        */
/* Where:    Q   = mdagm_resc * M^dagM                                    */
/*           mdagm_resc  = PolyArgResc/(1+4*Nd*Nd*KappaMC^2)              */  
/*           d(Q) = (Q+rho_1)*(Q+rho_2)*...*(Q+rho_m/2)               */
/*           n(Q) = (Q+nu_1)*(Q+nu_2)*...  *(Q+nu_m/2)                */
/*           m/2  = PolyDegHalf                                            */
/*									  */
/*  The positive real numbers rho_i and nu_i are the negative roots of    */
/*  the polynomials in the numerator and the denominator respectively.    */
/*  Up to rounding errors the order of the monomials does not             */
/*  matter and we thus do not stick to a  particular ordering.            */



include(types.mh)

SUBROUTINE(RatFrmNse, u, chi, Ncb, n_congrd)

multi1d<LatticeColorMatrix> u(Nd);
LatticeFermion chi;
int Ncb;
int n_congrd;

{
  include(COMMON_DECLARATIONS)
 
  START_CODE();
 
  QDP_error_exit("Algorithm not yet implemented for Wilson fermions", Algorithm);

  END_CODE();
}

