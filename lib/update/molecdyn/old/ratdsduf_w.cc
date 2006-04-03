/* $Id: ratdsduf_w.cc,v 3.0 2006-04-03 04:59:10 edwards Exp $ ($Date: 2006-04-03 04:59:10 $) */
/* ratdsduf								   */

/* This is a dummy routine for Wilson fermions!                            */

/* subroutine for pseudofermionic force calculation                        */
/* in the rational hybrid Monte Carlo and related algorithms.              */
/* (S. Sint, 4/98)  						           */

/* The pseudofermionic action is of the form  S_pf = chi^dag N*D^(-1) chi  */
/* where N(Q)=n^2(Q) and D(Q)=d^2(Q) are polynomials in the rescaled       */
/* staggered kernel: 							   */

/*     Q   = mdagm_resc * M^dagM                                           */     
/*     M^dagM      = 1 - Kappa^2 D'  D'                                    */
/*     mdagm_resc  = PolyArgResc/(1+4*Nd^2*KappaMD^2)                      */

/* Note that the leading coefficient of both N(Q) and D(Q) is taken        */
/* to be 1. The function x^(-alpha) is approximated by c*N(x)/D(x)         */
/* with c=PolyNorm^2. This has to be taken into account when comparing     */
/* the force with the standard HMC routine dsduf for alpha=1     	   */

/* ratdsduf -- computes the derivative of the pseudofermionic action with  */ 
/* respect to the link field                                               */


/* u -- gauge field           ( Read )   */

/*         |  dS      dS_pf              */
/* ds_u -- | ----   + -----   ( Modify ) */
/*         |  dU       dU                */

/* psi -- pseudofermion       ( Read )   */
/* eta -- n*d^(-1)*psi        ( Modify ) */
/*     =>  eta^dag eta = S_pf            */
/* (can be used for energy measurement)  */

void ratdsduf(multi1d<LatticeColorMatrix>& ds_u,
	      const multi1d<LatticeColorMatrix>& u,
	      const LatticeFermion& psi,
	      const LatticeFermion& eta,
	      int n_congrd)
{
  START_CODE("subroutine");;
 
  QDP_error_exit("Algorithm not yet implemented for Wilson fermions", Algorithm);


  END_CODE("subroutine");;
}
