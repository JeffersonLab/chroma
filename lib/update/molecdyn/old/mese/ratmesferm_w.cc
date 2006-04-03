/* $Id: ratmesferm_w.cc,v 3.0 2006-04-03 04:59:10 edwards Exp $ ($Date: 2006-04-03 04:59:10 $) */
/* RatMesFerm                                                            */
/* This is a dummy routine for Wilson fermions!                          */

/* Arguments:                                                             */

/* u         --  gauge field (Read)                                       */
/* chi       --  pseudofermion field     (Read)                           */
/* w_ferm    --  Pseudo fermion energy   (Write)                          */
/* n_congrd  --  number of CG iterations (Write)                          */

/* Functional description:                                                */

/* This measures the pseudofermion contribution to the Hamiltonian        */
/* for the case of rational evolution (with polynomials N(x) and D(x),    */
/* of degree 2*PolyDegHalf and leading coefficient 1)                     */
/*                                                                        */
/* w_ferm = chi_dag * N(A)*D(A)^(-1)* chi                                 */

/* where A is the rescaled staggered kernel c*MdagM with                  */
/* c=PolyArgResc/(1+4*KappaMD^2*Nd^2)			                  */
/*                                                                        */
/* In the end w_ferm is normalized per pseudofermion degree of freedom    */

include(types.mh)
SUBROUTINE(RatMesFerm, u, chi, Kappa, w_ferm, n_congrd)

multi1d<LatticeColorMatrix> u(Nd);
LatticeFermion chi;
Real Kappa;
Double w_ferm;
int n_congrd;

{ /* Local variables */
  include(COMMON_DECLARATIONS)


  START_CODE();
 
  QDP_error_exit("Algorithm not yet implemented for Wilson fermions", Algorithm);

  END_CODE();
}

