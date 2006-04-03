/*  $Id: polyfrmnse_s.cc,v 3.0 2006-04-03 04:59:10 edwards Exp $ ($Date: 2006-04-03 04:59:10 $)                                                             */
/* PolyFrmNse                                                             */
/* This routine is specific to staggered fermions!                        */

/* Arguments:                                                             */

/* u           --  gauge field ( Read )                                   */
/* chi         --  pseudofermion field - result of inversion ( Write )    */
/* Ncb         --  Number of checkerboards. Ignored for now     (Read)    */
/* n_congrd    --  number of CG iterations ( Write )                      */

/* Local variables:                                                       */

/* eta         --  Gaussian random noise                                  */
/* psi         --  pseudofermion field - source of inversion              */

/* Functional description:                                                */

/* This routine calculates the pseudofermion field (chi) for the case     */
/* of polynomial evolution                                                */
/*                                                                        */
/*           chi = [ Q^r p^s ]^(-1) * psi                                 */
/*                                                                        */
/* Where:    Q   = mdagm_resc * M^dagM                                    */
/*           mdagm_resc  = PolyArgResc/(1+4*Nd*Nd*KappaMC^2)              */  
/*           p = p(Q) = (p^tilde)^dagger p^tilde                          */
/*           r = PolyPowNum       s = PolyPowDen                          */
/*                                                                        */
/*           psi = [Q^r p^(s-1) (p^tilde)^dagger] eta                     */
/*                                                                        */

include(types.mh)

SUBROUTINE(PolyFrmNse, u, chi, Ncb, n_congrd)

multi1d<LatticeColorMatrix> u(Nd);
multi1d<LatticeFermion> chi(Ncb);
int Ncb;
int n_congrd;

{
  include(COMMON_DECLARATIONS)
  LINEAR_OPERATOR(Apoly);
  PROTOTYPE(`Apoly', `DESC', `DESC', `DESC', `VAL', `VAL')
  LINEAR_OPERATOR(A);

  multi1d<LatticeFermion> eta(Ncb);
  multi1d<LatticeFermion> psi(Ncb);

  int r;
  int s;
  int t;
  int h;
  int cb;

  Double psi_norm;
  
  START_CODE();

  phfctr (u, FORWARD);     /* Turn ON Staggered Phases */

          ConsLinOp (A, u, KappaMC, FermAct);

  /* Prepare arguments of Apoly */
  r = PolyPowNum;
  s = PolyPowDen - 1;
  t = 1;
  h = -1;

  /* Construct the Apoly for source computation */
  CONSTRUCT_LINEAR_OPERATOR(Apoly, lpoly, A, r, s, t, h);  
  
  /* Fill eta with random gaussian noise: < eta_dagger * eta > = 1 */
  gaussian(eta);

  /* Compute source for inversions */
  Apoly (Apoly, eta, psi, Ncb, 1);

  /* Compute the norm of psi:  psi_norm = | psi |   */
  psi_norm = 0;
  for(cb = 0; cb < Ncb; ++cb)
    psi_norm += norm2(psi[cb]);
  psi_norm = sqrt(psi_norm);

  FREE_LINEAR_OPERATOR(Apoly);

  /* Prepare arguments of Apoly for inversion */
  s = PolyPowDen;
  t = 0;

  /* Construct the new Apoly */
  CONSTRUCT_LINEAR_OPERATOR(Apoly, lpoly, A, r, s, t, h);  

  /* Calculate the pseudofermion field */
  chi = 0;
  InvCG1 (Apoly, psi, chi, psi_norm, RsdCGMC, Ncb, n_congrd);

  FREE_LINEAR_OPERATOR(Apoly);
  DestLinOp (A, FermAct);
        
  phfctr (u, BACKWARD);    /* Turn OFF Staggered Phases */

  END_CODE();
}

