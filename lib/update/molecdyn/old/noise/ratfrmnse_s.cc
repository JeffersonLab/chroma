/* RatFrmNse                                                              */
/* This is routine is specific for Staggered fermions!                    */

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
/*           chi =  n(Q)*[d(Q)]^(-1) * eta                                */
/*                                                                        */
/* Where:    Q   = mdagm_resc * M^dagM                                    */
/*           mdagm_resc  = PolyArgResc/(1+4*Nd*Nd*KappaMC^2)              */  
/*           d(Q) = (Q+rho_1)*(Q+rho_2)*...*(Q+rho_m)               */
/*           n(Q) = (Q+nu_1)*(Q+nu_2)*...  *(Q+nu_m)                */
/*           m  = HBRatDeg                                            */
/*	     							  */
/* The rational function n(x)/d(x) is the optimal rational */
/* approximation to the inverse square root of N(x)/D(x) which in */
/* turn is the optimal rational approximation to x^(-alpha).*/

/* To solve {n(Q)*[d(Q)]^(-1) * eta} the partial fraction expansion is */
/* used in combination with a multishift solver.*/

include(types.mh)

SUBROUTINE(RatFrmNse, u, psi, Ncb, n_congrd)

multi1d<LatticeColorMatrix> u(Nd);
multi1d<LatticeFermion> psi(Ncb);
int Ncb;
int n_congrd;

{
  include(COMMON_DECLARATIONS)

  LINEAR_OPERATOR(A);

  multi1d<LatticeFermion> eta(Ncb);
  multi2d<LatticeFermion> chi(Ncb, HBRatDeg);

  int i;
  int cb;
  Double eta_norm;
  
  multi1d<Real> RsdCG(HBRatDeg);

  START_CODE();

  phfctr (u, FORWARD);     /* Turn ON Staggered Phases */

  n_congrd=0;

/* Fill eta with random gaussian noise: < eta_dagger * eta > = 1 */

    gaussian(eta);

  eta_norm = 0;
  for (cb=0; cb<Ncb; ++cb) eta_norm += norm2(eta[cb]);
  eta_norm = sqrt(eta_norm);

      for (i=0; i<HBRatDeg; i++) {
    for (cb=0; cb<Ncb; ++cb) {
      chi[i][cb] = 0;
    }
    RsdCG[i] = RsdCGMD;
  }

    ConsLinOp (A, u, KappaMD, FermAct);

  /* Call the multishift solver.*/
  switch (InvType) {
  case CG_INVERTER:
    MInvCGm (A, eta, chi, eta_norm, HBDenRatRoots, HBRatDeg, 0, RsdCG, Ncb, n_congrd, OperEigVec, NOperEigDim);
    break;
  case MR_INVERTER:
    MInvMRm (A, eta, chi, eta_norm, HBDenRatRoots, HBRatDeg, 0, RsdCG, Ncb, n_congrd);
    break;
  case BICG_INVERTER:
    QDP_error_exit("Unsupported inverter type", InvType);
    break;
    
  default:
    QDP_error_exit("Unknown inverter type", InvType);
  }

  DestLinOp (A, FermAct);
  
  psi = 0;
  for (cb=0; cb<Ncb; ++cb)
    for (i=HBRatDeg-1; i>=0; i--) 
      psi[cb] += chi[cb][i] * HBNumRatRoots[i];

  /*Add final constant term.*/
  for (cb=0; cb<Ncb; ++cb)
    psi[cb] += eta[cb] * HBRatNorm;

        phfctr (u, BACKWARD);    /* Turn OFF Staggered Phases */
  END_CODE();
}




