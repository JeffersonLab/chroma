/* ratdsduf								   */

/* This routine is specific to staggered fermions!                         */

/* subroutine for pseudofermionic force calculation                        */
/* in the rational hybrid Monte Carlo and related algorithms.              */
/* (S. Sint, 4/98, revised M.Clark 02/03)				   */

/* The pseudofermionic action is of the form  S_pf = chi^dag N*D^(-1) chi  */
/* where N(Q) and D(Q) are polynomials in the rescaled staggered kernel:   */

/*     Q   = mdagm_resc * M^dagM                                           */  
/*     M^dagM      = 1 - Kappa^2 D'  D'                                    */
/*     mdagm_resc  = PolyArgResc/(1+4*Nd^2*KappaMD^2)                      */

/* The function x^(-alpha) is approximated by c*N(x)/D(x)         */
/* with c=FRatNorm. This has to be taken into account when comparing     */
/* the force with the standard HMC routine dsduf for alpha=1.  To solve */
/* ( N(Q)/D(Q) )psi the rational function is reformed as a partial fraction */
/* expansion and a multishift solver used to find the solution. */

/* ratdsduf -- computes the derivative of the pseudofermionic action with  */ 
/* respect to the link field.  Uses the partial fraction expansion of */
/* FRatRoots.  */


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
	      const multi1d<LatticeFermion>& psi,
	      const multi1d<LatticeFermion>& eta,
	      int n_congrd)
{
  multi1d<Real> RsdCG(FRatDeg);
  Real root_aux;
  int i;
  int cb;
  int mu;
  int Ncb;
  Real c;
  Double psi_norm;
  Real tmp;

  LINEAR_OPERATOR(A);
  multi2d<LatticeFermion> chi(Ncb, FRatDeg);
  multi1d<LatticeColorMatrix> ds_u_temp(Nd);

  START_CODE("subroutine");;

  // All fields live on 1 checkerboard only:
  Ncb=1;
    
  /*Construct linear operator.*/
  ConsLinOp (A, u, KappaMD, FermAct);
  tmp = 0;

  /*Construct psi norm.*/
  psi_norm = 0;
  for (cb=0; cb<Ncb; ++cb) psi_norm += norm2(psi[cb]);
  psi_norm = sqrt(psi_norm);

  for (i=0; i<FRatDeg; i++)
    RsdCG[i] = RsdCGMD;

  /* Call the multishift solver.*/
  switch (InvType) {
  case CG_INVERTER:
    MInvCGm (A, psi, chi, psi_norm, FDenRatRoots, FRatDeg, 0, RsdCG, Ncb, n_congrd, OperEigVec, NOperEigDim);
    break;
  case MR_INVERTER:
    MInvMRm (A, psi, chi, psi_norm, FDenRatRoots, FRatDeg, 0, RsdCG, Ncb, n_congrd);
    break;
  case BICG_INVERTER:
    QDP_error_exit("Unsupported inverter type", InvType);
    break;
    
  default:
    QDP_error_exit("Unknown inverter type", InvType);
  }
  DestLinOp (A, FermAct);
  
  /*Calculate the force and add it to the pure gauge contribution.*/
  for (cb=0; cb<Ncb; ++cb) 
  {
    for (i=FRatDeg-1; i>=0; i--) 
    {
      ds_u_temp = 0;
      StDsDuf (u, ds_u_temp, chi[i][cb]);
      for (mu=0; mu<Nd; mu++) 
      {
	ds_u[mu][0] += ds_u_temp[mu][0] * FNumRatRoots[i];
	ds_u[mu][1] += ds_u_temp[mu][1] * FNumRatRoots[i];
      } 
    }
  }
  
  END_CODE("subroutine");;
}
