/* RatMesFerm                                                            */
/* This routine is specific to staggered fermions!                         */

/* Arguments:                                                             */

/* u         --  gauge field (Read)                                       */
/* chi       --  pseudofermion field     (Read)                           */
/* w_ferm    --  Pseudo fermion energy   (Write)                          */
/* n_congrd  --  number of CG iterations (Write)                          */

/* Functional description:                                                */

/* This measures the pseudofermion contribution to the Hamiltonian        */
/* for the case of rational evolution (with polynomials N(x) and D(x),    */
/* of degree FRatDeg)                     */
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

  LINEAR_OPERATOR(A);

  multi1d<LatticeFermion> eta(FRatDeg);
  LatticeFermion psi;

  LatticeReal trace_aux;

  int Ncb;
  int i;

  multi1d<Real> RsdCG(FRatDeg);

  Double chi_norm;

  START_CODE();

  phfctr (u, FORWARD);     /* Turn ON Staggered Phases */

  Ncb=1;
  n_congrd=0;

  if ( FermiP == YES )
    {
                              
      ConsLinOp (A, u, Kappa, FermAct);
      
      /*Construct chi norm.*/
      chi_norm = norm2(chi);
      chi_norm = sqrt(chi_norm);
      
      for (i=0; i<FRatDeg; i++) RsdCG[i] = RsdCGMD;

      /* Call the multishift solver.*/
      switch (InvType) {
      case CG_INVERTER:
	MInvCGm (A, chi, eta, chi_norm, FDenRatRoots, FRatDeg, 0, RsdCG, Ncb, n_congrd, OperEigVec, NOperEigDim);
	break;
      case MR_INVERTER:
	MInvMRm (A, chi, eta, chi_norm, FDenRatRoots, FRatDeg, 0, RsdCG, Ncb, n_congrd);
	break;
      case BICG_INVERTER:
	QDP_error_exit("Unsupported inverter type", InvType);
	break;
	
      default:
	QDP_error_exit("Unknown inverter type", InvType);
      }
      DestLinOp (A, FermAct);
            
      psi = 0;
      for (i=FRatDeg-1; i>=0; i--) {
	psi += eta[i] * FNumRatRoots[i];
      }
      psi += chi * FRatNorm;

      /* Calculate: eta^dag * eta ---> w_ferm */
      trace_aux = real(trace(adj[chi] * psi));
      w_ferm = sum(trace_aux);
      
                        
      /* Calculate fermionic energy per pseudofermion degree of freedom */
      w_ferm = w_ferm / TO_DOUBLE(2*(vol/2)*Nc*Ns);
    }
  else
    {
      w_ferm = 0;
    }
  
  phfctr (u, BACKWARD);    /* Turn OFF Staggered Phases */
  
  END_CODE();
}

