namespace Chroma {


/* $Id: mespbp_s.cc,v 3.0 2006-04-03 04:59:04 edwards Exp $ ($Date: 2006-04-03 04:59:04 $) */

/* This routine is specific to staggered fermions! */

/* Calculates psi_bar_psi */

/* u           -- gauge field ( Read ) */
/* psi_bar_psi -- chiral condensate  ( Write ) */
/* n_congrd    -- Number of CG iteration ( Write ) */
/* ichiral     -- not used here! ( Read ) */

include(types.mh)

SUBROUTINE(MesPbp, u, psi_bar_psi, n_congrd, ichiral)

multi1d<LatticeColorMatrix> u(Nd);
int n_congrd;
int ichiral;
Double psi_bar_psi;
{ /* Local Variables */
  include(COMMON_DECLARATIONS)

  LatticeStaggeredFermion eta;
  LatticeStaggeredFermion aux;
  LatticeStaggeredFermion tmp;
  LatticeStaggeredFermion psi;
  
  Double aux_norm;
  
  START_CODE();
  
  phfctr (u, FORWARD);              /* ON */
  
    
  /* Fill aux with random gaussian noise such that < aux_dagger * aux > = 1 */
  /* for even sites */
  gaussian(aux);
  
    
  /* Fill eta with random gaussian noise such that < eta_dagger * eta > = 1 */
  /* for odd sites */
  gaussian(eta);
  
  /* For Schroedinger functional: mask out noise on boundaries */
  if ( SchrFun > 0 )
  {
    FILLMASK(aux, lSFmaskF(0), ZERO);
    FILLMASK(eta, lSFmaskF(1), ZERO);
  }

    
  /* tmp = D'_dag*eta_e */
  dslash (u, eta, tmp, MINUS, 1);
  
  /* aux = aux + KappaMC*tmp = M_dag * "Guassian noise" */
  aux += tmp * KappaMC;
  
      
  /* computing the norm of aux  aux_norm = | aux | */
  aux_norm = norm2(aux);
  aux_norm = sqrt(aux_norm);
  
    
  /* psi = (M_dag*M)^(-1) * aux */
  psi = 0;
  invert (u, aux, psi, aux_norm, KappaMC, RsdCGMC, 1, n_congrd);
  
  /* Chiral condensate = Tr [ psi * psi_dag ] = Sum | psi |^2 */
  psi_bar_psi = norm2(psi);
  
      
  psi_bar_psi = WORD_VALUE(WORD_psi_bar_psi,TWO) * TO_DOUBLE(KappaMC) *
                psi_bar_psi / TO_DOUBLE(vol_cb*Nc);
  
  phfctr (u, BACKWARD);              /* OFF */
  
  END_CODE();
}

}  // end namespace Chroma
