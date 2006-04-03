
namespace Chroma {


include(types.mh)

/* + */
/* $Id: mespbg5p_w.cc,v 3.0 2006-04-03 04:59:04 edwards Exp $ ($Date: 2006-04-03 04:59:04 $) */

/* MESPBG5P - Calculates noise estimator for psi_bar_psi, psi_bar gamma_5 psi, */
/*            <psi_bar gamma_5 psi psi_bar gamma_5 psi>,  */
/*            <psi_bar gamma_5 psi (psi_bar gamma_5 psi)^dag> */

/* This routine is specific to Wilson fermions! */

/* u           -- gauge field ( Read ) */
/* n_hit       -- number of Gaussian hits ( Read ) */
/* pbp_st      -- estimator of psi_bar_psi from first hit ( Write ) */
/* n_congrd    -- Number of CG iteration ( Write ) */
/* - */
SUBROUTINE(MesPbg5p, u, nhit, pbp_st, n_congrd)

multi1d<LatticeColorMatrix> u(Nd);
int n_congrd;
Double pbp_st;
int nhit;

{ /* Local variables */
  include(COMMON_DECLARATIONS)

  Double pbp_st_m;
  DComplex pbg5p_st;
  DComplex pbg5p_sq_st;
  DComplex pbg5p_mdsq_st;
  LatticeFermion psi;
  LatticeFermion tmp;
  LatticeFermion eta;
  LatticeComplex aux_0;
  LatticeComplex aux_1;
  LatticeReal lrtrace_aux;
  LatticeReal lrtrace_tmp;
  LatticeComplex lctrace_aux;

  Double ddummy1;
  Double ddummy2;
  int ihit;

  START_CODE();

  phfctr (u, FORWARD);              /* ON */

      
  pbp_st_m = 0;
  pbg5p_st = 0;
  pbg5p_sq_st = 0;
  pbg5p_mdsq_st = 0;
  n_congrd = 0;

      
  for(ihit = 1; ihit <= nhit; ++ihit)
  {
    /* Fill with random gaussian noise such that < eta_dagger * eta > = 1 */
    gaussian(eta);

    /* If using Schroedinger functional, zero out the boundaries */
    if ( SchrFun > 0 )
    {
      FILLMASK(eta(0), lSFmaskF(0), ZERO);
      FILLMASK(eta(1), lSFmaskF(1), ZERO);
    }

    /* First source */
    /* aux_0 = W^(-1) * eta */
    /* The inverse of the un-preconditioned matrix is needed here. This */
    /* is given by Qprop! */
            psi = 0;
    tmp = eta;
    Qprop (u, tmp, KappaMC, RsdCGMC, psi, n_congrd);
    
    /* Chiral condensate = Tr [ eta_dag * psi ] */
        
    lrtrace_aux = real(trace(adj[eta[0]] * psi[0]));
    lrtrace_tmp = real(trace(adj[eta[1]] * psi[1]));

    lrtrace_aux += lrtrace_tmp;
    if ( ihit == 1 )
    {
      pbp_st = sum(lrtrace_aux);
      pbp_st_m += pbp_st;
    }
    else
      pbp_st_m += sum(lrtrace_aux);

        
    /* Scalar condensate = Tr [ eta_dag * gamma_5 * psi ] */
        
    SPIN_PRODUCT(psi(0),15,tmp(0));
    SPIN_PRODUCT(psi(1),15,tmp(1));
    aux_0 = trace(adj[eta] * tmp);

    lctrace_aux = aux_0[0];
    lctrace_aux += aux_0[1];
    pbg5p_st += sum(lctrace_aux);

            

    /* Fill with random gaussian noise such that < eta_dagger * eta > = 1 */
    gaussian(eta);

    /* Second source */
    /* psi = W^(-1) * eta */
    /* The inverse of the un-preconditioned matrix is needed here. This */
    /* is given by Qprop! */
            psi = 0;
    tmp = eta;
    Qprop (u, tmp, KappaMC, RsdCGMC, psi, n_congrd);
    
    /* Chiral condensate = Tr [ eta_dag * psi ] */
        
    lrtrace_aux = real(trace(adj[eta[0]] * psi[0]));
    lrtrace_tmp = real(trace(adj[eta[1]] * psi[1]));

    lrtrace_aux += lrtrace_tmp;
    pbp_st_m += sum(lrtrace_aux);

        
    /* Scalar condensate = Tr [ eta_dag * gamma_5 * psi ] */
        
    SPIN_PRODUCT(psi(0),15,tmp(0));
    SPIN_PRODUCT(psi(1),15,tmp(1));
    aux_1 = trace(adj[eta] * tmp);

    lctrace_aux = aux_1[0];
    lctrace_aux += aux_1[1];
    pbg5p_st += sum(lctrace_aux);

            

    /* Construct   <psi_bar gamma_5 psi  psi_bar gamma_5 psi > */
    /* this is equivalent to  < tr( eta2^dag * psi2 * eta1^dag * psi1 )> */
    /* Also construct  <psi_bar gamma_5 psi  (psi_bar gamma_5 psi)^dag > */
    /* which is  equivalent to  < tr( eta2^dag * psi2 * (eta1^dag * psi1^dag) ) > */
    
    lctrace_aux = aux_1[0] * aux_0[0];
    lctrace_aux += aux_1[1] * aux_0[1];
    pbg5p_sq_st += sum(lctrace_aux);

    lctrace_aux = aux_1[0] * adj(aux_0[0]);
    lctrace_aux += aux_1[1] * adj(aux_0[1]);
    pbg5p_mdsq_st += sum(lctrace_aux);

      }

      
  phfctr (u, BACKWARD);              /* OFF */

  pbp_st = pbp_st / TO_DOUBLE(vol_cb*2*Nc*Ns);
  ddummy1 = TO_DOUBLE(1) / TO_DOUBLE(nhit * 2*vol*Nc*Ns);
  ddummy2 = TO_DOUBLE(1) / TO_DOUBLE(nhit * vol*Nc*Nc*Ns*Ns);
  pbp_st_m = pbp_st_m * ddummy1;
  pbg5p_st = pbg5p_st * ddummy1;
  pbg5p_sq_st = pbg5p_sq_st * ddummy2;
  pbg5p_mdsq_st = pbg5p_sq_st * ddummy2;

  push(xml_out,"Condensates");
write(xml_out, "nhit", nhit);
write(xml_out, "pbp_st_m", pbp_st_m);
write(xml_out, "pbg5p_st", pbg5p_st);
write(xml_out, "pbg5p_sq_st", pbg5p_sq_st);
write(xml_out, "pbg5p_mdsq_st", pbg5p_mdsq_st);
pop(xml_out);

      
  END_CODE();
}

}  // end namespace Chroma
