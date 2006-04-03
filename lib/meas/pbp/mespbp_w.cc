// $Id: mespbp_w.cc,v 3.0 2006-04-03 04:59:04 edwards Exp $
/*! \file
 *  \brief  Calculates noise estimator for psi_bar_psi
 */

#error "Converted but not tested"

#include "chromabase.h"
#include "meas/pbp/ovpbg5p_w.h"


namespace Chroma {

//! Calculates noise estimator for psi_bar_psi
/*!
 * This routine is specific to Wilson fermions!
 *
 * u           -- gauge field ( Read )
 * psi_bar_psi -- chiral condensate  ( Write )
 * n_congrd    -- Number of CG iteration ( Write )
 * ichiral     -- chirality with no zero modes ( Read ) 
 */

void MesPbp(XMLWriter& xml_out,
	    Double& psi_bar_psi;
	    const WilsonTypeFermAct< multi1d<LatticeFermion> >& S_f,
	    Handle<const ConnectState> state,
	    const multi1d<Real>& Mass, 
	    enum InvType invType,
	    const multi1d<Real>& RsdCG,
	    int MaxCG,
	    int n_congrd,
	    int ichiral)
{
  START_CODE();

  const int numMass = Mass.size();

  LatticeFermion eta;
  LatticeFermion chi;
  LatticeFermion psi;
  LatticeReal trace_aux0;
  LatticeReal trace_aux1;
  multi1d<Real> Mass(numMass);
  multi1d<Real> RsdCG(numMass);
  int numMass;
  int nhit;
  int n_zero;

  phfctr (u, FORWARD);              /* ON */

  if ((FermAct == OVERLAP_POLE) || (FermAct == ZOLOTAREV_4D))
  {
    /* Use special overlap psi-bar-psi routine */
    n_congrd = 0;
    psi_bar_psi = 0;

    numMass = 1;
    nhit = 1;
    FILL(Mass, KappaMC);
    FILL(RsdCG, RsdCGMC);
    n_zero = - ichiral;
    OvPbg5p (u, n_zero, Mass, numMass, nhit, RsdCG);
  }
  else
  {
    /* Use convenional psi-bar-psi measurement */
    
    /* fill with random gaussian noise such that < eta_dagger * eta > = 1 */
    gaussian(eta);

    /* If using Schroedinger functional, zero out the boundaries */
    if ( SchrFun > 0 )
    {
      FILLMASK(eta(0), lSFmaskF(0), ZERO);
      FILLMASK(eta(1), lSFmaskF(1), ZERO);
    }

        
    /* psi = W^(-1) * eta */
    /* The inverse of the un-preconditioned matrix is needed here. This */
    /* is given by Qprop! */
    /* Note: Qprop trashes the source! Hence copy eta onto chi and use chi */
    /* as the source. */
    chi = eta;
    psi = 0;
    n_congrd = 0;
    Qprop (u, chi, KappaMC, RsdCGMC, psi, n_congrd);

    /* Chiral condensate = Tr [ eta_dag * psi ] */
    trace_aux0 = real(trace(adj[eta[0]] * psi[0]));
    trace_aux1 = real(trace(adj[eta[1]] * psi[1]));
    trace_aux0 += trace_aux1;
    psi_bar_psi = sum(trace_aux0);
        
    psi_bar_psi = psi_bar_psi/TO_DOUBLE(vol_cb*2*Nc*Ns);

  }

  phfctr (u, BACKWARD);              /* OFF */

  END_CODE();
}

}  // end namespace Chroma
