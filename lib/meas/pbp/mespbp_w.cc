// $Id: mespbp_w.cc,v 3.1 2009-04-21 01:21:42 eneil Exp $
/*! \file
 *  \brief  Calculates noise estimator for psi_bar_psi
 */

#include "chromabase.h"
#include "fermact.h"
#include "meas/pbp/mespbp_w.h"


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

void MesPbp(
		Handle< SystemSolver<LatticeFermion> > qprop,
		Handle< FermState<LatticeFermion, multi1d<LatticeColorMatrix>,
			multi1d<LatticeColorMatrix> > > state,
		const multi1d<Real>& Mass,
		const int ichiral,
		XMLWriter& xml_out,
		const string& xml_group,
		const string& FermAct)
{

  START_CODE();

  const int numMass = Mass.size();

  // Grab the links from the state
  const multi1d<LatticeColorMatrix>& u = state->getLinks();

  LatticeFermion eta;
  LatticeFermion chi;
  LatticeFermion psi;
  LatticeReal trace_aux0;
  LatticeReal trace_aux1;
  Double psi_bar_psi;
  
  push(xml_out, xml_group);

// I don't know what phfctr is...must be in that mysterious header file!
//  phfctr (u, FORWARD);              /* ON */

  if ((FermAct == "OVERLAP_POLE") || (FermAct == "ZOLOTAREV_4D"))
  {
	QDPIO::cerr << "mespbp_w.cc: Action " << FermAct << " unsupported." << endl;
	QDP_abort(1);

    /* Use special overlap psi-bar-psi routine */
/*    n_congrd = 0;
    psi_bar_psi = 0;

    numMass = 1;
    nhit = 1;
    FILL(Mass, KappaMC);
    FILL(RsdCG, RsdCGMC);
    n_zero = - ichiral;
    OvPbg5p (u, n_zero, Mass, numMass, nhit, RsdCG);
 */
   }
  else
  {
    /* Use convenional psi-bar-psi measurement */
    
    /* fill with random gaussian noise such that < eta_dagger * eta > = 1 */
    gaussian(eta);

    /* psi = W^(-1) * eta */
    /* The inverse of the un-preconditioned matrix is needed here. This */
    /* is given by Qprop! */
    /* Note: Qprop trashes the source! Hence copy eta onto chi and use chi */
    /* as the source. */
    chi = eta;
    psi = zero;
    int n_congrd = 0;
//    Qprop (u, chi, KappaMC, RsdCGMC, psi, n_congrd);
	SystemSolverResults_t res = (*qprop)(psi, chi);
	n_congrd = res.n_count;

    /* Chiral condensate = Tr [ eta_dag * psi ] */
	  trace_aux0 = real(trace(adj(eta) * psi));
	  psi_bar_psi = sum(trace_aux0);
	
	// Normalization
	Real norm = 1.0 / Real(2 * QDP::Layout::vol());
	norm /= (Ns * Nc);
        
    psi_bar_psi = psi_bar_psi * norm;
	
	// Write out results
	push(xml_out, "elem");
	write(xml_out, "pbp", psi_bar_psi);
	pop(xml_out);

  }

//  phfctr (u, BACKWARD);              /* OFF */

  END_CODE();
}

}  // end namespace Chroma
