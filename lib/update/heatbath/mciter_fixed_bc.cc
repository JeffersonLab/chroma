/*! \file
 *  \brief One heatbath interation of updating the gauge field configuration
 */

#include "chromabase.h"
#include "util/gauge/reunit.h"
#include "update/heatbath/mciter_fixed_bc.h"
#include "update/heatbath/su3over.h"
#include "update/heatbath/su2_hb_update.h"

namespace Chroma 
{

  //! One heatbath interation of updating the gauge field configuration
  /*!
   * \ingroup heatbath
   *
   * Make one interation of updating the gauge field configuration:
   *      this consists of n_over overrelaxation sweeps followed
   *      by one heatbath sweep with nheat trials.
   * In the case of SU(3), for each link we loop over the 3 SU(2) subgroups.

   * Warning: this works only for Nc = 2 and 3 !

   * \param u        gauge field ( Modify )
   * \param S_g      gauge action ( Read )
   * \param hbp      heatbath parameters ( Read )
   * \param bc       the subset describing the boundary
   */

  void mciter_fixed_bc(multi1d<LatticeColorMatrix>& u, 
		       const LinearGaugeAction& S_g,
		       const HBParams& hbp,
		       const Subset& bc)
  {
    START_CODE();

    LatticeColorMatrix u0_mu ;
    LatticeColorMatrix u_mu_staple;
    int ntrials = 0;
    int nfails = 0;

    const Set& gauge_set = S_g.getSet();
    const int num_subsets = gauge_set.numSubsets();

    for(int iter = 0; iter <= hbp.nOver; ++iter)
    {
      for(int cb = 0; cb < num_subsets; ++cb)
      {
	for(int mu = 0; mu < Nd; ++mu)
	{
	  u0_mu = u[mu] ; // save the old gauge field
	  // Calculate the staple
	  {
	    typedef multi1d<LatticeColorMatrix>  P;
	    typedef multi1d<LatticeColorMatrix>  Q;

	    Handle< GaugeState<P,Q> > state(S_g.createState(u));

	    //staple 
	    S_g.staple(u_mu_staple, state, mu, cb);
	  }

	  if ( iter < hbp.nOver )
	  {
	    /* Do an overrelaxation step */
	    /*# Loop over SU(2) subgroup index */
	    for(int su2_index = 0; su2_index < Nc*(Nc-1)/2; ++su2_index)
	      {
		su3over(u[mu], u_mu_staple, su2_index, gauge_set[cb]);
		//u[mu][bc] = u0_mu
	      }
	  }
	  else
	  {
	    /* Do a heatbath step */
	    /*# Loop over SU(2) subgroup index */
	    for(int su2_index = 0; su2_index < Nc*(Nc-1)/2; ++su2_index)
	    {
	      int ntry = 0;
	      int nfail = 0;

	      su2_hb_update(u[mu], u_mu_staple,
			    Real(2.0/Nc),
			    su2_index, gauge_set[cb],
			    hbp.nmax());
	      //u[mu][bc] = u0_mu
	      ntrials += ntry;
	      nfails += nfail;
	    }
          
	    /* Reunitarize */
	    reunit(u[mu]);
	  }
	  // HERE I COPY BACK THE OLD FIELDS I NEED THE SETS
	  // I think I do not need the copy back after each SU(2) subgroup
	  // update.
	  u[mu][bc] = u0_mu ;

	  // If using Schroedinger functional, reset the boundaries
	  // NOTE: this routine resets all links and not just those under mu,cb
	  S_g.getGaugeBC().modify(u);

	}    // closes mu loop
      }      // closes cb loop
    }

  
    END_CODE();
  }

}  // end namespace Chroma
