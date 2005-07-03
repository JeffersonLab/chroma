// $Id: mciter.cc,v 1.5 2005-07-03 16:10:11 edwards Exp $
/*! \file
 *  \brief One heatbath interation of updating the gauge field configuration
 */

#include "chromabase.h"
#include "util/gauge/reunit.h"
#include "update/heatbath/mciter.h"
#include "update/heatbath/su3over.h"
#include "update/heatbath/su2_hb_update.h"
#include "update/heatbath/u_staple.h"

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
   * \param n_over   number of overrelaxation sweeps ( Read )
   * \param nheat    number of heatbath trials ( Read )
   * \param ntrials  total number of individual heatbath trials ( Modify )
   * \param nfails   total number of individual heatbath failures ( Modify ) 
   */

  void mciter(multi1d<LatticeColorMatrix>& u, 
	      const HBParams& hbp)
  {
    START_CODE();

    LatticeColorMatrix u_mu_staple;
    int ntrials = 0;
    int nfails = 0;

    for(int iter = 0; iter <= hbp.nOver; ++iter)
    {
      for(int cb = 0; cb < 2; ++cb)
	for(int mu = 0; mu < Nd; ++mu)
	{
	  //staple 
	  u_staple(u_mu_staple, u, mu, rb[cb], hbp);

	  if ( iter < hbp.nOver )
	  {
	    /* Do an overrelaxation step */
	    /*# Loop over SU(2) subgroup index */
	    for(int su2_index = 0; su2_index < Nc*(Nc-1)/2; ++su2_index)
	      su3over(u[mu], u_mu_staple, su2_index, rb[cb]);
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
			    (2.0/Nc*hbp.beta())/hbp.xi(),
			    su2_index, rb[cb],
			    hbp.nmax());

//	      su3hb(u[mu], u_mu_staple, su2_index, nheat, ntry, nfail, rb[cb]);
	      ntrials += ntry;
	      nfails += nfail;
	    }
          
	    /* Reunitarize */
	    reunit(u[mu]);

	  }

#if 0
	  /* If using Schroedinger functional, reset the boundaries */
	  if ( SchrFun > 0 )
	  {
	    copymask(u[mu][cb], lSFmask[mu][cb], SFBndFld[mu][cb], REPLACE);
	  }
#endif
	}        /* closes mu loop */
    }

  
    END_CODE();
  }

}  // end namespace Chroma
