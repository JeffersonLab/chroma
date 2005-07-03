// $Id: su3_hb_sweep.cc,v 1.3 2005-07-03 16:10:11 edwards Exp $
/*! \file
 *  \brief Do one su(3) sweep
 */

#include "chromabase.h"
#include "update/heatbath/hb_params.h"
#include "update/heatbath/u_staple.h"
#include "update/heatbath/su2_hb_update.h"
#include "util/gauge/reunit.h"

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
   *
   * \param u          field to be updated ( Modify )
   */
  void su3_hb_sweep(multi1d<LatticeColorMatrix>& u, 
		    const HBParams& hbp) 
  {
    LatticeColorMatrix u_mu_staple;

    for(int cb=0;cb<2;cb++) 
    {
      for(int mu=0; mu<Nd; mu++) 
      {
	//staple 
	u_staple(u_mu_staple,u,mu,rb[cb],hbp);
	//start heat bath
	//reweight BetaMC by 2/Nc = 2/3 for su(3) 
	for(int su2_index=0;su2_index<(Nc*(Nc-1))/2;su2_index++)
	  su2_hb_update(u[mu],u_mu_staple,
			(2.0/Nc*hbp.beta())/hbp.xi(),
			su2_index,rb[cb],
			hbp.nmax());
	reunit(u[mu]);
      } // close mu loop
    }	
  }

}  // end namespace Chroma
