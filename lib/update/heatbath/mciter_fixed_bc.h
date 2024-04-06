// -*- C++ -*-
/*! \file
 *  \brief One heatbath interation of updating the gauge field configuration
 */

#ifndef __mciter_fixed_bc_h__
#define __mciter_fixed_bc_h__

#include "actions/gauge/gaugeacts/wilson_gaugeact.h"
#include "update/heatbath/hb_params.h"

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
		       const Subset& bc);

}  // end namespace Chroma

#endif
