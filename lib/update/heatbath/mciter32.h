// -*- C++ -*-
// $Id: mciter32.h,v 3.1 2007-02-22 21:11:49 bjoo Exp $
/*! \file
 *  \brief One heatbath interation of updating the gauge field configuration
 */

#ifndef __mciter32_h__
#define __mciter32_h__

namespace Chroma {

//! One heatbath interation of updating the gauge field configuration
/*!
 * \ingroup heatbath
 *
 * Make one interation of updating the gauge field configuration:
 *      for Wilson or Symmanzik improved pure gauge action:
 *      this consists of n_over overrelaxation sweeps followed
 *      by one heatbath sweep with nheat trials.
 * In the case of SU(3), for each link we loop over the 3 SU(2) subgroups.
 *
 * Because of the option for the Symanzik improved action we need 2^(d+1)
 *      sublattices, deviding first into 2^d hypercubes and then
 *       checkerboarding those. 
 *
 *
 * Warning: this works only for Nc = 2 and 3 !
 *
 * \param u        gauge field ( Modify )
 * \param n_over   number of overrelaxation sweeps ( Read )
 * \param nheat    number of heatbath trials ( Read )
 * \param ntrials  total number of individual heatbath trials ( Modify )
 * \param nfails   total number of individual heatbath failures ( Modify ) 
 */

void mciter32(multi1d<LatticeColorMatrix>& u, 
	      int n_over, int nheat,
	      int& ntrials, int& nfails,
	      const Set& ss,
	      const multi3d<int>& neighsubl);

}  // end namespace Chroma

#endif

