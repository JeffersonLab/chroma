// -*- C++ -*-
// $Id: fuzwilp.h,v 1.2 2004-04-27 20:29:07 edwards Exp $
/*! \file
 *  \brief Calculate ape-fuzzed Wilson loops
 */

#ifndef __fuzwilp_h__
#define __fuzwilp_h__

//! Calculate ape-fuzzed Wilson loops
/*!
 * \ingroup glue
 *
 * Computes time-like APE_fuzzed Wilson loops, including non-planar loops,
 *
 * This version makes APE-smeared links with no blocking as required
 *          for potential clculations
 *
 * Warning: This version is VERY Slow as it has non-recursive shifting
 *          of some link products
 *          Search for 'cap' on loop values to control no of loops
 * 	    calculated
 * Warning: this works only for Nc = 2 and 3 ! (Projection of
 *                                              smeared/blocked links)
 *
 * Warning: this version assumes the space-like directions (perpendicular
 *          to j_decay) to have equal length.
 *
 * \param u         gauge field ( Read )
 * \param j_decay   'time' direction for 'fuzzed' Wilson loops ( Read )
 * \param n_smear   number of applying smearing to the gauge links ( Read )
 * \param sm_fact   "smearing" factor = weight of old link w. r. to staples ( Read )
 * \param BlkAccu   accuracy in fuzzy link projection ( Read )
 * \param BlkMax    maximum number of iterations in fuzzy link projection ( Read ) 
 */

void fuzwilp(const multi1d<LatticeColorMatrix>& u, 
        int j_decay, int n_smear,
        const Real& sm_fact, const Real& BlkAccu, int BlkMax,
	XMLWriter& xml, const string& xml_group);

#endif
