// -*- C++ -*-
// $Id: fuzwilp.h,v 1.1 2004-04-26 16:12:49 mcneile Exp $
/*! \file
 *  \brief Calculate ape-fuzzed Wilson loops
 */

#ifndef __fuzwilp_h__
#define __fuzwilp_h__

//! 
/*!
 * \ingroup glue
 *
*/
/* Computes time-like APE_fuzzed Wilson loops, including non-planar loops, */
/* Based on RE's SZIN code fuzwilp.m                                       */
/* This version makes APE-smeared links with no blocking                   */

/* Warning: this works only for Nc = 2 and 3 ! (Projection of */
/*                                              smeared/blocked links) */

/* Warning: this version assumes the space-like directions (perpendicular */
/*          to j_decay) to have equal length. */

/* u       -- gauge field ( Read ) */
/* j_decay -- 'time' direction for 'fuzzed' Wilson loops ( Read ) */
/* n_smear -- number of applying smearing to the gauge links ( Read ) */
/* sm_fact -- "smearing" factor = weight of old link w. r. to staples ( Read ) */
/* BlkAccu -- accuracy in fuzzy link projection ( Read ) */
/* BlkMax  -- maximum number of iterations in fuzzy link projection ( Read ) */

void fuzwilp(const multi1d<LatticeColorMatrix>& u, 
        int j_decay, int n_smear,
        const Real& sm_fact, const Real& BlkAccu, int BlkMax,
	XMLWriter& xml, const string& xml_group);

#endif
