// -*- C++ -*-
// $Id: barseqsrc_w.h,v 1.3 2005-01-14 18:42:36 edwards Exp $
/*! \file
 *  \brief Construct baryon sequential sources.
 */

#ifndef __barseqsrc_w_h__
#define __barseqsrc_w_h__

namespace Chroma {

//! Construct baryon sequential sources
/*!
 * \ingroup hadron
 *
 * This routine is specific to Wilson fermions!
 *
 * Construct baryon sequential sources.
 *
 * Note: only equal baryons at source and sink, for now.
 *
 * \param quark_propagator_1   first (u) quark propagator ( Read )
 * \param quark_propagator_2   second (d) quark propagator ( Read )
 * \param seq_src_prop         sequential source as propagator ( Write )
 * \param t_sink               time coordinate of the sink ( Read )
 * \param sink_mom             sink baryon momentum ( Read )
 * \param j_decay              direction of the exponential decay ( Read )
 * \param seq_src              flag indicating the type of the sequential source ( Read ) 
 */

void barSeqSource(const LatticePropagator& quark_propagator_1, 
		  const LatticePropagator& quark_propagator_2,
		  LatticePropagator& seq_src_prop, 
		  int t_sink, const multi1d<int>& sink_mom, 
		  int j_decay, int seq_src);

}  // end namespace Chroma


#endif
