// -*- C++ -*-
// $Id: mesonseqsrc_w.h,v 1.3 2005-01-14 18:42:36 edwards Exp $
/*! \file
 *  \brief Construct meson sequential sources.
 */

#ifndef __mesonseqsrc_w_h__
#define __mesonseqsrc_w_h__

namespace Chroma {

//! Construct a meson sequential source.
/*!
 * \ingroup hadron
 *
 * This routine is specific to Wilson fermions!
 *
 *  delta(tz-tx) exp(i p.z) \gamma_5 G \gamma_5
 *
 * \param quark_propagator   input quark propagator ( Read ) 
 * \param seq_src_prop       sequential source as propagator ( Write ) 
 * \param t_sink             time coordinate of the sink ( Read ) 
 * \param sink_mom           sink pion momentum ( Read ) 
 * \param j_decay            direction of the exponential decay ( Read ) 
 * \param seq_src            the particular type of source ( Read )
 */
void mesonSeqSource(const LatticePropagator& quark_propagator,
		    LatticePropagator& seq_src_prop, 
		    int t_sink, multi1d<int>& sink_mom, 
		    int j_decay, int seq_src);

}  // end namespace Chroma

#endif
