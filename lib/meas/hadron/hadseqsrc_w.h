// -*- C++ -*-
// $Id: hadseqsrc_w.h,v 1.1 2005-03-07 02:55:20 edwards Exp $
/*! \file
 *  \brief Construct baryon sequential sources.
 */

#ifndef __hadseqsrc_w_h__
#define __hadseqsrc_w_h__

namespace Chroma 
{
  
  //! Construct hadron sequential sources
  /*!
   * \ingroup hadron
   *
   * This routine is specific to Wilson fermions!
   *
   * Construct hadronic sequential sources.
   *
   * \param quark_propagator_1   first (u) quark propagator ( Read )
   * \param quark_propagator_2   second (d) quark propagator ( Read )
   * \param quark_propagator_3   second (s) quark propagator ( Read )
   * \param seq_src_prop         sequential source as propagator ( Write )
   * \param t_sink               time coordinate of the sink ( Read )
   * \param sink_mom             sink baryon momentum ( Read )
   * \param j_decay              direction of the exponential decay ( Read )
   * \param seq_src_name         string name of the sequential source ( Read )
   *
   * \return Sequential source propagator
   */

  LatticePropagator hadSeqSource(const LatticePropagator& quark_propagator_1, 
				 const LatticePropagator& quark_propagator_2,
				 const LatticePropagator& quark_propagator_3,
				 int t_sink, const multi1d<int>& sink_mom, 
				 int j_decay, 
				 const string& seq_src_name);

}  // end namespace Chroma


#endif
