// -*- C++ -*-
// $Id: hadseqsrc_w.h,v 2.1 2005-09-26 04:48:35 edwards Exp $
/*! \file
 *  \brief Construct hadron sequential sources.
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
   * \param quark_propagators    array of quark propagators ( Read )
   * \param seq_src_prop         sequential source as propagator ( Write )
   * \param t_sink               time coordinate of the sink ( Read )
   * \param sink_mom             sink baryon momentum ( Read )
   * \param j_decay              direction of the exponential decay ( Read )
   * \param seq_src_name         string name of the sequential source ( Read )
   *
   * \return Sequential source propagator
   */

  LatticePropagator hadSeqSource(const multi1d<LatticePropagator>& quark_propagators, 
				 int t_sink, const multi1d<int>& sink_mom, 
				 int j_decay, 
				 const string& seq_src_name);

}  // end namespace Chroma


#endif
