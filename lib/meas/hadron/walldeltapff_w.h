// -*- C++ -*-
// $Id: walldeltapff_w.h,v 1.1 2004-05-31 05:12:32 edwards Exp $
/*! \file
 *  \brief Wall-sink delta^+ -> gamma+proton form-factors 
 *
 *  Form factors constructed from a quark and a wall sink
 */

#ifndef __walldeltapff_h__
#define __walldeltapff_h__

#include "util/ft/sftmom.h"

//! Compute contractions for current insertion 3-point functions.
/*!
 * \ingroup hadron
 *
 * This routine is specific to Wilson fermions!
 *
 * \param xml                buffer for writing the data ( Write )
 * \param u                  gauge fields (used for non-local currents) ( Read )
 * \param forw_u_prop        forward U quark propagator ( Read )
 * \param back_u_prop        backward D quark propagator ( Read )
 * \param forw_d_prop        forward U quark propagator ( Read )
 * \param back_d_prop        backward D quark propagator ( Read )
 * \param phases             fourier transform phase factors ( Read )
 * \param t0                 time coordinates of the source ( Read )
 * \param t_sink             time coordinates of the sink ( Read )
 */

void wallDeltaPFormFac(XMLWriter& xml,
		       const multi1d<LatticeColorMatrix>& u, 
		       const LatticePropagator& forw_u_prop,
		       const LatticePropagator& back_u_prop, 
		       const LatticePropagator& forw_d_prop,
		       const LatticePropagator& back_d_prop, 
		       const SftMom& phases,
		       int t0, int t_sink);

#endif
