// -*- C++ -*-
// $Id: wallnuclff_w.h,v 1.3 2004-06-04 21:13:15 edwards Exp $
/*! \file
 *  \brief Wall-sink nucleon form-factors 
 *
 *  Form factors constructed from a quark and a wall sink
 */

#ifndef __wallnuclff_h__
#define __wallnuclff_h__

#include "util/ft/sftmom.h"
#include "meas/hadron/wallff_w.h"

//! Wall-sink nucleon-> gamma+nucleon form-factors
/*!
 * \ingroup hadron
 *
 * This routine is specific to Wilson fermions!
 *
 * \param form               Mega-structure holding form-factors ( Write )
 * \param u                  gauge fields (used for non-local currents) ( Read )
 * \param forw_u_prop        forward U quark propagator ( Read )
 * \param back_u_prop        backward D quark propagator ( Read )
 * \param forw_d_prop        forward U quark propagator ( Read )
 * \param back_d_prop        backward D quark propagator ( Read )
 * \param phases             fourier transform phase factors ( Read )
 * \param t0                 time coordinates of the source ( Read )
 * \param t_sink             time coordinates of the sink ( Read )
 */

void wallNuclFormFac(WallFormFac_formfacs_t& form,
		     const multi1d<LatticeColorMatrix>& u, 
		     const LatticePropagator& forw_u_prop,
		     const LatticePropagator& back_u_prop, 
		     const LatticePropagator& forw_d_prop,
		     const LatticePropagator& back_d_prop, 
		     const SftMom& phases,
		     int t0, int t_sink);

#endif
