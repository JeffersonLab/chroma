// -*- C++ -*-
// $Id: wallpionff_w.h,v 1.2 2004-04-11 05:35:12 edwards Exp $
/*! \file
 *  \brief Wall-sink pion form-factors 
 *
 *  Form factors constructed from a quark and a wall sink
 */

#ifndef __wallpionff_h__
#define __wallpionff_h__

#include "util/ft/sftmom.h"

//! Compute contractions for current insertion 3-point functions.
/*!
 * \ingroup hadron
 *
 * This routine is specific to Wilson fermions!
 *
 * \param xml                xml output ( Modify )
 * \param u                  gauge fields (used for non-local currents) ( Read )
 * \param forw_prop          forward quark propagator ( Read )
 * \param back_prop          backward quark propagator ( Read )
 * \param phases             fourier transform phase factors ( Read )
 * \param t0                 cartesian coordinates of the source ( Read )
 */

void wallPionFormFac(XMLWriter& xml,
		     const multi1d<LatticeColorMatrix>& u, 
		     const LatticePropagator& forw_prop,
		     const LatticePropagator& back_prop, 
		     const SftMom& phases,
		     int t0, int t_sink);

#endif
