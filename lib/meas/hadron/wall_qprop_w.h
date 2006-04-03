// -*- C++ -*-
// $Id: wall_qprop_w.h,v 3.0 2006-04-03 04:59:01 edwards Exp $
/*! \file
 *  \brief Construct a wall-sink propagator
 */

#ifndef __wall_qprop_w_h__
#define __wall_qprop_w_h__

#include "util/ft/sftmom.h"

namespace Chroma {

//! Construct a wall-sink propagator:
/*!
 * \ingroup hadron
 *
 * This routine is specific to Wilson fermions!
 *
 * Each time slice will have one non-zero entry,
 * with the propagator summed over the entire time slice.
 *
 * \param wall_quark_prop         wall-sink quark propagator ( Write )
 * \param quark_propagator        quark propagator ( Read )
 * \param phases                  object holds list of momenta and Fourier phases ( Read )
 */

void wall_qprop(LatticePropagator& wall_quark_prop, 
		const LatticePropagator& quark_propagator, 
		const SftMom& phases);

}  // end namespace Chroma

#endif
