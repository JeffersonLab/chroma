// -*- C++ -*-
/*! \file
 *  \brief Wall-sink rho-> gamma+rho form-factors 
 *
 *  Form factors constructed from a quark and a wall sink
 */

#ifndef __wallrhoff_h__
#define __wallrhoff_h__

#include "util/ft/sftmom.h"
#include "meas/hadron/wallff_w.h"

namespace Chroma {

//! Wall-sink rho-> gamma+rho form-factors 
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
 * \param u_x2               forward U quark propagator evaluated at sink  ( Read )
 * \param d_x2               forward D quark propagator evaluated at sink  ( Read )
 * \param phases             fourier transform phase factors ( Read )
 * \param t0                 time slice of the source ( Read )
 * \param wall_source        true if using a wall source ( Read )
 */

void wallRhoFormFac(WallFormFac_formfacs_t& form,
		    const multi1d<LatticeColorMatrix>& u, 
		    const LatticePropagator& forw_u_prop,
		    const LatticePropagator& back_u_prop, 
		    const LatticePropagator& forw_d_prop,
		    const LatticePropagator& back_d_prop, 
		    const Propagator& u_x2,
		    const Propagator& d_x2,
		    const SftMom& phases,
		    int t0,
		    bool wall_source);

}  // end namespace Chroma

#endif
