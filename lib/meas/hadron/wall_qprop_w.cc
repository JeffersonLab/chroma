// $Id: wall_qprop_w.cc,v 3.0 2006-04-03 04:59:01 edwards Exp $
/*! \file
 *  \brief Construct a wall-sink propagator
 */

#include "chromabase.h"
#include "util/ft/sftmom.h"
#include "meas/hadron/wall_qprop_w.h"

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
		const SftMom& phases)
{
  START_CODE();

  // Length of lattice in decay direction
  int length  = phases.numSubsets();
  int j_decay = phases.getDir();

  // Project propagator onto zero momentum: Do a slice-wise sum.
  multi1d<DPropagator> dprop_slice = sumMulti(quark_propagator, phases.getSet());
  
  // Now create the mask for 1 site per time slice
  LatticeBoolean lbmask = true;
  for(int mu = 0; mu < Nd; ++mu)
    if( mu != j_decay )
      lbmask &= (Layout::latticeCoordinate(mu) == 0);

  // Now copy onto the lattice
  LatticeInteger my_coord = Layout::latticeCoordinate(j_decay);

  wall_quark_prop = zero;
  for(int t = 0; t < length; ++t)
    copymask(wall_quark_prop, 
	     LatticeBoolean(lbmask & (my_coord == t)),
	     LatticePropagator(dprop_slice[t]));
            
  END_CODE();
}

}  // end namespace Chroma

