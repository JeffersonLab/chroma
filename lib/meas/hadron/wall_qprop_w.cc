// $Id: wall_qprop_w.cc,v 1.1 2003-12-16 04:34:32 edwards Exp $
/*! \file
 *  \brief Construct a wall-sink propagator
 */

#include "chromabase.h"
#include "util/ft/sftmom.h"
#include "meas/hadron/wall_qprop_w.h"

using namespace QDP;

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
  START_CODE("wall_qprop");

  // Length of lattice in decay direction
  int length  = ss.numSubsets();
  int j_decay = ss.getDir();

  // Project propagator onto zero momentum: Do a slice-wise sum.
  multi1d<DPropagator> dprop_slice = sumMulti(quark_propagator, ss.getSubset());
  
  // Now create the mask for 1 site per time slice
  LatticeBoolean lbtmp1 = true;
  for(int mu = 0; mu < Nd; ++mu)
    if( mu != j_decay )
      lbtmp1 &= (Layout::latticeCoordinate(mu) == 0);

  for(int t = 0; t < length; ++t)
    wall_quark_prop = where(lbtmp1 & (Layout::latticeCoordinate()[j_decay] == t),
			    LatticePropagator(dprop_slice[t]),
			    LatticePropagator(zero));
            
  END_CODE("wall_qprop");
}

