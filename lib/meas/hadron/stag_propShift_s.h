
#ifndef STAG_PROPSHIFT_S_H
#define STAG_PROPSHIFT_S_H

/* Forward declarations for the shifting of the staggered propagator
 * index to co-incide with the corect delta for meson spectroscopy
 */

#include "chromabase.h"

int deltaToPropIndex(multi1d<int>& delta);

LatticeStaggeredPropagator shiftDeltaProp(multi1d<int>& delta, 
					  const LatticeStaggeredPropagator& src);

#endif
