#ifndef STAG_SCALARS_S_H
#define STAG_SCALARS_S_H

#define NUM_STAG_PROPS   8
#define NUM_STAG_PIONS   16


#include "chroma.h"

void
staggeredScalars(multi1d<LatticePropagator>& quark_props,
		      multi2d<DComplex>& scalar_corr_fn,
		      int j_decay);

#endif
