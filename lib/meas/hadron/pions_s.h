#ifndef PIONS_FOLLANA_S_H
#define PIONS_FOLLANA_S_H

#define NUM_STAG_PROPS   8
#define NUM_STAG_PIONS   16


#include "chroma.h"


void 
pions(multi1d<LatticeStaggeredPropagator>& quark_props,
      multi2d<DComplex>& pion_corr_fn,
      int j_decay) ;

#endif
