/*! This function, converts a set of shifts stored in delta, into
 *  an index for the array of propagators. 
 */

#include "meas/hadron/stag_propShift_s.h"

int deltaToPropIndex(multi1d<int>& delta)
{

  if( Nd != 4 ) {
    QDPIO::cerr << "This routine only works for Nd=4, not Nd="<<Nd << endl;
    QDP_abort(1);
  }

  /* Paranoia -- Delta has to have Nd elements */
  if( delta.size() != Nd ) {
    QDPIO::cerr << "Delta has to have Nd elements as opposed to " << delta.size();
    QDP_abort(1);
  }

  /* Paranoia2 -- Delta can contain only 0 (no shift) and 1s (shift) */
  for( int mu = 0; mu < Nd; mu++) {
    if( delta[mu] != 0 && delta[mu] != 1 ) {
      QDPIO::cerr << "Delta must be made up of either zeros or ones. Element " << mu << " has value " << delta[mu] << endl;
      QDP_abort(1);
    }
  }

  /* The delta array is turned into a lexicographic index within
     the 4d hypercube */
  return(delta[0] + 2*( delta[1] + 2*(delta[2] + 2*delta[3])));
}

/*! Given an array of forward shifts (up to 1 in each dimension) and a
 *  propagator. This routine will carry out up to 4 shifts on the input,
 *  possibly 1 in each dimension
 */

LatticeStaggeredPropagator shiftDeltaProp(multi1d<int>& delta,
                                 const LatticeStaggeredPropagator& src)
{

  int mu;

  if( delta.size() != Nd ) {
    QDPIO::cerr << "Delta has to have Nd elements as opposed to " << delta.size();
    QDP_abort(1);
  }

  for( mu = 0; mu < Nd; mu++) {
    if( delta[mu] != 0 && delta[mu] != 1 ) {
      QDPIO::cerr << "Delta must be made up of either zeros or ones. Element " << mu << " has value " << delta[mu] << endl;
      QDP_abort(1);
    }
  }

  LatticeStaggeredPropagator ret_val = src;
  LatticeStaggeredPropagator tmp;

  for( mu = 0; mu < Nd; mu++) {
    if( delta[mu] == 1 ) {
      // This at the moment cannot occur without a temporary
      tmp = shift(ret_val, FORWARD, mu);
      ret_val = tmp;
    }
  }

  return ret_val;
}

