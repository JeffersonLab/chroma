#ifndef ABS_MOM_REFRESH_H
#define ABS_MOM_REFRESH_H

#include "chromabase.h"
#include "update/field_state.h"
#include "update/molecdyn/abs_hamiltonian.h"

using namespace QDP;

// This is just a filler for the Abstract HMC to call
template< typename P, typename Q>
class AbsMomGaussianHeatbath {
public:
  virtual void operator()(AbsFieldState<P,Q>& state,
			  const AbsHamiltonian<P,Q>& H) const = 0;
};

#endif

