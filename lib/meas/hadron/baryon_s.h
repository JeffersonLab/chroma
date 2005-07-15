#ifndef baryon_s_h
#define baryon_s_h

#include "chromabase.h"

namespace Chroma {

void baryon_s(LatticeStaggeredPropagator & quark_propagator_in, 
	      multi1d<Complex> & barprop,
	      multi1d<int> & t_source,
	      int j_decay, int bc_spec) ;


  void baryon_s(
		LatticeStaggeredPropagator & quark_propagator_in_a, 
		LatticeStaggeredPropagator & quark_propagator_in_b, 
		LatticeStaggeredPropagator & quark_propagator_in_c, 
		multi1d<Complex> & barprop,
		multi1d<int> & t_source,
		int j_decay, int bc_spec) ;


  void baryon_local_s(LatticeStaggeredPropagator & quark_propagator_in, 
		multi1d<Complex> & barprop,
		multi1d<int> & t_source,
		      int j_decay, int bc_spec) ;

}  // end namespace Chroma

#endif
