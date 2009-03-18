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
  //
  // spin j = 3/2
  //

  void baryon_class7_s(
		LatticeStaggeredPropagator & quark_propagator_in_a,
		LatticeStaggeredPropagator & quark_propagator_in_b,
		LatticeStaggeredPropagator & quark_propagator_in_c,
		multi1d<Complex> & barprop,
		multi1d<int> & t_source,
		int j_decay, int bc_spec)  ; 


  void baryon_class7_NLT_s(
		LatticeStaggeredPropagator & quark_propagator_in_a,
		LatticeStaggeredPropagator & quark_propagator_in_b,
		LatticeStaggeredPropagator & quark_propagator_in_c,
		multi1d<LatticeColorMatrix> & u , 
		multi1d<Complex> & barprop,
		multi1d<int> & t_source,
		int j_decay, int bc_spec)  ; 

  void baryon_class4_s(
		LatticeStaggeredPropagator & quark_propagator_in_a,
		LatticeStaggeredPropagator & quark_propagator_in_b,
		LatticeStaggeredPropagator & quark_propagator_in_c,
		multi1d<Complex> & barprop,
		multi1d<int> & t_source,
		int j_decay, int bc_spec)  ;



}  // end namespace Chroma

#endif
