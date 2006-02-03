#ifndef UTIL_COMPUTE_MESON_S__H 
#define UTIL_COMPUTE_MESON_S__H 

namespace Chroma { 

  void compute_vary_ps(
		      LatticeStaggeredPropagator & quark_propagator_Lsink_Lsrc,
		      LatticeStaggeredPropagator & quark_propagator_Fsink_Lsrc,
		      LatticeStaggeredPropagator & quark_propagator_Lsink_Fsrc,
		      LatticeStaggeredPropagator & quark_propagator_Fsink_Fsrc,
		      const multi1d<LatticeColorMatrix> & u , 
		      bool gauge_shift, bool sym_shift ,
		      XMLWriter & xml_out,
		      int j_decay,int t_length, int t_source);

  void compute_8_pions(	multi1d<LatticeStaggeredPropagator> & stag_prop,
			const multi1d<LatticeColorMatrix> & u ,
			bool gauge_shift, bool sym_shift ,
			XMLWriter & xml_out, 
			int j_decay, int t_length, int t_source);

}


#endif

