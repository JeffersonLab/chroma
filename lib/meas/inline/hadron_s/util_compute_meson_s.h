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

  void compute_vary_scalar(LatticeStaggeredPropagator& prop_Lsn_Lsr,
			   LatticeStaggeredPropagator& prop_Fsn_Lsr,
			   LatticeStaggeredPropagator& prop_Lsn_Fsr,
			   LatticeStaggeredPropagator& prop_Fsn_Fsr,
			   const multi1d<LatticeColorMatrix>& u,
			   bool gauge_shift, bool sym_shift,
			   XMLWriter& xml_out, int j_decay, int t_length, int t_source);

  void compute_8_pions(	multi1d<LatticeStaggeredPropagator> & stag_prop,
			const multi1d<LatticeColorMatrix> & u ,
			bool gauge_shift, bool sym_shift ,
			XMLWriter & xml_out, 
			int j_decay, int t_length, int t_source);

  void compute_8_scalars(multi1d<LatticeStaggeredPropagator> & stag_prop,
			 const multi1d<LatticeColorMatrix> & u ,
			 bool gauge_shift, bool sym_shift,
			 XMLWriter & xml_out,
			 int j_decay, int t_length, int t_source);

  void compute_8_vectors(multi1d<LatticeStaggeredPropagator> & stag_prop,
			 const multi1d<LatticeColorMatrix> & u ,
			 bool gauge_shift, bool sym_shift,
			 XMLWriter & xml_out,
			 int j_decay, int t_length, int t_source);
 void compute_8_pions(	multi1d<LatticeStaggeredPropagator> & stag_prop,
			const multi1d<LatticeColorMatrix> & u ,
			bool gauge_shift, bool sym_shift ,
			XMLWriter & xml_out, 
			int j_decay, int t_length, int t_source,
			bool binary_meson_dump, 
			std::string binary_name);

  void compute_8_scalars(multi1d<LatticeStaggeredPropagator> & stag_prop,
			 const multi1d<LatticeColorMatrix> & u ,
			 bool gauge_shift, bool sym_shift,
			 XMLWriter & xml_out,
			 int j_decay, int t_length, int t_source,
			 bool binary_meson_dump, 
			 std::string binary_name);

  void compute_8_vectors(multi1d<LatticeStaggeredPropagator> & stag_prop,
			 const multi1d<LatticeColorMatrix> & u ,
			 bool gauge_shift, bool sym_shift,
			 XMLWriter & xml_out,
			 int j_decay, int t_length, int t_source,
			 bool binary_meson_dump,
			 std::string binary_name);
}


#endif

