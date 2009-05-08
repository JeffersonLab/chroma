/*! \file
 * \brief Inline construction of staggered spectrum
 *
 * Spectrum calculations
 */

#include "meas/inline/hadron_s/inline_spectrum_s.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "meas/smear/ape_smear.h"
#include "meas/hadron/stag_propShift_s.h"

#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "io/param_io.h"
#include "io/enum_io/enum_stochsrc_io.h"

// staggered stuff
#include "handle.h"
#include "state.h"
#include "actions/ferm/fermbcs/simple_fermbc_s.h"
#include "actions/ferm/fermacts/fermacts_s.h"

#include "meas/inline/io/named_objmap.h"

#include "util/ferm/transf.h"
#include "meas/hadron/ks_local_loops.h"
#include "meas/smear/fuzz_smear.h"
#include "meas/hadron/vector_meson_s.h"


#include "meas/inline/make_xml_file.h"

#include "util_compute_quark_prop_s.h"
#include "util_compute_meson_s.h"
#include "meas/hadron/pion_local_s.h"
#include "meas/hadron/g4g5_x_g4g5_local.h"

namespace Chroma 
{ 
  int build_basic_8_props(multi1d<LatticeStaggeredPropagator> stag_prop,
			  stag_src_type type_of_src,
			  bool gauge_shift, bool sym_shift,
			  int fuzz_width,
			  const multi1d<LatticeColorMatrix> & u,
			  const multi1d<LatticeColorMatrix> & u_smr,
			  Handle< SystemSolver<LatticeStaggeredFermion> >
			  & qprop,
			  XMLWriter & xml_out,
			  Real RsdCG, Real Mass, int j_decay ) ;

}


// ------------------------


namespace Chroma { 

  int compute_quark_propagator_s(LatticeStaggeredFermion & psi,
				 stag_src_type type_of_src,
				 bool gauge_shift,
				 bool sym_shift,
				 int fuzz_width,
				 const multi1d<LatticeColorMatrix> & u ,
				 multi1d<LatticeColorMatrix> & u_smr,
				 Handle< SystemSolver<LatticeStaggeredFermion> > & qprop,
				 XMLWriter & xml_out,
				 Real RsdCG, Real Mass, 
				 int j_decay,
				 int src_ind, int color_source);

  int compute_quark_propagator_s(LatticeStaggeredFermion & psi,
				 stag_src_type type_of_src,
				 bool gauge_shift,
				 bool sym_shift,
				 const multi1d<LatticeColorMatrix> & u ,
				 Handle< SystemSolver<LatticeStaggeredFermion> > & qprop,
				 XMLWriter & xml_out,
				 Real RsdCG, Real Mass, 
				 int j_decay,
				 int src_ind, int color_source, int t_source=0);


  int compute_quark_propagator_s(LatticeStaggeredFermion & psi1,
				 LatticeStaggeredFermion & psi2,
				 stag_src_type type_of_src,
				 bool gauge_shift,
				 bool sym_shift,
				 const multi1d<LatticeColorMatrix> & u ,
				 Handle< SystemSolver<LatticeStaggeredFermion> > & qprop1,
				 Handle< SystemSolver<LatticeStaggeredFermion> > & qprop2,
				 XMLWriter & xml_out,
				 Real RsdCG, Real Mass1, Real Mass2, 
				 int j_decay,
				 int src_ind, int color_source, int t_source) ; 

  int compute_quark_propagator_s(LatticeStaggeredFermion & psi,
				 stag_src_type type_of_src,
				 bool gauge_shift,
				 bool sym_shift,
				 int fuzz_width,
				 const multi1d<LatticeColorMatrix> & u ,
				 multi1d<LatticeColorMatrix> & u_smr,
				 Handle< SystemSolver<LatticeStaggeredFermion> > & qprop,
				 XMLWriter & xml_out,
				 Real RsdCG, Real Mass, 
				 int j_decay,
				 int src_ind, int color_source, 
                                 LatticeStaggeredFermion & q_source_in ) ;


  int compute_quark_propagator_s(LatticeStaggeredFermion & psi,
				 stag_src_type type_of_src,
				 bool gauge_shift,
				 bool sym_shift,
				 const multi1d<LatticeColorMatrix> & u ,
				 Handle< SystemSolver<LatticeStaggeredFermion> > & qprop,
				 XMLWriter & xml_out,
				 Real RsdCG, Real Mass, 
				 int j_decay,
				 int src_ind, int color_source, int t_source,
				 LatticeStaggeredFermion & q_source_in ) ; 



  int compute_quark_propagator_s(LatticeStaggeredFermion & psi1,
				 LatticeStaggeredFermion & psi2,
				 stag_src_type type_of_src,
				 bool gauge_shift,
				 bool sym_shift,
				 const multi1d<LatticeColorMatrix> & u ,
				 Handle< SystemSolver<LatticeStaggeredFermion> > & qprop1,
				 Handle< SystemSolver<LatticeStaggeredFermion> > & qprop2,
				 XMLWriter & xml_out,
				 Real RsdCG, Real Mass1, Real Mass2, 
				 int j_decay,
				 int src_ind, int color_source, int t_source, 
                                 LatticeStaggeredFermion & q_source_in) ;


  void ks_compute_baryon(string name,
			 LatticeStaggeredPropagator & quark_propagator,
			 XMLWriter & xml_out,
			 int j_decay, int tlength) ; 

  void ks_compute_baryon(string name,
			 LatticeStaggeredPropagator & quark_propagator_a,
			 LatticeStaggeredPropagator & quark_propagator_b,
			 LatticeStaggeredPropagator & quark_propagator_c,
			 XMLWriter & xml_out,
			 int j_decay, int tlength);

  void ks_compute_baryon(string name,
                         LatticeStaggeredPropagator & quark_propagator_a,
                         LatticeStaggeredPropagator & quark_propagator_b,
                         LatticeStaggeredPropagator & quark_propagator_c,
                         XMLWriter & xml_out,
                         int j_decay, int tlength,
                         bool binary_baryon_dump, std::string binary_name);


  void write_smearing_info(string name, stag_src_type type_of_src,
			   XMLWriter &xml_out, int fuzz_width ) ;


  void compute_vary_baryon_s(XMLWriter &xml_out, int t_source, int fuzz_width,
			     int j_decay, int t_len,
			     LatticeStaggeredPropagator & 
			     quark_propagator_Lsink_Lsrc,
			     LatticeStaggeredPropagator & 
			     quark_propagator_Fsink_Lsrc,
			     LatticeStaggeredPropagator & 
			     quark_propagator_Lsink_Fsrc,
			     LatticeStaggeredPropagator & 
			     quark_propagator_Fsink_Fsrc) ;


  void compute_vary_baryon_s(XMLWriter &xml_out, int t_source, int fuzz_width,
                             int j_decay, int t_len,
                             LatticeStaggeredPropagator &
                             quark_propagator_Lsink_Lsrc,
                             LatticeStaggeredPropagator &
                             quark_propagator_Fsink_Lsrc,
                             LatticeStaggeredPropagator &
                             quark_propagator_Lsink_Fsrc,
                             LatticeStaggeredPropagator &
                             quark_propagator_Fsink_Fsrc,
                             bool binary_baryon_dump,
                             std::string binary_name) ;




  int compute_singlet_ps(LatticeStaggeredFermion & psi,
			 LatticeStaggeredPropagator quark_propagator,
			 stag_src_type type_of_src,
			 bool gauge_shift,
			 bool sym_shift,
			 const multi1d<LatticeColorMatrix> & u ,
			 Handle< SystemSolver<LatticeStaggeredFermion> > &
			 qprop,
			 XMLWriter & xml_out,
			 Real RsdCG, Real Mass, 
			 int j_decay, int t_source, int t_length);

  int compute_vary_singlet_ps(LatticeStaggeredFermion & psi,
       	           LatticeStaggeredPropagator & quark_propagator_Lsink_Lsrc,
	       	   LatticeStaggeredPropagator & quark_propagator_Fsink_Lsrc,
		   LatticeStaggeredPropagator & quark_propagator_Lsink_Fsrc,
		   LatticeStaggeredPropagator & quark_propagator_Fsink_Fsrc,
		   stag_src_type type_of_src,
		   bool gauge_shift,
		   bool sym_shift,
		   const multi1d<LatticeColorMatrix> & u ,
		   const multi1d<LatticeColorMatrix> & u_smr ,
		   Handle< SystemSolver<LatticeStaggeredFermion> > & qprop,
		   XMLWriter & xml_out,
		   Real RsdCG, Real Mass, 
		   int j_decay, int t_source, int t_length, 
			      int fuzz_width);

}

// ----------
/***************************************************************************/
/***************************************************************************/

namespace Chroma { 

/***************************************************************************/

  namespace InlineStaggeredSpectrumEnv 
  { 
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) {
	return new InlineStaggeredSpectrum(InlineStaggeredSpectrumParams(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "SPECTRUM_S";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, 
									  createMeasurement);
	registered = true;
      }
      return success;
    }
  }

/***************************************************************************/


  //! Reader for parameters
  void read(XMLReader& xml, const string& path,
	    InlineStaggeredSpectrumParams::Param_t& param) 
  {
    XMLReader paramtop(xml, path);

    int version;
    read(paramtop, "version", version);

    //    switch (version) 

    read(paramtop, "Meson_local", param.Meson_local);
    read(paramtop, "Meson_charm_local", param.Meson_charm_local);
    read(paramtop, "Meson_charm_noisy_local", param.Meson_charm_noisy_local);
    read(paramtop, "Wilson_loops", param.Wilson_loops);
    read(paramtop, "Pion_nondegen_noisy_local", param.Pion_nondegen_noisy_local);

    param.Pion_nondegen_noisy_local2 = false ;
    if (paramtop.count("Pion_nondegen_noisy_local2") == 1)
    read(paramtop, "Pion_nondegen_noisy_local2", param.Pion_nondegen_noisy_local2);

    param.Pion_nondegen_noisy_local3 = false ;
    if (paramtop.count("Pion_nondegen_noisy_local3") == 1)
    read(paramtop, "Pion_nondegen_noisy_local3", param.Pion_nondegen_noisy_local3);

    param.Pion_nondegen_noisy_local4 = false ;
    if (paramtop.count("Pion_nondegen_noisy_local4") == 1)
    read(paramtop, "Pion_nondegen_noisy_local4", param.Pion_nondegen_noisy_local4);



    read(paramtop, "Baryon_local", param.Baryon_local);
    read(paramtop, "Baryon_vary", param.Baryon_vary);
    read(paramtop, "LocalPion_vary", param.LocalPion_vary);

    if( paramtop.count("LocalScalar_vary") > 0 ) {
      read(paramtop, "LocalScalar_vary", param.LocalScalar_vary);
    }
    else
      {
	param.LocalScalar_vary = false ; 
      }

    read(paramtop, "disconnected_local", param.disconnected_local);
    read(paramtop, "disconnected_fuzz", param.disconnected_fuzz);
    read(paramtop, "singletPs_Conn_local", param.ps4link_singlet_conn);
    read(paramtop, "singletPs_Conn_local_fuzz", param.ps4link_singlet_conn_fuzz);

    read(paramtop, "eight_pions", param.eight_pions);
    read(paramtop, "eight_scalars", param.eight_scalars);
    read(paramtop, "eight_rhos", param.eight_rhos);

    read(paramtop, "t_srce", param.t_srce);
    read(paramtop, "nrow", param.nrow);
    read(paramtop, "sym_shift_oper", param.sym_shift_oper);
    read(paramtop, "gauge_invar_oper", param.gauge_invar_oper);
    read(paramtop, "loop_checkpoint", param.loop_checkpoint);

    param.binary_name = "dump_"  ; 
    param.binary_loop_checkpoint = false ;
    param.binary_meson_dump = false ;
    param.binary_baryon_dump = false ;

    if (paramtop.count("binary_loop_checkpoint") == 1){
      read(paramtop, "binary_loop_checkpoint", 
	   param.binary_loop_checkpoint) ;
    }

    if (paramtop.count("binary_meson_dump") == 1){
      read(paramtop, "binary_meson_dump",
           param.binary_meson_dump) ;
    }

    if (paramtop.count("binary_baryon_dump") == 1){
      read(paramtop, "binary_baryon_dump",
           param.binary_baryon_dump) ;
    }

      if( param.binary_loop_checkpoint || param.binary_baryon_dump  ||
      param.binary_meson_dump )
	{
	  read(paramtop, "binary_name", param.binary_name) ;
	}




    param.fermact = readXMLGroup(paramtop, "FermionAction", "FermAct");

    if( param.Pion_nondegen_noisy_local )
      param.fermact2 = readXMLGroup(paramtop, "FermionAction2", "FermAct");


    if( param.Pion_nondegen_noisy_local2 )
      param.fermact3 = readXMLGroup(paramtop, "FermionAction3", "FermAct");
    if( param.Pion_nondegen_noisy_local )
      param.fermact4 = readXMLGroup(paramtop, "FermionAction4", "FermAct");
    if( param.Pion_nondegen_noisy_local )
      param.fermact5 = readXMLGroup(paramtop, "FermionAction5", "FermAct");


    read(paramtop, "src_seperation", param.src_seperation);

    if( !param.gauge_invar_oper ){
      //read gauge-fixing parameters
      read(paramtop, "GFAccu", param.GFAccu);
      read(paramtop, "OrPara", param.OrPara );
      read(paramtop, "GFMax", param.GFMax );
    }


//    if(( param.disconnected_local )||(param.disconnected_fuzz)){

    if(( param.disconnected_local )||(param.disconnected_fuzz) || param.Meson_charm_noisy_local || param.Pion_nondegen_noisy_local   ){

      read(paramtop, "Number_sample", param.Nsamp);
      read(paramtop, "CFGNO", param.CFGNO);
      // more work
      param.volume_source = GAUSSIAN ;
      read(paramtop, "volume_src", param.volume_source);

    }else{
      param.Nsamp = 10 ;
      param.CFGNO = 10 ;
      param.volume_source = GAUSSIAN ;
    }

    if(((param.Baryon_vary)||(param.ps4link_singlet_conn_fuzz))||
       ((param.disconnected_fuzz)||(param.LocalPion_vary) || param.LocalScalar_vary )){
      read(paramtop, "fuzz_width", param.fuzz_width);
    }else{
      param.fuzz_width = 0 ;
    }
  }

/***************************************************************************/


  //! Writer for parameters
  void write(XMLWriter& xml, const string& path, 
	     const InlineStaggeredSpectrumParams::Param_t& param) {
    push(xml, path);

    int version = 1;
    write(xml, "version", version);

    write(xml, "Meson_local", param.Meson_local);
    write(xml, "Meson_charm_local", param.Meson_charm_local);
    write(xml, "Meson_charm_noisy_local", param.Meson_charm_noisy_local);
    write(xml, "Pion_nondegen_noisy_local", param.Pion_nondegen_noisy_local);

    write(xml, "Pion_nondegen_noisy_local2", param.Pion_nondegen_noisy_local2);
    write(xml, "Pion_nondegen_noisy_local3", param.Pion_nondegen_noisy_local3);
    write(xml, "Pion_nondegen_noisy_local4", param.Pion_nondegen_noisy_local4);

    write(xml, "Wilson_loops", param.Wilson_loops);
    write(xml, "Baryon_local", param.Baryon_local);
    write(xml, "Baryon_vary", param.Baryon_vary);
    write(xml, "disconnected_local", param.disconnected_local);
    write(xml, "disconnected_fuzz", param.disconnected_fuzz);
    write(xml, "nrow", param.nrow);
    write(xml, "t_srce", param.t_srce);
    write(xml, "volume_source",param.volume_source);
    write(xml, "fuzz_width", param.fuzz_width);
    pop(xml);
  }

/***************************************************************************/

  //! Propagator generation params input
  void read(XMLReader& xml, const string& path, 
	    InlineStaggeredSpectrumParams::Quark_Prop_t& input) {
    XMLReader inputtop(xml, path);


    //    read(inputtop, "RsdCG", input.invParam.RsdCG);
    //  read(inputtop, "MaxCG", input.invParam.MaxCG);
    read(inputtop, "invParam", input.invParam); // inverter parameters


    if( inputtop.count("invParam/MaxCGRestart") > 0 ) {
      read(inputtop, "invParam/MaxCGRestart", input.invParam.MaxCGRestart);
    }
    else {
      input.invParam.MaxCGRestart = 0 ; 
    }

  }

/***************************************************************************/


  //! Propagator output
  void write(XMLWriter& xml, const string& path, 
	     const InlineStaggeredSpectrumParams::Quark_Prop_t& input) {
    push(xml, path);
    write(xml, "Mass", input.Mass);
//    write(xml,"Inverter",input.invParam.invType);
    write(xml,"RsdCG", input.invParam.RsdCG);
    write(xml,"MaxCG", input.invParam.MaxCG);
    pop(xml);
  }


/***************************************************************************/
  //! Named object input
  void read(XMLReader& xml, const string& path, InlineStaggeredSpectrumParams::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "gauge_id", input.gauge_id);
  }

/***************************************************************************/

  //! Named object output
  void write(XMLWriter& xml, const string& path, const InlineStaggeredSpectrumParams::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "gauge_id", input.gauge_id);

    pop(xml);
  }

/***************************************************************************/

  // Param stuff
  InlineStaggeredSpectrumParams::InlineStaggeredSpectrumParams()
  { 
    frequency = 0; 
  }

/***************************************************************************/

  InlineStaggeredSpectrumParams::InlineStaggeredSpectrumParams(XMLReader& xml_in, 
							       const std::string& path) {
    try {
      XMLReader paramtop(xml_in, path);

      if (paramtop.count("Frequency") == 1){
	read(paramtop, "Frequency", frequency);
      }else {
	frequency = 1;
      }

      // Parameters for source construction
      read(paramtop, "Param", param);

      // Read in the output propagator/source configuration info
      read(paramtop, "Inversion", prop_param);

      // Read in the gauge field id
      read(paramtop, "NamedObject", named_obj);

      // Possible alternate XML file pattern
      if (paramtop.count("xml_file") != 0) {
	read(paramtop, "xml_file", xml_file);
      }
    } catch(const std::string& e) {
      QDPIO::cerr << "Caught Exception reading XML: " << e << endl;
      QDP_abort(1);
    }
  }

/***************************************************************************/

  void
  InlineStaggeredSpectrumParams::write(XMLWriter& xml_out, 
				       const std::string& path) {
    push(xml_out, path);
    
    Chroma::write(xml_out, "Param", param);
    Chroma::write(xml_out, "Inversion", prop_param);
    Chroma::write(xml_out, "NamedObject", named_obj);
    QDP::write(xml_out, "xml_file", xml_file);

    pop(xml_out);
  }

/***************************************************************************/


  int
  build_basic_8_props(multi1d<LatticeStaggeredPropagator> &stag_prop,
		      stag_src_type type_of_src, 
		      bool gauge_shift, bool sym_shift,
		      int fuzz_width,
		      const multi1d<LatticeColorMatrix> & u,
		      multi1d<LatticeColorMatrix> & u_smr,
		      Handle< SystemSolver<LatticeStaggeredFermion> > & 
		      qprop,
		      XMLWriter & xml_out,
		      Real RsdCG, Real Mass, int j_decay){


    LatticeStaggeredFermion q_source, psi ;
    LatticeStaggeredFermion q_source_fuzz ;
    LatticeStaggeredPropagator quark_propagator;

    int ncg_had = 0;                         // count CG iterations

    psi = zero;                              // note this is ``zero'' and not 0



    for(int src_ind = 0; src_ind < 8 ; ++src_ind){
      stag_prop[src_ind] = zero ;
    
      for(int color_source = 0; color_source < Nc; ++color_source) {

	ncg_had += compute_quark_propagator_s(psi,type_of_src,
					      fuzz_width, 
					      gauge_shift, sym_shift,
					      u, u_smr, qprop, xml_out,
					      RsdCG, Mass, j_decay,
					      src_ind, color_source) ;


	/*
         * Move the solution to the appropriate components
         * of quark propagator.
         */
	FermToProp(psi, quark_propagator, color_source);
      }  //color_source
    
      stag_prop[src_ind] = quark_propagator;

    }// end src_ind

    return ncg_had;

  }
/***************************************************************************/

  // no fuzzing version
  int
  build_basic_8_props(multi1d<LatticeStaggeredPropagator> &stag_prop,
		      stag_src_type type_of_src, 
		      bool gauge_shift, bool sym_shift,
		      const multi1d<LatticeColorMatrix> & u,
		      Handle< SystemSolver<LatticeStaggeredFermion> > & 
		      qprop,
		      XMLWriter & xml_out,
		      Real RsdCG, Real Mass, int j_decay){


    LatticeStaggeredFermion q_source, psi ;
    LatticeStaggeredFermion q_source_fuzz ;
    LatticeStaggeredPropagator quark_propagator;

    int ncg_had = 0;                         // count CG iterations

    psi = zero;                              // note this is ``zero'' and not 0



    for(int src_ind = 0; src_ind < 8 ; ++src_ind){
      stag_prop[src_ind] = zero ;
    
      for(int color_source = 0; color_source < Nc; ++color_source) {

	ncg_had += compute_quark_propagator_s(psi,type_of_src,
					      gauge_shift, sym_shift,
					      u, qprop, xml_out,
					      RsdCG, Mass, j_decay,
					      src_ind, color_source) ;


	/*
         * Move the solution to the appropriate components
         * of quark propagator.
         */
	FermToProp(psi, quark_propagator, color_source);
      }  //color_source
    
      stag_prop[src_ind] = quark_propagator;

    }// end src_ind

    return ncg_had;
  }

/***************************************************************************/
  int
  MakeFuzzedCornerProp(LatticeStaggeredFermion & psi,
		       int fuzz_width,
		       bool gauge_shift, bool sym_shift,
		       const multi1d<LatticeColorMatrix> & u , 
		       multi1d<LatticeColorMatrix> & u_smr,
		       Handle< SystemSolver<LatticeStaggeredFermion> > & 
		       qprop,
		       XMLWriter & xml_out,
		       Real RsdCG, Real Mass,
		       int j_decay,
		       bool do_fuzzing,
		       LatticeStaggeredFermion &psi_fuzz,
		       LatticeStaggeredPropagator &quark_propagator_Lsink_Lsrc,
		       LatticeStaggeredPropagator &quark_propagator_Fsink_Lsrc,
		       LatticeStaggeredPropagator &quark_propagator_Lsink_Fsrc,
		       LatticeStaggeredPropagator &quark_propagator_Fsink_Fsrc
    ){

    //    stag_src_type type_of_src = LOCAL_SRC ;
    stag_src_type type_of_src = GAUGE_INVAR_LOCAL_SOURCE;

    QDPIO::cout << "LOCAL INVERSIONS"  << endl;
    int ncg_had = 0 ;

    for(int color_source = 0; color_source < Nc; ++color_source){
      psi = zero;    // note this is ``zero'' and not 0

      const int src_ind = 0 ;
      ncg_had += compute_quark_propagator_s(psi,type_of_src, 
					    gauge_shift, sym_shift,
					    fuzz_width,
					    u, u_smr, qprop, xml_out,
					    RsdCG, Mass,
					    j_decay, src_ind, color_source) ;

      /*
       * Move the solution to the appropriate components
       * of quark propagator.
       */
      FermToProp(psi, quark_propagator_Lsink_Lsrc, color_source);


      //  fuzz at the sink 
      if( do_fuzzing ) {
	fuzz_smear(u_smr, psi, psi_fuzz, fuzz_width, j_decay) ;
	FermToProp(psi_fuzz, quark_propagator_Fsink_Lsrc, color_source);
      }
    } // color_source


    if( do_fuzzing ){
      type_of_src = FUZZED_SRC ;

      QDPIO::cout << "FUZZED SOURCE INVERSIONS"  << endl;

      for(int color_source = 0; color_source < Nc; ++color_source) {
	psi = zero;            // note this is ``zero'' and not 0

	const int src_ind = 0 ;

	ncg_had += compute_quark_propagator_s(psi,type_of_src, 
					      gauge_shift, sym_shift, 
					      fuzz_width,
					      u, u_smr, qprop, xml_out,
					      RsdCG, Mass, j_decay, src_ind, 
					      color_source);

	/*
	 * Move the solution to the appropriate components
         * of quark propagator.
         */
	FermToProp(psi, quark_propagator_Lsink_Fsrc, color_source);


        //  fuzz at the sink 
	fuzz_smear(u_smr, psi, psi_fuzz, fuzz_width, j_decay) ;
	FermToProp(psi_fuzz, quark_propagator_Fsink_Fsrc, color_source);

      }  //color_source

    }  // end of compute fuzzed correlator

    return ncg_had;
  }
/***************************************************************************/
   int
  MakeCornerProp(LatticeStaggeredFermion & psi,
		 bool gauge_shift, bool sym_shift,
		 const multi1d<LatticeColorMatrix> & u ,
		 Handle< SystemSolver<LatticeStaggeredFermion> > &  qprop,
		 XMLWriter & xml_out,
		 Real RsdCG, Real Mass,
		 int j_decay,
		 LatticeStaggeredPropagator &quark_propagator_Lsink_Lsrc,
		 stag_src_type type_of_src, int t_source = 0  ){

    //    stag_src_type type_of_src = LOCAL_SRC ;
    //    stag_src_type type_of_src = GAUGE_INVAR_LOCAL_SOURCE;

     //    QDPIO::cout << "LOCAL INVERSIONS"  << endl;
    int ncg_had = 0 ;

    for(int color_source = 0; color_source < Nc; ++color_source){
      psi = zero;    // note this is ``zero'' and not 0

      const int src_ind = 0 ;
      ncg_had += compute_quark_propagator_s(psi,type_of_src, 
					    gauge_shift, sym_shift,
					    u, qprop, xml_out,
					    RsdCG, Mass,
					    j_decay, src_ind, 
					    color_source,t_source) ;

      /*
       * Move the solution to the appropriate components
       * of quark propagator.
       */
      FermToProp(psi, quark_propagator_Lsink_Lsrc, color_source);

    } // color_source

    return ncg_had;
  }



   int
  MakeCornerProp(LatticeStaggeredFermion & psi,
		 bool gauge_shift, bool sym_shift,
		 const multi1d<LatticeColorMatrix> & u ,
		 Handle< SystemSolver<LatticeStaggeredFermion> > &  qprop,
		 XMLWriter & xml_out,
		 Real RsdCG, Real Mass,
		 int j_decay,
		 LatticeStaggeredPropagator &quark_propagator_Lsink_Lsrc,
		 stag_src_type type_of_src, int t_source,
		 LatticeStaggeredPropagator &qsource_out){

    //    stag_src_type type_of_src = LOCAL_SRC ;
    //    stag_src_type type_of_src = GAUGE_INVAR_LOCAL_SOURCE;
       LatticeStaggeredFermion q_source ;


     //    QDPIO::cout << "LOCAL INVERSIONS"  << endl;
    int ncg_had = 0 ;

    for(int color_source = 0; color_source < Nc; ++color_source){
      psi = zero;    // note this is ``zero'' and not 0

      // load the source
      PropToFerm(qsource_out,q_source,  color_source);

      const int src_ind = 0 ;
      ncg_had += compute_quark_propagator_s(psi,type_of_src, 
					    gauge_shift, sym_shift,
					    u, qprop, xml_out,
					    RsdCG, Mass,
					    j_decay, src_ind, 
					    color_source,t_source,q_source) ;

      /*
       * Move the solution to the appropriate components
       * of quark propagator.
       */
      FermToProp(psi, quark_propagator_Lsink_Lsrc, color_source);



    } // color_source

    return ncg_had;
  }





   int
  MakeCornerProp(LatticeStaggeredFermion & psi_1,
		 bool gauge_shift, bool sym_shift,
		 const multi1d<LatticeColorMatrix> & u ,
		 Handle< SystemSolver<LatticeStaggeredFermion> > &  qprop_1,
		 Handle< SystemSolver<LatticeStaggeredFermion> > &  qprop_2,
		 XMLWriter & xml_out,
		 Real RsdCG, Real Mass_1,Real Mass_2,
		 int j_decay,
		 LatticeStaggeredPropagator &quark_propagator_Lsink_Lsrc_1,
		 LatticeStaggeredPropagator &quark_propagator_Lsink_Lsrc_2,
		 stag_src_type type_of_src, int t_source = 0  ){

     LatticeStaggeredFermion psi_2 ; 
    int ncg_had = 0 ;

    for(int color_source = 0; color_source < Nc; ++color_source){
      psi_1 = zero;    // note this is ``zero'' and not 0
      psi_2 = zero;    // note this is ``zero'' and not 0

      const int src_ind = 0 ;
      ncg_had += compute_quark_propagator_s(psi_1,psi_2,type_of_src, 
					    gauge_shift, sym_shift,
					    u, qprop_1,qprop_2, xml_out,
					    RsdCG, Mass_1,Mass_2,
					    j_decay, src_ind, 
					    color_source,t_source) ;

      /*
       * Move the solution to the appropriate components
       * of quark propagator.
       */
      FermToProp(psi_1, quark_propagator_Lsink_Lsrc_1, color_source);
      FermToProp(psi_2, quark_propagator_Lsink_Lsrc_2, color_source);

    } // color_source

    return ncg_had;
  }


/*  Also store the random source for additional inversions


*/

   int
  MakeCornerProp(LatticeStaggeredFermion & psi_1,
		 bool gauge_shift, bool sym_shift,
		 const multi1d<LatticeColorMatrix> & u ,
		 Handle< SystemSolver<LatticeStaggeredFermion> > &  qprop_1,
		 Handle< SystemSolver<LatticeStaggeredFermion> > &  qprop_2,
		 XMLWriter & xml_out,
		 Real RsdCG, Real Mass_1,Real Mass_2,
		 int j_decay,
		 LatticeStaggeredPropagator &quark_propagator_Lsink_Lsrc_1,
		 LatticeStaggeredPropagator &quark_propagator_Lsink_Lsrc_2,
		 stag_src_type type_of_src, int t_source, 
		 LatticeStaggeredPropagator &qsource_out  ){

     LatticeStaggeredFermion psi_2 ; 
     LatticeStaggeredFermion q_source ; 


    int ncg_had = 0 ;

    for(int color_source = 0; color_source < Nc; ++color_source){
      psi_1 = zero;    // note this is ``zero'' and not 0
      psi_2 = zero;    // note this is ``zero'' and not 0

      // load the source 
      PropToFerm( qsource_out, q_source,color_source);


      const int src_ind = 0 ;
      ncg_had += compute_quark_propagator_s(psi_1,psi_2,type_of_src, 
					    gauge_shift, sym_shift,
					    u, qprop_1,qprop_2, xml_out,
					    RsdCG, Mass_1,Mass_2,
					    j_decay, src_ind, 
					    color_source,t_source,q_source) ;

      /*
       * Move the solution to the appropriate components
       * of quark propagator.
       */
      FermToProp(psi_1, quark_propagator_Lsink_Lsrc_1, color_source);
      FermToProp(psi_2, quark_propagator_Lsink_Lsrc_2, color_source);


      // save the source 
      FermToProp(q_source, qsource_out, color_source);

    } // color_source

    return ncg_had;
  }



/***************************************************************************/

  void
  DoFuzzing(const multi1d<LatticeColorMatrix> & u,
	    multi1d<LatticeColorMatrix> & u_smr,
	    int j_decay){

    QDPIO::cout << "Starting to APE smear the gauge configuration" << endl;
  
    Real sm_fact = 2.5;   // typical parameter
    int sm_numb = 10;     // number of smearing hits
	
    int BlkMax = 100;   // max iterations in max-ing trace
    Real BlkAccu = 1.0e-5;  // accuracy of max-ing
	
    u_smr = u;
    for(int i=0; i < sm_numb; ++i){
      multi1d<LatticeColorMatrix> u_tmp(Nd);
	    
      for(int mu = 0; mu < Nd; ++mu){
	if ( mu != j_decay ){
	  APE_Smear(u_smr, u_tmp[mu], mu, 0, sm_fact, BlkAccu, BlkMax, 
		    j_decay);
	}else{
	  u_tmp[mu] = u_smr[mu];
	}
      }
      u_smr = u_tmp;
    }



  }


/***************************************************************************/

  // Function call
  void 
  InlineStaggeredSpectrum::operator()(unsigned long update_no,
				      XMLWriter& xml_out) {
    // If xml file not empty, then use alternate
    if (params.xml_file != ""){
      string xml_file = makeXMLFileName(params.xml_file, update_no);

      push(xml_out, "spectrum_s");
      write(xml_out, "update_no", update_no);
      write(xml_out, "xml_file", xml_file);
      pop(xml_out);

      XMLFileWriter xml(xml_file);
      func(update_no, xml);
    }
    else
    {
      func(update_no, xml_out);
    }
  }


  /***  **/

  void meson_charm(LatticeStaggeredPropagator & quark_prop, 
		   XMLWriter& xml_out, 
		   const multi1d<LatticeColorMatrix> & u,
		   int  t_source, int j_decay, int t_length)
  {
  push(xml_out, "local_meson_charm_correlators");

  // local pseudoscalar pion
  staggered_local_pion pion(t_length,u) ;
  pion.compute(quark_prop,quark_prop,j_decay) ;
  pion.dump(t_source,xml_out) ;

  bool avg_equiv_mom = true ;
  int mom2_max = 2 ; 
  // non-zero momentum
  multi1d<int> tsrc(4) ; 
  tsrc[0] =  tsrc[1] =  tsrc[2] =  tsrc[3] =  0 ; 
  tsrc[j_decay] = t_source ;

  SftMom phases(mom2_max, tsrc, avg_equiv_mom,j_decay);


  pion.compute_and_dump(quark_prop,quark_prop,j_decay,t_source,
	       phases,xml_out) ;

  // ( gamma_4 gamma_5 cross gamma_4 gamma_5 ) pion
  g4g5_x_g4g5_local_meson pion_g4g5(t_length,u)  ;  
  pion_g4g5.compute(quark_prop,quark_prop,j_decay) ;
  pion_g4g5.dump(t_source,xml_out) ;

  // vector mesons
  vector_meson vector(t_length,  u) ;
  vector.compute(quark_prop,j_decay);
  vector.dump(t_source,xml_out);

  pop(xml_out);

}


  /**
    Non-degenerate noisy heavy-light pseudoscalars

  **/

  void noisy_pion_nondegen(LatticeStaggeredPropagator & quark_prop_1, 
			   Real Mass1,
			   LatticeStaggeredPropagator & quark_prop_2, 
			   Real Mass2,
			   XMLWriter& xml_out, 
			   const multi1d<LatticeColorMatrix> & u,
			   int  t_source, int j_decay, int t_length)
  {
  push(xml_out, "non_degenerate_noisy_pseudoscalar_correlators");

  // local pseudoscalar pion
  push(xml_out, "meson_11");
  write(xml_out,"Mass",Mass1); 
  staggered_local_pion pion_1(t_length,u) ;
  pion_1.compute(quark_prop_1,quark_prop_1,j_decay) ;
  pion_1.dump(t_source,xml_out) ;
  pop(xml_out);

  push(xml_out, "meson_22");
  write(xml_out,"Mass",Mass2); 
  staggered_local_pion pion_2(t_length,u) ;
  pion_2.compute(quark_prop_2,quark_prop_2,j_decay) ;
  pion_2.dump(t_source,xml_out) ;
  pop(xml_out);

  push(xml_out, "meson_12");
  write(xml_out,"Mass1",Mass1); 
  write(xml_out,"Mass2",Mass2); 
  staggered_local_pion pion_12(t_length,u) ;
  pion_12.compute(quark_prop_1,quark_prop_2,j_decay) ;
  pion_12.dump(t_source,xml_out) ;
  pop(xml_out);


  //-----------------------
  pop(xml_out);

}






/***************************************************************************/

  // Real work done here
  void 
  InlineStaggeredSpectrum::func(unsigned long update_no,
				XMLWriter& xml_out) 
  {
    START_CODE();

    StopWatch snoop;
    snoop.reset();
    snoop.start();

    // Test and grab a reference to the gauge field
    XMLBufferWriter gauge_xml;
    try
    {
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);
      TheNamedObjMap::Instance().get(params.named_obj.gauge_id).getRecordXML(gauge_xml);
    }
    catch( std::bad_cast ) 
    {
      QDPIO::cerr << InlineStaggeredSpectrumEnv::name << ": caught dynamic cast error" 
		  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << InlineStaggeredSpectrumEnv::name << ": map call failed: " << e 
		  << endl;
      QDP_abort(1);
    }
    const multi1d<LatticeColorMatrix>& u = 
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

    QDPIO::cout << InlineStaggeredSpectrumEnv::name << ": Spectroscopy for Staggered-like fermions" 
		<< endl;
    QDPIO::cout << "Gauge group: SU(" << Nc << ")" << endl;

    push(xml_out, "spectrum_s");
    write(xml_out, "update_no", update_no);

    /*  set some flags **/

    bool do_Baryon_local         = params.param.Baryon_local;
    bool do_ps4_singlet          = params.param.ps4link_singlet_conn;
    bool do_ps4_singlet_fuzz     = params.param.ps4link_singlet_conn_fuzz;
    bool Meson_local             = params.param.Meson_local ;
    bool Meson_charm_local       = params.param.Meson_charm_local ;
    bool Meson_charm_noisy_local = params.param.Meson_charm_noisy_local ;
    bool Wilson_loops             = params.param.Wilson_loops ;
    bool Pion_nondegen_noisy_local = params.param.Pion_nondegen_noisy_local ;

    bool Pion_nondegen_noisy_local2 = params.param.Pion_nondegen_noisy_local2 ;
    bool Pion_nondegen_noisy_local3 = params.param.Pion_nondegen_noisy_local3 ;
    bool Pion_nondegen_noisy_local4 = params.param.Pion_nondegen_noisy_local4 ;

    bool do_Baryon_vary          = params.param.Baryon_vary ;
    bool do_LocalPion_vary       = params.param.LocalPion_vary;
    bool do_LocalScalar_vary     = params.param.LocalScalar_vary;
    bool do_8_pions              = params.param.eight_pions;
    bool do_8_scalars            = params.param.eight_scalars;
    bool do_8_rhos               = params.param.eight_rhos;
    bool do_fuzzed_disc_loops    = params.param.disconnected_fuzz  ;
    bool do_local_disc_loops     = params.param.disconnected_local  ;

    //change this to do stochastic connected loops
    bool do_stoch_conn_corr       = false;

    bool do_fuzzing               = false;
    bool do_variational_spectra   = false;

    bool need_basic_8             = false;  

    bool need_fuzzed_corner_prop  = false;

    bool done_ps4_singlet         = false;
    bool done_ps4_singlet_fuzz    = false;
    bool done_local_baryons       = false;
    bool done_fuzzed_baryons      = false;
    bool done_local_disc_loops    = false;
    bool done_fuzzed_disc_loops   = false;
    bool done_meson_corr          = false;
    bool done_meson_charm_corr    = false;

    //shouldnt be hard-coded
    stag_src_type type_of_src = GAUGE_INVAR_LOCAL_SOURCE ;


    if( ( do_Baryon_vary || do_LocalPion_vary  ||  do_LocalScalar_vary )|| 
	(do_fuzzed_disc_loops||do_ps4_singlet_fuzz)){
      // need smeared links
      do_fuzzing = true ;
    }

    if (( do_Baryon_vary || do_LocalPion_vary)||(do_ps4_singlet_fuzz) || do_LocalScalar_vary  ) {

      // make the fuzzed corner props
      // (LsrcLsnk, LsrcFsnk,FsrcLsink,FsrcFsnk)

      do_variational_spectra = true;
    }

    if( do_8_pions || do_8_scalars || do_8_rhos){
      need_basic_8 = true;
    }


    int          ncg_had = 0;                 // tally of all CG iterations


    bool         gauge_shift     = true;
    bool         sym_shift       = true;
    bool         loop_checkpoint = false;
    const int    j_decay         = Nd-1;
    int          Nsamp;
    int          CFGNO;
    VolSrc_type  volume_source;
    int          src_seperation;
    int          t_length;


    gauge_shift        = params.param.gauge_invar_oper;   // use covar shift?
    sym_shift          = params.param.sym_shift_oper;     // use symm shift?
    loop_checkpoint    = params.param.loop_checkpoint;    // write all meas?
    Nsamp              = params.param.Nsamp;              // # of stoch srcs
    CFGNO              = params.param.CFGNO;              // seed
    volume_source      = params.param.volume_source;      // type of vol src
    src_seperation     = params.param.src_seperation;     // sep of dilute srcs
    t_length           = params.param.nrow[j_decay];      // length of t dir
    Real RsdCG         = params.prop_param.invParam.RsdCG;// CG residual
    //    Real Mass          = params.prop_param.Mass;          // fermion mass
    Real Mass ; 
    int  fuzz_width    = params.param.fuzz_width;         // fuzzing width
    int  t_source      = params.param.t_srce[j_decay];    // source t coord

    /*
     * Sanity checks
     */


    QDPIO::cout << "Volume: " << params.param.nrow[0];
    for (int i=1; i<Nd; ++i) {
      QDPIO::cout << " x " << params.param.nrow[i];
    }
    QDPIO::cout << endl;

    proginfo(xml_out);    // Print out basic program info


    // Write out the input
    params.write(xml_out, "Input");

    // Write out the config info
    write(xml_out, "Config_info", gauge_xml);

    push(xml_out, "Output_version");
    write(xml_out, "out_version", 1);
    pop(xml_out);

    // First calculate some gauge invariant observables just for info.
    MesPlq(xml_out, "Observables", u);

    // Typedefs to save typing
    typedef LatticeStaggeredFermion      T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

#if 0 
    // Create a fermion state
    // boundary has been deleted
    Handle< CreateFermState<T,P,Q> > cfs(new CreateSimpleFermState<T,P,Q>(params.param.boundary));

    // Initialize fermion action
    AsqtadFermActParams asq_param;
    asq_param.Mass = params.prop_param.Mass;
    asq_param.u0   = params.prop_param.u0;
    AsqtadFermAct S_f(cfs, asq_param);
    Handle< FermState<T,P,Q> > state(S_f.createState(u));
#endif

    //
    // Initialize fermion action
    //

    XMLReader fermact_reader ;
    // Make a memory 'input stream' out of the XML, so we can open an
    // XML Reader on it.
    try{
      std::istringstream is(params.param.fermact.xml);

      // Open a reader on the memory stream.
      //  XMLReader fermact_reader(is);
      fermact_reader.open(is);
    }
    catch (...)
      {
	QDPIO::cerr << "Error reading action name " << endl;
	throw;
      }

    Handle< StaggeredTypeFermAct< T,P,Q> > fermact(
TheStagTypeFermActFactory::Instance().createObject(params.param.fermact.id, 
fermact_reader, 
params.param.fermact.path));
    // Cast of a pointer to a reference?
    StaggeredTypeFermAct<T,P,Q>& S_f= *(fermact);

    Mass = S_f.getQuarkMass() ;
    params.prop_param.Mass = S_f.getQuarkMass() ;

    // Set up a state for the current u,
    // (compute fat & triple links)
    // Use S_f.createState so that S_f can pass in u0

    Handle< FermState<T,P,Q> > state(S_f.createState(u));

    // Jiggery-pokery to turn a CG struct into a GroupXML_t for the qprops
    GroupXML_t inv_param;
    {
      XMLBufferWriter xml_buf;
      write(xml_buf, "InvertParam", params.prop_param.invParam);
      XMLReader xml_in(xml_buf);
      inv_param = readXMLGroup(xml_in, "/InvertParam", "invType");
    }
    Handle< SystemSolver<LatticeStaggeredFermion> > 
      qprop(S_f.qprop(state, inv_param));



    //
    //  local inversions
    // 
    LatticeStaggeredFermion psi ;
    LatticeStaggeredFermion q_source ; 


    if( Wilson_loops ) 
      {
	//	push(xml_out, "Wilson_loops") ;
        Wloop(xml_out,"Wilson_loops", u) ;
	//  pop(xml_out);
      }



    // Do local disconnected loops if we are not going to do fuzzed 
    //   disconnected loops.

    if (((  do_local_disc_loops ) && (do_stoch_conn_corr )) && 
    	(!do_fuzzed_disc_loops)){
      push(xml_out, "disconnected_loops");

      StopWatch swatch;
      swatch.start();
      ks_local_loops_and_stoch_conn(qprop, q_source, psi , u, xml_out, 
				    gauge_shift, sym_shift, loop_checkpoint,
				    t_length, Mass, Nsamp, RsdCG, CFGNO, 
				    volume_source, src_seperation, j_decay);
      swatch.stop();
      double time_in_sec  = swatch.getTimeInSeconds();
      QDPIO::cout << "PROF1:ks_local_loops_and_stoch_conn" << time_in_sec << " sec" << endl;

      pop(xml_out);

      done_local_disc_loops = true;
    }


    if((( do_local_disc_loops  )  && (!done_local_disc_loops))&&
       (!do_fuzzed_disc_loops)){
      push(xml_out, "disconnected_loops");

      StopWatch swatch;
      swatch.start();
      ks_local_loops(qprop, q_source, psi , u, xml_out, 
		     gauge_shift, sym_shift, loop_checkpoint,
		     t_length, Mass, Nsamp, RsdCG, CFGNO, volume_source, 
		     src_seperation, j_decay);

      swatch.stop();
      double time_in_sec  = swatch.getTimeInSeconds();
      QDPIO::cout << "PROF2:ks_local_loops " << time_in_sec << " sec" << endl;

      done_local_disc_loops = true;

      pop(xml_out);
    }



    if (need_basic_8 ){
      //Dont need to allocate u_smr here

      multi1d<LatticeStaggeredPropagator> stag_prop(8);

      StopWatch swatch;
      swatch.start();
      ncg_had += build_basic_8_props(stag_prop, type_of_src, gauge_shift,
				     sym_shift, u, qprop, xml_out, RsdCG,  
				     Mass, j_decay);
      swatch.stop();
      double time_in_sec  = swatch.getTimeInSeconds();
      QDPIO::cout << "PROF3:build_basic_8_props " << time_in_sec << " sec" << endl;


      // put the spectrum calls here 
      if(do_8_pions){
	// do pion stuff

	StopWatch swatch;
	swatch.start();
	compute_8_pions( stag_prop, u , gauge_shift, sym_shift,
			 xml_out, j_decay, t_length, t_source,
			 params.param.binary_meson_dump,
			 params.param.binary_name);

	swatch.stop();
	double time_in_sec  = swatch.getTimeInSeconds();
	QDPIO::cout << "PROF4:compute_8_pions " << time_in_sec << " sec" << endl;

      }
      if(do_8_scalars){
	// do scalar stuff
	compute_8_scalars( stag_prop, u,  gauge_shift, sym_shift,
			   xml_out, j_decay, t_length, t_source,
			   params.param.binary_meson_dump,
			   params.param.binary_name);
	
      }
      if(do_8_rhos){
	// do vector stuff
	compute_8_vectors( stag_prop, u,  gauge_shift, sym_shift,
			   xml_out, j_decay, t_length, t_source,
			   params.param.binary_meson_dump,
			   params.param.binary_name);
      }


      // if we need to do 4-link singlets and local baryons, do them here
      // so we can re-use the stag_pro[0] as the local corner prop

      if((( do_ps4_singlet ) && (!done_ps4_singlet))&&(!do_ps4_singlet_fuzz)) {
	//only do the non-fuzzed ps singlet if we are not doing the fuzzed 
	// ps singlet connected, as the latter does both

	StopWatch swatch;
	swatch.start();

	type_of_src = GAUGE_INVAR_LOCAL_SOURCE ;
	ncg_had += 
	  compute_singlet_ps(psi,stag_prop[0], type_of_src, gauge_shift, 
			     sym_shift, u, qprop, xml_out, RsdCG,Mass,
			     j_decay, t_source, t_length);

	done_ps4_singlet = true;
	swatch.stop();
	double time_in_sec  = swatch.getTimeInSeconds();
	QDPIO::cout << "PROF5:compute_singlet_ps " << time_in_sec << " sec" << endl;


      }


      if((( do_Baryon_local) &&(!done_local_baryons)) && (!do_Baryon_vary)) {

	push(xml_out, "baryon_correlators");

	// describe the source
	string NN ;
	write(xml_out, "source_time", t_source);
	push(xml_out, "smearing_info");
	NN = "L" ;
	write_smearing_info(NN, LOCAL_SRC,xml_out, fuzz_width) ;
	pop(xml_out);
	string b_tag("srcLLL_sinkLLL_nucleon") ;
	string b_filename(params.param.binary_name+b_tag);
	ks_compute_baryon(b_tag,stag_prop[0],stag_prop[0],stag_prop[0],
			  xml_out, j_decay, t_length,
			  params.param.binary_baryon_dump,b_filename) ;

	pop(xml_out);
	done_local_baryons = true;
      }

 
    } // basic 8 



    // fuzz the gauge configuration. if needed

    if( do_fuzzing  ){

      LatticeStaggeredFermion psi_fuzz ;
      multi1d<LatticeColorMatrix> u_smr(Nd);

      u_smr = u ;

      {
	StopWatch swatch;
	swatch.start();
	DoFuzzing(u, u_smr, j_decay);
	swatch.stop();
	double time_in_sec  = swatch.getTimeInSeconds();
	QDPIO::cout << "PROF6: DoFuzzing " << time_in_sec << " sec" << endl;
      }

      if ( do_variational_spectra ){

	// only allocate these if needed
	LatticeStaggeredPropagator quark_propagator_Lsink_Lsrc;
	LatticeStaggeredPropagator quark_propagator_Fsink_Lsrc ;
	LatticeStaggeredPropagator quark_propagator_Lsink_Fsrc ;
	LatticeStaggeredPropagator quark_propagator_Fsink_Fsrc ;

	{
	  StopWatch swatch;
	  swatch.start();

	  ncg_had+=  MakeFuzzedCornerProp(psi, fuzz_width, gauge_shift, 
					sym_shift, u, u_smr, qprop, 
					xml_out, RsdCG, Mass, j_decay,
					do_fuzzing, psi_fuzz,
					quark_propagator_Lsink_Lsrc,
					quark_propagator_Fsink_Lsrc,
					quark_propagator_Lsink_Fsrc,
					quark_propagator_Fsink_Fsrc );
	  swatch.stop();
	  double time_in_sec  = swatch.getTimeInSeconds();
	  QDPIO::cout << "PROF7:MakeFuzzedCornerProp " << time_in_sec << " sec" << endl;

	}


	  {

	    //testing crap
	    DComplex productLL;
	    DComplex productLF;
	    DComplex productFL;
	    DComplex productFF;

	    LatticeComplex test_corr_fnLL;
	    LatticeComplex test_corr_fnLF;
	    LatticeComplex test_corr_fnFL;
	    LatticeComplex test_corr_fnFF;

	    test_corr_fnLL = trace(adj(quark_propagator_Lsink_Lsrc)
				 *quark_propagator_Lsink_Lsrc);
	    test_corr_fnFL = trace(adj(quark_propagator_Fsink_Lsrc)
				 *quark_propagator_Fsink_Lsrc);
	    test_corr_fnLF = trace(adj(quark_propagator_Lsink_Fsrc)
				 *quark_propagator_Lsink_Fsrc);
	    test_corr_fnFF = trace(adj(quark_propagator_Fsink_Fsrc)
				 *quark_propagator_Fsink_Fsrc);

	    productLL = sum(test_corr_fnLL);
	    productFL = sum(test_corr_fnFL);
	    productLF = sum(test_corr_fnLF);
	    productFF = sum(test_corr_fnFF);

	    QDPIO::cout << "product of LL= " << productLL << std::endl; 
	    QDPIO::cout << "product of FL= " << productFL << std::endl; 
	    QDPIO::cout << "product of LF= " << productLF << std::endl; 
	    QDPIO::cout << "product of FF= " << productFF << std::endl; 
	    // one copy on output
	  }




	if( do_Baryon_vary ) {

	  compute_vary_baryon_s(xml_out,t_source,
				fuzz_width,
				j_decay,t_length,
				quark_propagator_Lsink_Lsrc,
				quark_propagator_Fsink_Lsrc,
				quark_propagator_Lsink_Fsrc,
				quark_propagator_Fsink_Fsrc,
				params.param.binary_baryon_dump,
				params.param.binary_name) ;

	  done_local_baryons = true;
	  done_fuzzed_baryons = true;

	}


	if( do_LocalPion_vary ) {
	  compute_vary_ps(quark_propagator_Lsink_Lsrc,
			  quark_propagator_Fsink_Lsrc,
			  quark_propagator_Lsink_Fsrc,
			  quark_propagator_Fsink_Fsrc,
			  u, gauge_shift, sym_shift,
			  xml_out,j_decay,
			  t_length,t_source);
	}


	if( do_LocalScalar_vary ) {
	  compute_vary_scalar(quark_propagator_Lsink_Lsrc,
			  quark_propagator_Fsink_Lsrc,
			  quark_propagator_Lsink_Fsrc,
			  quark_propagator_Fsink_Fsrc,
			  u, gauge_shift, sym_shift,
			  xml_out,j_decay,
			  t_length,t_source);
	}

	if(( do_ps4_singlet_fuzz ) && (!done_ps4_singlet_fuzz)) {


	  ncg_had += 
	    compute_vary_singlet_ps(psi, quark_propagator_Lsink_Lsrc,
				    quark_propagator_Fsink_Lsrc,
				    quark_propagator_Lsink_Fsrc,
				    quark_propagator_Fsink_Fsrc,
				    type_of_src, gauge_shift, sym_shift, u ,
				    u_smr , qprop, xml_out, RsdCG, Mass, 
				    j_decay, t_source, t_length,fuzz_width);

	  done_ps4_singlet = true;
	  done_ps4_singlet_fuzz = true;

	}


	// if we need to do local 4-link singlets and local baryons, 
	// do them here so we can re-use quark_propagator_Lsink_Lsrc 
	// as the local corner prop

	if(( do_ps4_singlet ) && (!done_ps4_singlet)) {
	  type_of_src = GAUGE_INVAR_LOCAL_SOURCE ;

	  ncg_had += 
	    compute_singlet_ps(psi,quark_propagator_Lsink_Lsrc, type_of_src,
			       gauge_shift, sym_shift, u, qprop, xml_out, 
			       RsdCG, Mass, j_decay, t_source, t_length);

	  done_ps4_singlet = true;
	}


	if(( do_Baryon_local) &&(!done_local_baryons)) {

	  push(xml_out, "baryon_correlators");

	  // describe the source
	  string NN ;
	  write(xml_out, "source_time", t_source);
	  push(xml_out, "smearing_info");
	  NN = "L" ;
	  write_smearing_info(NN, LOCAL_SRC,xml_out, fuzz_width) ;
	  pop(xml_out);
	  string b_tag("srcLLL_sinkLLL_nucleon") ;
	  ks_compute_baryon(b_tag,quark_propagator_Lsink_Lsrc,
			    quark_propagator_Lsink_Lsrc,
			    quark_propagator_Lsink_Lsrc,
			    xml_out, j_decay, t_length,
			    params.param.binary_baryon_dump,
			    params.param.binary_name) ;

	  pop(xml_out);
	  done_local_baryons = true;
	}

	
      }// do_variational_spectra



      // still have the smeared links if needed for fuzzed loops

      if(( do_fuzzed_disc_loops  )&&(!done_fuzzed_disc_loops)) {
	push(xml_out, "disconnected_loops");

	StopWatch swatch;
	swatch.start();

	ks_fuzz_loops(qprop,q_source, psi ,psi_fuzz , u, u_smr,xml_out, 
		      gauge_shift, sym_shift, loop_checkpoint, t_length, Mass, 
		      Nsamp, RsdCG, CFGNO, volume_source, fuzz_width, 
		      src_seperation, j_decay, 
		      params.param.binary_loop_checkpoint,
		      params.param.binary_name);
	swatch.stop();
	double time_in_sec  = swatch.getTimeInSeconds();
	QDPIO::cout << "PROF9:ks_fuzz_loops " << time_in_sec << " sec" << endl;


/* for testing

	ks_fuzz_loops_stoch_conn(qprop,q_source,psi ,psi_fuzz ,u,u_smr,
				      xml_out, gauge_shift, sym_shift,
				      loop_checkpoint, t_length, Mass, Nsamp,
				      RsdCG, CFGNO, volume_source, fuzz_width, 
				      src_seperation, j_decay);
*/
	pop(xml_out);
	done_local_disc_loops = true;
	done_fuzzed_disc_loops  = true;
      }


    }  // do fuzzing




    //might still need to local pions or baryons

    if((( do_ps4_singlet ) && (!done_ps4_singlet)) ||
       (( do_Baryon_local) &&(!done_local_baryons))  || 
       ( Meson_local  && ! done_meson_corr  )  ||
       (  Meson_charm_local && ! done_meson_charm_corr)
       )
      {

      //still need to compute local corner propagator
      LatticeStaggeredPropagator local_corner_prop;

      ncg_had += 
	MakeCornerProp(psi, gauge_shift, sym_shift, u , qprop, xml_out, RsdCG, 
		       Mass, j_decay, 
                       local_corner_prop, type_of_src, t_source);



      if(( do_ps4_singlet ) && (!done_ps4_singlet)) {
	type_of_src = GAUGE_INVAR_LOCAL_SOURCE ;

	ncg_had += 
	  compute_singlet_ps(psi,local_corner_prop, type_of_src, gauge_shift, 
			     sym_shift, u, qprop, xml_out, RsdCG, Mass, 
			     j_decay, t_source, t_length);

	done_ps4_singlet = true;
      }


      if(( do_Baryon_local) &&(!done_local_baryons)) {

	push(xml_out, "baryon_correlators");

	// describe the source
	string NN ;
	write(xml_out, "source_time", t_source);
	push(xml_out, "smearing_info");
	NN = "L" ;
	write_smearing_info(NN, LOCAL_SRC,xml_out, fuzz_width) ;
	pop(xml_out);
	string b_tag("srcLLL_sinkLLL_nucleon") ;
	ks_compute_baryon(b_tag,local_corner_prop, local_corner_prop,
			  local_corner_prop, xml_out, j_decay, t_length,
			  params.param.binary_baryon_dump,
			  params.param.binary_name) ;

	pop(xml_out);
	done_local_baryons = true;
      }

      // local meson correlators 
      if( !done_meson_corr && Meson_local )
	{
	  push(xml_out, "local_meson_correlators");
	  staggered_local_pion pion(t_length,u) ;
	  pion.compute(local_corner_prop,local_corner_prop,j_decay) ;
	  pion.dump(t_source,xml_out) ;
	  pop(xml_out);
	  done_meson_corr = true ;
	}


      // local meson correlators 
      if( !done_meson_charm_corr && Meson_charm_local )
	{
	  meson_charm(local_corner_prop,
		      xml_out, u,t_source,j_decay,t_length) ;

	  done_meson_charm_corr = true ;
	}


      } // end if-then compute quark propagator


    if( Meson_charm_noisy_local 
	|| Pion_nondegen_noisy_local
	|| Pion_nondegen_noisy_local2
	|| Pion_nondegen_noisy_local3  
	|| Pion_nondegen_noisy_local4  
	)
    {
	int seed = t_source +  t_length * CFGNO  ;
	RNG::setrn(seed);
	QDPIO::cout << "Seeded RNG for noisy time slice source " << seed << endl;
	QDPIO::cout << "Set from t_source " << t_source << " CFG " << CFGNO << "\n" ;
    }


    // compute wall source
    if( Meson_charm_noisy_local  )
      {
	  push(xml_out, "noisy_local_meson_correlators");

      LatticeStaggeredPropagator noisy_corner_prop;

      QDPIO::cout << "Starting inversion for noisy wall source \n" ;
      ncg_had += 
	MakeCornerProp(psi, gauge_shift, sym_shift, u , qprop, xml_out, RsdCG, 
		       Mass, j_decay, 
                       noisy_corner_prop, NOISY_LOCAL_SOURCE , t_source);


      meson_charm(noisy_corner_prop,xml_out, u,t_source,j_decay,t_length) ;

      pop(xml_out); 
      }

    LatticeStaggeredPropagator quark_source_nondegen ; 
    LatticeStaggeredPropagator noisy_corner_prop_strange;


    if( Pion_nondegen_noisy_local  )
      {
	  push(xml_out, "noisy_local_nondegen_meson_correlators");
      QDPIO::cout << "Setting up second inverter for nondegen\n" ;

	// create second inverter for different mass
	XMLReader fermact_reader2 ;
	try{
	  std::istringstream is(params.param.fermact2.xml);
	  fermact_reader2.open(is);
	}
	catch (...)
	  {
	    QDPIO::cerr << "Error reading SECOND action name " << endl;
	    throw;
	  }

	Handle< StaggeredTypeFermAct< T,P,Q> > fermact2(
						       TheStagTypeFermActFactory::Instance().createObject(params.param.fermact.id,
fermact_reader2,params.param.fermact2.path));

	StaggeredTypeFermAct<T,P,Q>& S_f2= *(fermact2);
	Handle< FermState<T,P,Q> > state2(S_f2.createState(u));

	Handle< SystemSolver<LatticeStaggeredFermion> > 
	  qprop2(S_f2.qprop(state, inv_param));

	Real Mass2 = S_f2.getQuarkMass() ;
	Real Mass1 = S_f.getQuarkMass() ;
	

	LatticeStaggeredPropagator noisy_corner_prop_2;
	
	type_of_src = GAUGE_INVAR_LOCAL_SOURCE ;
//	type_of_src = NOISY_LOCAL_SOURCE ; 
    
	QDPIO::cout << "Starting non-degen inversions for local source\n" ;
    ncg_had += 
      MakeCornerProp(psi, gauge_shift, sym_shift, u , qprop, qprop2, 
		       xml_out, RsdCG, 
		       Mass1, Mass2, j_decay, 
                       noisy_corner_prop_strange, noisy_corner_prop_2, 
		       type_of_src , t_source,quark_source_nondegen);

    noisy_pion_nondegen(noisy_corner_prop_strange,Mass1,
			noisy_corner_prop_2,Mass2,
			xml_out, u,t_source,j_decay,t_length) ;

    pop(xml_out); 
      }




    if( Pion_nondegen_noisy_local2  )
      {
	  push(xml_out, "noisy_local_nondegen_meson_correlators_2");
      QDPIO::cout << "Setting up third inverter for nondegen\n" ;

	// create second inverter for different mass
	XMLReader fermact_reader2 ;
	try{
	  std::istringstream is(params.param.fermact3.xml);
	  fermact_reader2.open(is);
	}
	catch (...)
	  {
	    QDPIO::cerr << "Error reading THIRD action name " << endl;
	    throw;
	  }

	Handle< StaggeredTypeFermAct< T,P,Q> > fermact2(
						       TheStagTypeFermActFactory::Instance().createObject(params.param.fermact3.id,
fermact_reader2,params.param.fermact3.path));

	StaggeredTypeFermAct<T,P,Q>& S_f2= *(fermact2);
	Handle< FermState<T,P,Q> > state2(S_f2.createState(u));

	Handle< SystemSolver<LatticeStaggeredFermion> > 
	  qprop2(S_f2.qprop(state, inv_param));

	Real Mass2 = S_f2.getQuarkMass() ;
	Real Mass1 = S_f.getQuarkMass() ;
	
	LatticeStaggeredPropagator noisy_corner_prop_2;
	
	type_of_src = GAUGE_INVAR_LOCAL_SOURCE ;
//	type_of_src = LOAD_IN_SOURCE ;
    
	QDPIO::cout << "Starting 2nd inversions for local source\n" ;
	ncg_had += 
      MakeCornerProp(psi, gauge_shift, sym_shift, u ,qprop2, 
		       xml_out, RsdCG, 
		       Mass2, j_decay, 
                       noisy_corner_prop_2, 
		       type_of_src , t_source,quark_source_nondegen);




    noisy_pion_nondegen(noisy_corner_prop_strange,Mass1,
			noisy_corner_prop_2,Mass2,
			xml_out, u,t_source,j_decay,t_length) ;

    pop(xml_out); 
      }






    if( Pion_nondegen_noisy_local3  )
      {
	  push(xml_out, "noisy_local_nondegen_meson_correlators_3");
      QDPIO::cout << "Setting up fourth inverter for nondegen\n" ;

	// create second inverter for different mass
	XMLReader fermact_reader4 ;
	try{
	  std::istringstream is(params.param.fermact4.xml);
	  fermact_reader4.open(is);
	}
	catch (...)
	  {
	    QDPIO::cerr << "Error reading FOURTH action name " << endl;
	    throw;
	  }

	Handle< StaggeredTypeFermAct< T,P,Q> > fermact2(
						       TheStagTypeFermActFactory::Instance().createObject(params.param.fermact4.id,
fermact_reader4,params.param.fermact4.path));

	StaggeredTypeFermAct<T,P,Q>& S_f2= *(fermact2);
	Handle< FermState<T,P,Q> > state2(S_f2.createState(u));

	Handle< SystemSolver<LatticeStaggeredFermion> > 
	  qprop2(S_f2.qprop(state, inv_param));

	Real Mass2 = S_f2.getQuarkMass() ;
	Real Mass1 = S_f.getQuarkMass() ;
	
	LatticeStaggeredPropagator noisy_corner_prop_2;
	
//	type_of_src = LOAD_IN_SOURCE ;
	type_of_src = GAUGE_INVAR_LOCAL_SOURCE ;
    
	QDPIO::cout << "Starting 3rd inversions for local source\n" ;
	ncg_had += 
      MakeCornerProp(psi, gauge_shift, sym_shift, u ,qprop2, 
		       xml_out, RsdCG, 
		       Mass2, j_decay, 
                       noisy_corner_prop_2, 
		       type_of_src , t_source,quark_source_nondegen);


    noisy_pion_nondegen(noisy_corner_prop_strange,Mass1,
			noisy_corner_prop_2,Mass2,
			xml_out, u,t_source,j_decay,t_length) ;

    pop(xml_out); 
      }





    if( Pion_nondegen_noisy_local4  )
      {
	  push(xml_out, "noisy_local_nondegen_meson_correlators_4");
      QDPIO::cout << "Setting up fifth inverter for nondegen\n" ;

	// create second inverter for different mass
	XMLReader fermact_reader5 ;
	try{
	  std::istringstream is(params.param.fermact5.xml);
	  fermact_reader5.open(is);
	}
	catch (...)
	  {
	    QDPIO::cerr << "Error reading FIFTH action name " << endl;
	    throw;
	  }

	Handle< StaggeredTypeFermAct< T,P,Q> > fermact5(
						       TheStagTypeFermActFactory::Instance().createObject(params.param.fermact5.id,
fermact_reader5,params.param.fermact5.path));

	StaggeredTypeFermAct<T,P,Q>& S_f2= *(fermact5);
	Handle< FermState<T,P,Q> > state2(S_f2.createState(u));

	Handle< SystemSolver<LatticeStaggeredFermion> > 
	  qprop2(S_f2.qprop(state, inv_param));

	Real Mass2 = S_f2.getQuarkMass() ;
	Real Mass1 = S_f.getQuarkMass() ;
	
	LatticeStaggeredPropagator noisy_corner_prop_2;
	
//	type_of_src = LOAD_IN_SOURCE ;
	type_of_src = GAUGE_INVAR_LOCAL_SOURCE ;
	QDPIO::cout << "Starting 4th inversions for local source\n" ;

	ncg_had += 
      MakeCornerProp(psi, gauge_shift, sym_shift, u ,qprop2, 
		       xml_out, RsdCG, 
		       Mass2, j_decay, 
                       noisy_corner_prop_2, 
		       type_of_src , t_source,quark_source_nondegen);

    noisy_pion_nondegen(noisy_corner_prop_strange,Mass1,
			noisy_corner_prop_2,Mass2,
			xml_out, u,t_source,j_decay,t_length) ;

    pop(xml_out); 
      }




  
    pop(xml_out); // spectrum_s

    snoop.stop();
    QDPIO::cout << InlineStaggeredSpectrumEnv::name << ": total time = "
		<< snoop.getTimeInSeconds() 
		<< " secs" << endl;

    QDPIO::cout << InlineStaggeredSpectrumEnv::name << ": ran successfully" << endl;

    END_CODE();
  }  // end of InlineStaggeredSpectrum

}
