/*! \file
 * \brief Inline construction of staggered spectrum
 *
 * Spectrum calculations
 */

#include "meas/inline/hadron_s/inline_spectrum_s.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "meas/smear/ape_smear.h"
#include "meas/smear/sink_smear2.h"
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



#include "util/ferm/transf.h"
#include "meas/hadron/ks_local_loops.h"
#include "meas/smear/fuzz_smear.h"

#include "meas/inline/make_xml_file.h"

#include "util_compute_quark_prop_s.h"
#include "util_compute_meson_s.h"

namespace Chroma { 
  int build_basic_8_props(multi1d<LatticeStaggeredPropagator> stag_prop,
			   stag_src_type type_of_src,
			   bool gauge_shift, bool sym_shift,
			   int fuzz_width,
			   const multi1d<LatticeColorMatrix> & u,
			   const multi1d<LatticeColorMatrix> & u_smr,
			   Handle<const SystemSolver<LatticeStaggeredFermion> >
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
				 Handle<const SystemSolver<LatticeStaggeredFermion> > & qprop,
				 XMLWriter & xml_out,
				 Real RsdCG, Real Mass, 
				 int j_decay,
				 int src_ind, int color_source);

  int compute_quark_propagator_s(LatticeStaggeredFermion & psi,
				 stag_src_type type_of_src,
				 bool gauge_shift,
				 bool sym_shift,
				 const multi1d<LatticeColorMatrix> & u ,
				 Handle<const SystemSolver<LatticeStaggeredFermion> > & qprop,
				 XMLWriter & xml_out,
				 Real RsdCG, Real Mass, 
				 int j_decay,
				 int src_ind, int color_source);

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


  int compute_singlet_ps(LatticeStaggeredFermion & psi,
			 LatticeStaggeredPropagator quark_propagator,
			 stag_src_type type_of_src,
			 bool gauge_shift,
			 bool sym_shift,
			 const multi1d<LatticeColorMatrix> & u ,
			 Handle<const SystemSolver<LatticeStaggeredFermion> > &
   			     qprop,
			 XMLWriter & xml_out,
			 Real RsdCG, Real Mass, 
			 int j_decay, int t_source, int t_length);

}

// ----------
/***************************************************************************/
/***************************************************************************/

namespace Chroma { 

/***************************************************************************/

  namespace InlineStaggeredSpectrumEnv { 
    AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					    const std::string& path) {
      return new InlineSpectrum_s(InlineSpectrumParams_s(xml_in, path));
    }

    const std::string name = "SPECTRUM_S";
    const bool registered = 
        TheInlineMeasurementFactory::Instance().registerObject(name, 
							       createMeasurement);
  };

/***************************************************************************/


  //! Reader for parameters
  void read(XMLReader& xml, const string& path,
	    InlineSpectrumParams_s::Param_t& param) {
    XMLReader paramtop(xml, path);

    int version;
    read(paramtop, "version", version);

    //    switch (version) 


    read(paramtop, "Meson_local", param.Meson_local);
    read(paramtop, "Baryon_local", param.Baryon_local);
    read(paramtop, "Baryon_vary", param.Baryon_vary);
    read(paramtop, "LocalPion_vary", param.LocalPion_vary);
    read(paramtop, "disconnected_local", param.disconnected_local);
    read(paramtop, "disconnected_fuzz", param.disconnected_fuzz);
    read(paramtop, "singletPs_Conn_local", param.ps4link_singlet_conn);

    read(paramtop, "eight_pions", param.eight_pions);
    read(paramtop, "eight_scalars", param.eight_scalars);
    read(paramtop, "eight_rhos", param.eight_rhos);

    read(paramtop, "boundary", param.boundary);
    read(paramtop, "t_srce", param.t_srce);
    read(paramtop, "nrow", param.nrow);
    read(paramtop, "sym_shift_oper", param.sym_shift_oper);
    read(paramtop, "gauge_invar_oper", param.gauge_invar_oper);
    read(paramtop, "loop_checkpoint", param.loop_checkpoint);

    read(paramtop, "src_seperation", param.src_seperation);

    if( !param.gauge_invar_oper ){
      //read gauge-fixing parameters
      read(paramtop, "GFAccu", param.GFAccu);
      read(paramtop, "OrPara", param.OrPara );
      read(paramtop, "GFMax", param.GFMax );
    }


    if( param.disconnected_local ){
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

    if(param.Baryon_vary){
      read(paramtop, "fuzz_width", param.fuzz_width);
    }else{
      param.fuzz_width = 0 ;
    }
  }

/***************************************************************************/


  //! Writer for parameters
  void write(XMLWriter& xml, const string& path, 
	     const InlineSpectrumParams_s::Param_t& param) {
    push(xml, path);

    int version = 1;
    write(xml, "version", version);

    write(xml, "Meson_local", param.Meson_local);
    write(xml, "Baryon_local", param.Baryon_local);
    write(xml, "Baryon_vary", param.Baryon_vary);
    write(xml, "disconnected_local", param.disconnected_local);
    write(xml, "disconnected_fuzz", param.disconnected_fuzz);
    write(xml, "boundary", param.boundary);
    write(xml, "nrow", param.nrow);
    write(xml, "t_srce", param.t_srce);

    pop(xml);
  }

/***************************************************************************/

  //! Propagator generation params input
  void read(XMLReader& xml, const string& path, 
	    InlineSpectrumParams_s::Quark_Prop_t& input) {
    XMLReader inputtop(xml, path);

    //    read(inputtop, "prop_inversion", input.prop_param);
    read(inputtop, "Mass", input.Mass);
    read(inputtop, "u0", input.u0);
    input.invParam.invType = CG_INVERTER;   
    read(inputtop, "RsdCG", input.invParam.RsdCG);
    read(inputtop, "MaxCG", input.invParam.MaxCG);

  }

/***************************************************************************/


  //! Propagator output
  void write(XMLWriter& xml, const string& path, 
           const InlineSpectrumParams_s::Quark_Prop_t& input) {
    push(xml, path);
    write(xml, "Mass", input.Mass);
    write(xml,"Inverter",input.invParam.invType);
    write(xml,"RsdCG", input.invParam.RsdCG);
    write(xml,"MaxCG", input.invParam.MaxCG);
    pop(xml);
  }

/***************************************************************************/

  // Param stuff
  InlineSpectrumParams_s::InlineSpectrumParams_s() { frequency = 0; }

/***************************************************************************/

  InlineSpectrumParams_s::InlineSpectrumParams_s(XMLReader& xml_in, 
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
  InlineSpectrumParams_s::write(XMLWriter& xml_out, 
				const std::string& path) {
    push(xml_out, path);
    
    Chroma::write(xml_out, "Param", param);
    Chroma::write(xml_out, "Inversion", prop_param);
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
		      Handle<const SystemSolver<LatticeStaggeredFermion> > & 
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
		      Handle<const SystemSolver<LatticeStaggeredFermion> > & 
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
		       Handle<const SystemSolver<LatticeStaggeredFermion> > & 
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
      //      type_of_src = FUZZED_SRC ;

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
		 Handle<const SystemSolver<LatticeStaggeredFermion> > &  qprop,
		 XMLWriter & xml_out,
		 Real RsdCG, Real Mass,
		 int j_decay,
		 LatticeStaggeredPropagator &quark_propagator_Lsink_Lsrc,
		 stag_src_type type_of_src ){

    //    stag_src_type type_of_src = LOCAL_SRC ;
    //    stag_src_type type_of_src = GAUGE_INVAR_LOCAL_SOURCE;

    QDPIO::cout << "LOCAL INVERSIONS"  << endl;
    int ncg_had = 0 ;

    for(int color_source = 0; color_source < Nc; ++color_source){
      psi = zero;    // note this is ``zero'' and not 0

      const int src_ind = 0 ;
      ncg_had += compute_quark_propagator_s(psi,type_of_src, 
					    gauge_shift, sym_shift,
					    u, qprop, xml_out,
					    RsdCG, Mass,
					    j_decay, src_ind, color_source) ;

     /*
       * Move the solution to the appropriate components
       * of quark propagator.
       */
      FermToProp(psi, quark_propagator_Lsink_Lsrc, color_source);

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
  InlineSpectrum_s::operator()(const multi1d<LatticeColorMatrix>& u,
			     XMLBufferWriter& gauge_xml,
			     unsigned long update_no,
			     XMLWriter& xml_out) {
    // If xml file not empty, then use alternate
    if (params.xml_file != ""){
      string xml_file = makeXMLFileName(params.xml_file, update_no);

      push(xml_out, "spectrum_s");
      write(xml_out, "update_no", update_no);
      write(xml_out, "xml_file", xml_file);
      pop(xml_out);

      XMLFileWriter xml(xml_file);
      func(u, gauge_xml, update_no, xml);

    }else{

      func(u, gauge_xml, update_no, xml_out);
    }
  }

/***************************************************************************/

  // Real work done here
  void 
  InlineSpectrum_s::func(const multi1d<LatticeColorMatrix>& u,
		       XMLBufferWriter& gauge_xml,
		       unsigned long update_no,
		       XMLWriter& xml_out) {

    QDPIO::cout << "SPECTRUM_S: Spectroscopy for Staggered-like fermions" 
		<< endl;
    QDPIO::cout << "Gauge group: SU(" << Nc << ")" << endl;

    push(xml_out, "spectrum_s");
    write(xml_out, "update_no", update_no);

    /*  set some flags **/

    bool do_Baryon_local         = params.param.Baryon_local;
    bool do_ps4_singlet          = params.param.ps4link_singlet_conn;

    bool do_Baryon_vary          = params.param.Baryon_vary ;
    bool do_LocalPion_vary       = params.param.LocalPion_vary;
    bool do_8_pions              = params.param.eight_pions;
    bool do_8_scalars            = params.param.eight_scalars;
    bool do_8_rhos               = params.param.eight_rhos;
    bool do_fuzzed_disc_loops    = params.param.disconnected_fuzz  ;
    bool do_local_disc_loops     = params.param.disconnected_local  ;

    bool do_stoch_conn_corr      = true;

    bool do_fuzzing              = false;
    bool do_variational_spectra  = false;

    bool need_basic_8            = false;  
    bool need_local_corner_prop  = false;
    bool need_fuzzed_corner_prop = false;

    bool done_ps4_singlet        = false;
    bool done_local_baryons      = false;
    bool done_local_disc_loops   = false;

    bool have_basic_8            = false;
    bool have_local_corner_prop  = false;
    bool have_fuzzed_corner_prop = false;

    //shouldnt be hard-coded
    stag_src_type type_of_src = GAUGE_INVAR_LOCAL_SOURCE ;


    if( do_Baryon_vary || do_LocalPion_vary){
      do_fuzzing = true ;
    }

    if ( do_Baryon_vary || do_LocalPion_vary){
      do_variational_spectra = true;
    }

    if( do_8_pions || do_8_scalars || do_8_rhos){
      need_basic_8 = true;
    }

    if ((do_Baryon_local) || (do_ps4_singlet)){
      need_local_corner_prop = true;
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
    Real Mass          = params.prop_param.Mass;          // fermion mass
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



    // Create a fermion BC
    Handle< FermBC<LatticeStaggeredFermion> >  
      fbc(new SimpleFermBC<LatticeStaggeredFermion>(params.param.boundary));


    // Initialize fermion action

    AsqtadFermAct S_f(fbc, params.prop_param.Mass,
		      params.prop_param.u0);
    Handle<const ConnectState > state(S_f.createState(u));
    Handle<const SystemSolver<LatticeStaggeredFermion> > 
      qprop(S_f.qprop(state, params.prop_param.invParam));



    //
    //  local inversions
    // 
    LatticeStaggeredFermion psi ;
    LatticeStaggeredFermion q_source ; 


    // Do disconnected loops --- independent of other measurements

    if ((  do_local_disc_loops ) && (do_stoch_conn_corr )){
      push(xml_out, "disconnected_loops");
      ks_local_loops_and_stoch_conn(qprop, q_source, psi , u, xml_out, 
		     gauge_shift, sym_shift, loop_checkpoint,
		     t_length, Mass, Nsamp, RsdCG, CFGNO, volume_source, 
		     src_seperation, j_decay);
      pop(xml_out);

      done_local_disc_loops = true;
    }


    if(( do_local_disc_loops  )  && (!done_local_disc_loops)){
      push(xml_out, "disconnected_loops");

      ks_local_loops(qprop, q_source, psi , u, xml_out, 
		     gauge_shift, sym_shift, loop_checkpoint,
		     t_length, Mass, Nsamp, RsdCG, CFGNO, volume_source, 
		     src_seperation, j_decay);

      done_local_disc_loops = true;

      pop(xml_out);
    }



    if (need_basic_8 ){
      //Dont need to allocate u_smr here

      multi1d<LatticeStaggeredPropagator> stag_prop(8);

      ncg_had += build_basic_8_props(stag_prop, type_of_src, gauge_shift,
				     sym_shift, u, qprop, xml_out, RsdCG,  
				     Mass, j_decay);

      // put the spectrum calls here 
      if(do_8_pions){
	// do pion stuff

	compute_8_pions( stag_prop, u , gauge_shift, sym_shift,
			 xml_out, j_decay, t_length, t_source);
      }
      if(do_8_scalars){
	// do scalar stuff
	compute_8_scalars( stag_prop, u,  gauge_shift, sym_shift,
			   xml_out, j_decay, t_length, t_source);
      }
      if(do_8_rhos){
	// do vector stuff
	compute_8_vectors( stag_prop, u,  gauge_shift, sym_shift,
			   xml_out, j_decay, t_length, t_source);
      }


      // if we need to do 4-link singlets and local baryons, do them here
      // so we can re-use the stag_pro[0] as the local corner prop

      if(( do_ps4_singlet ) && (!done_ps4_singlet)) {
	type_of_src = GAUGE_INVAR_LOCAL_SOURCE ;

	ncg_had += 
	  compute_singlet_ps(psi,stag_prop[0], type_of_src, gauge_shift, 
			     sym_shift, u, qprop, xml_out, RsdCG,Mass,
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
	ks_compute_baryon(b_tag,stag_prop[0],stag_prop[0],stag_prop[0],
			  xml_out, j_decay, t_length) ;

	pop(xml_out);
	done_local_baryons = true;
      }

 
    } // basic 8 and no fuzzing



    // fuzz the gauge configuration. if needed

    if( do_fuzzing  ){

      LatticeStaggeredFermion psi_fuzz ;
      multi1d<LatticeColorMatrix> u_smr(Nd);

      u_smr = u ;
      DoFuzzing(u, u_smr, j_decay);
    

      if ( do_variational_spectra ){

	// only allocate these if needed
	LatticeStaggeredPropagator quark_propagator_Lsink_Lsrc;
	LatticeStaggeredPropagator quark_propagator_Fsink_Lsrc ;
	LatticeStaggeredPropagator quark_propagator_Lsink_Fsrc ;
	LatticeStaggeredPropagator quark_propagator_Fsink_Fsrc ;


	ncg_had+=  MakeFuzzedCornerProp(psi, fuzz_width, gauge_shift, 
					sym_shift, u, u_smr, qprop, 
					xml_out, RsdCG, Mass, j_decay,
					do_fuzzing, psi_fuzz,
					quark_propagator_Lsink_Lsrc,
					quark_propagator_Fsink_Lsrc,
					quark_propagator_Lsink_Fsrc,
					quark_propagator_Fsink_Fsrc );


	//re-use quark_propagator_Lsink_Lsrc here for 
	// 4-link singlets and local baryons.

	if(( do_ps4_singlet ) && (!done_ps4_singlet)) {

	  type_of_src = GAUGE_INVAR_LOCAL_SOURCE ;

	  ncg_had += 
	    compute_singlet_ps(psi,quark_propagator_Lsink_Lsrc, type_of_src, 
			       gauge_shift, sym_shift, u, qprop, xml_out, 
			       RsdCG, Mass, j_decay, t_source, t_length);

	  done_ps4_singlet = true;
	}


	if(( do_Baryon_local ) && (!done_local_baryons)) {

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
			    xml_out, j_decay,
			    t_length) ;

	  pop(xml_out);
	  done_local_baryons = true;
	}



	if( do_Baryon_vary ) {
	  compute_vary_baryon_s(xml_out,t_source,
				fuzz_width,
				j_decay,t_length,
				quark_propagator_Lsink_Lsrc,
				quark_propagator_Fsink_Lsrc,
				quark_propagator_Lsink_Fsrc,
				quark_propagator_Fsink_Fsrc) ;
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



	// if we need to do 4-link singlets and local baryons, do them here so
	// we can re-use quark_propagator_Lsink_Lsrc as the local corner prop

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
			    xml_out, j_decay, t_length) ;

	  pop(xml_out);
	  done_local_baryons = true;
	}

	
      }// do_variational_spectra



      if(( do_fuzzed_disc_loops  )) {
	push(xml_out, "disconnected_loops");
	ks_fuzz_loops(qprop,q_source, psi ,psi_fuzz , u, u_smr,xml_out, 
		      gauge_shift, sym_shift, loop_checkpoint, t_length, Mass, 
		      Nsamp, RsdCG, CFGNO, volume_source, fuzz_width, 
		      src_seperation, j_decay);

	pop(xml_out);
      }


    }  // do fuzzing




    //might still need to local pions or baryons

    if((( do_ps4_singlet ) && (!done_ps4_singlet)) ||
       (( do_Baryon_local) &&(!done_local_baryons))){

      //still need to compute local corner propagator
      LatticeStaggeredPropagator local_corner_prop;

      ncg_had += 
	MakeCornerProp(psi, gauge_shift, sym_shift, u , qprop, xml_out, RsdCG, 
		       Mass, j_decay, local_corner_prop, type_of_src);



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
			  local_corner_prop, xml_out, j_decay, t_length) ;

	pop(xml_out);
	done_local_baryons = true;
      }
    }

  
  
    pop(xml_out); // spectrum_s
    QDPIO::cout << "Staggered spectroscopy ran successfully" << endl;

    END_CODE();
  }  // end of InlineSpectrum_s

}
