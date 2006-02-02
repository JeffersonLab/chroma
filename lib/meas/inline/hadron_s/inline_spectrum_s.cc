// $Id: inline_spectrum_s.cc,v 2.1 2006-02-02 16:22:17 egregory Exp $
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
#include "actions/ferm/fermbcs/fermbcs.h"
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


namespace Chroma 
{ 

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

    read(paramtop, "src_seperation", param.src_seperation);

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
    }


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
    bool do_fuzzing = false ;
    bool do_variational_spectra = false ;
    bool need_basic_8 = false;  // slap this here for now

    //shouldnt be hard-coded
    stag_src_type type_of_src = GAUGE_INVAR_LOCAL_SOURCE ;


    if( params.param.Baryon_vary || params.param.LocalPion_vary){
      do_fuzzing = true ;
    }

    if ( params.param.Baryon_vary || params.param.LocalPion_vary){
      do_variational_spectra = true;
    }

    if( params.param.eight_pions || params.param.eight_scalars || 
	params.param.eight_rhos){
      need_basic_8 = true;
    }



    int ncg_had = 0;    // tally of all the inversion iterations


    bool gauge_shift;
    bool sym_shift;

    gauge_shift = params.param.gauge_invar_oper;
    sym_shift = params.param.sym_shift_oper;


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

    const int j_decay = Nd-1;

    // GAUGE FIX (but what about const)
    // perhaps that should be done outside ?????



    Real RsdCG      = params.prop_param.invParam.RsdCG; // CG residual
    Real Mass       = params.prop_param.Mass;           // fermion mass
    int  fuzz_width = params.param.fuzz_width;          // fuzzing width
    int  t_source   =  params.param.t_srce[Nd-1] ;      // source t coordinate

    // Create a fermion BC
    Handle< FermBC<LatticeStaggeredFermion> >  
      fbc(new SimpleFermBC<LatticeStaggeredFermion>(params.param.boundary));

    //
    // Initialize fermion action
    //
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


    LatticeStaggeredPropagator quark_propagator_Lsink_Lsrc;


    if( params.param.ps4link_singlet_conn ) {
      type_of_src = GAUGE_INVAR_LOCAL_SOURCE ;

      ncg_had += 
	compute_singlet_ps(psi,quark_propagator_Lsink_Lsrc,
			   type_of_src, gauge_shift, sym_shift, 
			   u,qprop,xml_out,
			   params.prop_param.invParam.RsdCG,
			   params.prop_param.Mass,
			   j_decay, t_source, params.param.nrow[3]) ;

    }


    if( params.param.Baryon_local ) {

      push(xml_out, "baryon_correlators");

      // describe the source
      string NN ;
      write(xml_out, "source_time", t_source);
      push(xml_out, "smearing_info");
      NN = "L" ;
      write_smearing_info(NN, LOCAL_SRC,xml_out, params.param.fuzz_width) ;
      pop(xml_out);
      string b_tag("srcLLL_sinkLLL_nucleon") ;
      ks_compute_baryon(b_tag,quark_propagator_Lsink_Lsrc,
			quark_propagator_Lsink_Lsrc,
			quark_propagator_Lsink_Lsrc,
			xml_out, j_decay,
			params.param.nrow[3]) ;

      pop(xml_out);

    }

      LatticeStaggeredFermion psi_fuzz ;


    if( params.param.disconnected_local  ) {
      push(xml_out, "disconnected_loops");
      ks_local_loops(qprop,q_source,psi,u,xml_out,
		     params.param.gauge_invar_oper,
		     params.param.sym_shift_oper,
		     params.param.nrow[3],params.prop_param.Mass,
		     params.param.Nsamp,
		     params.prop_param.invParam.RsdCG,
		     params.param.CFGNO,
		     params.param.volume_source,
		     params.param.src_seperation,
		     j_decay ) ;

      pop(xml_out);
    }


    if (need_basic_8 && !do_fuzzing){
      //Dont need to allocate u_smr here

      multi1d<LatticeStaggeredPropagator> stag_prop(8);

      ncg_had += build_basic_8_props(stag_prop, type_of_src, gauge_shift,
				     sym_shift, u, qprop, xml_out, RsdCG,  
				     Mass, j_decay);

      // put the spectrum calls here 
      if(params.param.eight_pions){
	// do pion stuff
      }
      if(params.param.eight_scalars){
	// do scalar stuff
      }
      if(params.param.eight_rhos){
	// do vector stuff
      }


    }



    // fuzz the gauge configuration. if needed

    if( do_fuzzing  ){

      multi1d<LatticeColorMatrix> u_smr(Nd);

      u_smr = u ;

      DoFuzzing(u, u_smr, j_decay);
    




      // Make basic 8 quark props if needed
      if ( need_basic_8 ){

	// This is inefficient for memory   ???
	multi1d<LatticeStaggeredPropagator> stag_prop(8);

	ncg_had += build_basic_8_props(stag_prop, type_of_src, gauge_shift,
				       sym_shift, fuzz_width, u, u_smr, qprop, 
				       xml_out, RsdCG, Mass, j_decay);

       // put the spectrum calls here 
	if(params.param.eight_pions){
	  // do pion stuff
	}
	if(params.param.eight_scalars){
	  // do scalar stuff
	}
	if(params.param.eight_rhos){
	// do vector stuff
	}
      }




      if( params.param.disconnected_fuzz  ) {
	push(xml_out, "disconnected_loops");
	ks_fuzz_loops(qprop,q_source,psi,psi_fuzz,
		      u,u_smr,xml_out,
		      params.param.gauge_invar_oper,
		      params.param.sym_shift_oper,
		      params.param.nrow[3],params.prop_param.Mass,
		      params.param.Nsamp,
		      params.prop_param.invParam.RsdCG,
		      params.param.CFGNO,
		      params.param.volume_source,
		      params.param.fuzz_width,j_decay) ;


	pop(xml_out);
      }



      if ( do_variational_spectra ){


	// only allocate these if needed

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




	if( params.param.Baryon_vary  ) {
	  compute_vary_baryon_s(xml_out,t_source,
				params.param.fuzz_width,
				j_decay,params.param.nrow[3],
				quark_propagator_Lsink_Lsrc,
				quark_propagator_Fsink_Lsrc,
				quark_propagator_Lsink_Fsrc,
				quark_propagator_Fsink_Fsrc) ;
	}


	if( params.param.LocalPion_vary ) {
	  compute_vary_ps(quark_propagator_Lsink_Lsrc,
			  quark_propagator_Fsink_Lsrc,
			  quark_propagator_Lsink_Fsrc,
			  quark_propagator_Fsink_Fsrc,
			  u,xml_out,j_decay,params.param.nrow[3],t_source) ;
	}


      } // do_variational_spectra

    }  // do fuzzing
    pop(xml_out);  // spectrum_s
    QDPIO::cout << "Staggered spectroscopy ran successfully" << endl;

    END_CODE();
  }

}

