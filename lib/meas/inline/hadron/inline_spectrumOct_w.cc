// $Id: inline_spectrumOct_w.cc,v 1.2 2005-05-11 17:27:58 kostas Exp $
/*! \file
 * \brief Inline construction of Octet spectrum
 *
 * Octet spectrum calculations
 */

#include "meas/inline/hadron/inline_spectrumOct_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "meas/smear/ape_smear.h"
#include "meas/smear/sink_smear2.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "io/param_io.h"
#include "io/qprop_io.h"
#include "meas/hadron/mesons_w.h"
#include "meas/hadron/baryon_w.h"
#include "meas/hadron/barhqlq_w.h"
#include "meas/hadron/curcor2_w.h"
#include "meas/hadron/wall_qprop_w.h"
#include "meas/glue/mesfield.h"
#include "util/gauge/taproj.h"
#include "meas/inline/make_xml_file.h"

namespace Chroma 
{ 
  namespace InlineSpectrumOctEnv 
  { 
    AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					    const std::string& path) 
    {
      return new InlineSpectrumOct(InlineSpectrumOctParams(xml_in, path));
    }

    const std::string name = "SPECTRUM_OCT";
    const bool registered = TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
  };



  //! Reader for parameters
  void read(XMLReader& xml, const string& path, InlineSpectrumOctParams::Param_t& param)
  {
    XMLReader paramtop(xml, path);

    int version;
    read(paramtop, "version", version);

    param.HybMesP = false;

    switch (version) 
    {
    case 9:
      param.Wl_snk = false;
      break;

    case 10:
      read(paramtop, "Wl_snk", param.Wl_snk);
      break;

    default:
      QDPIO::cerr << "Input parameter version " << version << " unsupported." << endl;
      QDP_abort(1);
    }

    read(paramtop, "Pt_snk", param.Pt_snk);
    read(paramtop, "Sl_snk", param.Sl_snk);

    read(paramtop, "MesonP", param.MesonP);
    read(paramtop, "CurrentP", param.CurrentP);
    read(paramtop, "BaryonP", param.BaryonP);
    read(paramtop, "time_rev", param.time_rev);

    read(paramtop, "mom2_max", param.mom2_max);
    read(paramtop, "avg_equiv_mom", param.avg_equiv_mom);

    read(paramtop, "wvf_kind", param.wvf_kind);
    read(paramtop, "wvf_param", param.wvf_param);
    read(paramtop, "wvfIntPar", param.wvfIntPar);

    if (param.wvf_param.size() != param.wvfIntPar.size())
    {
      QDPIO::cerr << "wvf_param size inconsistent with wvfintpar size" << endl;
      QDP_abort(1);
    }

    read(paramtop, "nrow", param.nrow);
  }


  //! Writer for parameters
  void write(XMLWriter& xml, const string& path, const InlineSpectrumOctParams::Param_t& param)
  {
    push(xml, path);

    int version = 10;
    write(xml, "version", version);

    write(xml, "Pt_snk", param.Pt_snk);
    write(xml, "Sl_snk", param.Sl_snk);
    write(xml, "Wl_snk", param.Wl_snk);

    write(xml, "MesonP", param.MesonP);
    write(xml, "CurrentP", param.CurrentP);
    write(xml, "BaryonP", param.BaryonP);

    write(xml, "time_rev", param.time_rev);

    write(xml, "mom2_max", param.mom2_max);
    write(xml, "avg_equiv_mom", param.avg_equiv_mom);

    write(xml, "wvf_kind", param.wvf_kind);
    write(xml, "wvf_param", param.wvf_param);
    write(xml, "wvfIntPar", param.wvfIntPar);

    write(xml, "nrow", param.nrow);

    pop(xml);
  }


  //! Propagator input
  void read(XMLReader& xml, const string& path, InlineSpectrumOctParams::Prop_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "prop_files", input.prop_files);
  }

  //! Propagator output
  void write(XMLWriter& xml, const string& path, const InlineSpectrumOctParams::Prop_t& input)
  {
    push(xml, path);

    write(xml, "prop_files", input.prop_files);

    pop(xml);
  }


  // Param stuff
  InlineSpectrumOctParams::InlineSpectrumOctParams() { frequency = 0; }

  InlineSpectrumOctParams::InlineSpectrumOctParams(XMLReader& xml_in, const std::string& path) 
  {
    try 
    {
      XMLReader paramtop(xml_in, path);

      if (paramtop.count("Frequency") == 1)
	read(paramtop, "Frequency", frequency);
      else
	frequency = 1;

      // Parameters for source construction
      read(paramtop, "Param", param);

      // Read in the output propagator/source configuration info
      read(paramtop, "Prop", prop);

      // Possible alternate XML file pattern
      if (paramtop.count("xml_file") != 0) 
      {
	read(paramtop, "xml_file", xml_file);
      }
    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << "Caught Exception reading XML: " << e << endl;
      QDP_abort(1);
    }
  }


  void
  InlineSpectrumOctParams::write(XMLWriter& xml_out, const std::string& path) 
  {
    push(xml_out, path);
    
    Chroma::write(xml_out, "Param", param);
    Chroma::write(xml_out, "Prop", prop);
    QDP::write(xml_out, "xml_file", xml_file);

    pop(xml_out);
  }


  // Function call
  void 
  InlineSpectrumOct::operator()(const multi1d<LatticeColorMatrix>& u,
			     XMLBufferWriter& gauge_xml,
			     unsigned long update_no,
			     XMLWriter& xml_out) 
  {
    // If xml file not empty, then use alternate
    if (params.xml_file != "")
    {
      string xml_file = makeXMLFileName(params.xml_file, update_no);

      push(xml_out, "spectrumOct_w");
      write(xml_out, "update_no", update_no);
      write(xml_out, "xml_file", xml_file);
      pop(xml_out);

      XMLFileWriter xml(xml_file);
      func(u, gauge_xml, update_no, xml);
    }
    else
    {
      func(u, gauge_xml, update_no, xml_out);
    }
  }


  // Real work done here
  void 
  InlineSpectrumOct::func(const multi1d<LatticeColorMatrix>& u,
		       XMLBufferWriter& gauge_xml,
		       unsigned long update_no,
		       XMLWriter& xml_out) 
  {
    push(xml_out, "spectrumOct_w");
    write(xml_out, "update_no", update_no);

    QDPIO::cout << " SPECTRUM_OCT: Spectroscopy for Wilson-like fermions" << endl;


    /*
     * Sanity checks
     */
    if (params.param.wvf_param.size() != params.prop.prop_files.size())
    {
      QDPIO::cerr << "wvf_param size inconsistent with prop_files size" << endl;
      QDP_abort(1);
    }

    QDPIO::cout << endl << "     Gauge group: SU(" << Nc << ")" << endl;

    QDPIO::cout << "     volume: " << params.param.nrow[0];
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
    write(xml_out, "out_version", 11);
    pop(xml_out);


    // First calculate some gauge invariant observables just for info.
    MesPlq(xml_out, "Observables", u);

    // Keep an array of all the xml output buffers
    // XMLArrayWriter xml_array(xml_out,params.prop.prop_files.size());
    // This spectrum code does only one measurement using two propagators
    XMLArrayWriter xml_array(xml_out,1); 
    push(xml_array, "Wilson_hadron_measurements");


    // Fist propagator is the light quark second is the strange quark
    if(params.prop.prop_files.size()!=2){
      QDPIO::cerr << "SPECTRUM_OCT needs two propagators. Found " ;
      QDPIO::cerr <<params.prop.prop_files.size()<< endl;
      QDP_abort(2);
    }

    multi1d<ChromaProp_t> prop_header(params.prop.prop_files.size());
    multi1d<PropSource_t> source_header(params.prop.prop_files.size());
    multi1d<LatticePropagator> qprop(params.prop.prop_files.size());
    multi1d<Real> Mass(params.prop.prop_files.size());
    
    multi2d<int> bc(params.prop.prop_files.size(), 4); 
    
    // Now loop over the various fermion masses
    for (int loop=0; loop < params.prop.prop_files.size(); ++loop)
    {
      // Read the quark propagator and extract headers
      {
	XMLReader prop_file_xml, prop_record_xml;
	readQprop(prop_file_xml, 
		  prop_record_xml, qprop[loop],
		  params.prop.prop_files[loop], QDPIO_SERIAL);

	// Try to invert this record XML into a ChromaProp struct
	// Also pull out the id of this source
	try
	{
	  read(prop_record_xml, "/Propagator/ForwardProp", prop_header[loop]);
	  read(prop_record_xml, "/Propagator/PropSource", source_header[loop]);
	}
	catch (const string& e) 
	{
	  QDPIO::cerr << "Error extracting forward_prop header: " << e << endl;
	  throw;
	}
      }

      // Derived from input prop
       // Hunt around to find the mass
      // NOTE: this may be problematic in the future if actions are used with no
      // clear def. of a Mass
      std::istringstream  xml_s(prop_header[loop].fermact);
      XMLReader  fermacttop(xml_s);
      const string fermact_path = "/FermionAction";
      string fermact;
      
      QDPIO::cout << "Try action and mass" << endl;
      try
      {
	XMLReader top(fermacttop, fermact_path);

	read(top, "FermAct", fermact);

	// Yuk - need to hop some hoops. This should be isolated.
	if (top.count("Mass") != 0) 
	{
	  read(top, "Mass", Mass[loop]);
	}
	else if (top.count("Kappa") != 0)
	{
	  Real Kappa;
	  read(top, "Kappa", Kappa);
	  Mass[loop] = kappaToMass(Kappa);    // Convert Kappa to Mass
	}
	else if (top.count("m_q") != 0) 
	{
	  read(top, "m_q", Mass[loop]);
	}
	else
	{
	  QDPIO::cerr << "Neither Mass nor Kappa found" << endl;
	  throw std::string("Neither Mass nor Kappa found");
	}
      }
      catch (const string& e) 
      {
	QDPIO::cerr << "Error reading fermact or mass: " << e << endl;
	QDP_abort(1);
      }
    
      QDPIO::cout << "FermAct = " << fermact << endl;
      QDPIO::cout << "Mass = " << Mass[loop] << endl;
    }
    // Derived from input prop
    int j_decay = source_header[0].j_decay;
    multi1d<int> t_source = source_header[0].t_source;
    int t0      = t_source[j_decay];
    int bc_spec = bc[0][j_decay] ;
    for (int loop(0); loop < params.prop.prop_files.size(); ++loop){
      if(source_header[loop].j_decay!=j_decay){
	QDPIO::cerr << "Error!! j_decay must be the same for all propagators " << endl;
	throw;
      }
      if(bc[loop][j_decay]!=bc_spec){
	QDPIO::cerr << "Error!! bc must be the same for all propagators " << endl;
	throw;
      }
      for(int d(0);d<Nd;d++)
	if(source_header[loop].t_source[d]!=t_source[d]){
	  QDPIO::cerr << "Error!! t_source must be the same for all propagators " << endl;
	  throw;
	}
    }
  
    // Initialize the slow Fourier transform phases
    SftMom phases(params.param.mom2_max, params.param.avg_equiv_mom, j_decay);

    // Keep a copy of the phases with NO momenta
    SftMom phases_nomom(0, true, j_decay);

    // Next array element - name auto-written
    push(xml_array);
    //write(xml_array, "loop", loop);
    push(xml_out,"Propagator_info") ;
    write(xml_array, "Masses", Mass);
    write(xml_array, "t_source", t_source);
    for (int loop=0; loop < params.prop.prop_files.size(); ++loop){
      push(xml_out, "Propagator");
      // Save prop input
      write(xml_array, "ForwardProp", prop_header[loop]);
      write(xml_array, "PropSource", source_header[loop]);
      multi1d<Double> qp_corr= sumMulti(localNorm2(qprop[loop]),phases.getSet());
      push(xml_out, "Qprop_correlator");
      write(xml_out, "qp_corr", qp_corr);
      pop(xml_out); //Qprop_correlator
      pop(xml_out); //Propagator
    }
    pop(xml_out); //Propagator_info

    // Determine what kind of source to use
    bool Pt_src = false;
    bool Sl_src = false;
    bool Wl_src = false;
    
    //light quark determines the source type....
    //this is a problem that has to be fixed!
    switch (source_header[0].source_type)
      {
      case SRC_TYPE_POINT_SOURCE:
	Pt_src = true;
	break;
	
      case SRC_TYPE_SHELL_SOURCE:
	Sl_src = true;
	break;
	
      case SRC_TYPE_WALL_SOURCE:
	Wl_src = true;
	break;
	
      default:
	QDPIO::cerr << "Unsupported source type" << endl;
	QDP_abort(1);
      }
    
    
    // Do the mesons first
    if (params.param.MesonP) 
      {
	// Construct {Point|Shell}-Point mesons, if desired
	if (params.param.Pt_snk) 
	  {
	    if (Pt_src)
	      mesons(qprop[0], qprop[1], phases, t0,
		     xml_array, "Point_Point_Wilson_Mesons");
	    
	    if (Sl_src)
	      mesons(qprop[0], qprop[1], phases, t0,
		     xml_array, "Shell_Point_Wilson_Mesons");
	    
	    if (Wl_src)
	      mesons(qprop[0], qprop[1],phases_nomom, t0,
		     xml_array, "Wall_Point_Wilson_Mesons");
	  } // end if (Pt_snk)

	// Convolute the quark propagator with the sink smearing function.
	// Make a copy of the quark propagator and then overwrite it with
	// the convolution. 
	if (params.param.Sl_snk) 
	{
	  multi1d<LatticePropagator> quark_prop_smr(params.prop.prop_files.size());
	  for (int loop(0); loop < params.prop.prop_files.size(); ++loop){
	    quark_prop_smr[loop] = qprop[loop];
	    sink_smear2(u, quark_prop_smr[loop], 
			params.param.wvf_kind, 
			params.param.wvf_param[loop],
			params.param.wvfIntPar[loop], 
			j_decay);
	  }

	  if (Pt_src)
	    mesons(quark_prop_smr[0], quark_prop_smr[1], phases, t0,
		   xml_array, "Point_Shell_Wilson_Mesons");

	  if (Sl_src)
	    mesons(quark_prop_smr[0], quark_prop_smr[1], phases, t0,
		   xml_array, "Shell_Shell_Wilson_Mesons");

	  if (Wl_src)
	    mesons(quark_prop_smr[0], quark_prop_smr[1], phases_nomom, t0,
		   xml_array, "Wall_Shell_Wilson_Mesons");
	} // end if (Sl_snk)

	// Wall sink
	if (params.param.Wl_snk) 
	{
	  multi1d<LatticePropagator> wall_quark_prop(params.prop.prop_files.size());
	  for (int loop(0); loop < params.prop.prop_files.size(); ++loop){
	    wall_qprop(wall_quark_prop[loop], qprop[loop], phases_nomom);
	  } 
	  if (Pt_src)
	    mesons(wall_quark_prop[0], wall_quark_prop[1], phases_nomom, t0,
		   xml_array, "Point_Wall_Wilson_Mesons");

	  if (Sl_src)
	    mesons(wall_quark_prop[0], wall_quark_prop[1], phases_nomom, t0,
		   xml_array, "Shell_Wall_Wilson_Mesons");

	  if (Wl_src)
	    mesons(wall_quark_prop[0], wall_quark_prop[1], phases_nomom, t0,
		   xml_array, "Wall_Wall_Wilson_Mesons");
	} // end if (Wl_snk)

      } // end if (MesonP)


      // Do the currents next
      if (params.param.CurrentP) 
      {
	// Construct the rho vector-current and the pion axial current divergence
	if (Pt_src)
	  curcor2(u, qprop[0], qprop[1], phases, 
		  t0, 3,
		  xml_array, "Point_Point_Meson_Currents");
        
	if (Sl_src)
	  curcor2(u, qprop[0], qprop[1], phases, 
		  t0, 3,
		  xml_array, "Shell_Point_Meson_Currents");
        
	if (Wl_src)
	  curcor2(u, qprop[0], qprop[1], phases_nomom, 
		  t0, 3,
		  xml_array, "Wall_Point_Meson_Currents");
      } // end if (CurrentP)


      // Do the baryons
      if (params.param.BaryonP) 
      {
	// Construct {Point|Shell}-Point mesons, if desired
	if (params.param.Pt_snk) 
	{
	  if (Pt_src)
	     barhqlq(qprop[1],qprop[0], phases, 
		   t0, bc_spec, params.param.time_rev, 
		   xml_array, "Point_Point_Wilson_Baryons");
        
	  if (Sl_src)
	     barhqlq(qprop[1],qprop[0], phases, 
		   t0, bc_spec, params.param.time_rev, 
		   xml_array, "Shell_Point_Wilson_Baryons");
        
	  if (Wl_src)
	     barhqlq(qprop[1],qprop[0], phases_nomom, 
		   t0, bc_spec, params.param.time_rev, 
		   xml_array, "Wall_Point_Wilson_Baryons");
	} // end if (Pt_snk)

	// Convolute the quark propagator with the sink smearing function.
	// Make a copy of the quark propagator and then overwrite it with
	// the convolution. 
	if (params.param.Sl_snk) 
	{
	  multi1d<LatticePropagator> quark_prop_smr(params.prop.prop_files.size());
	  for (int loop(0); loop < params.prop.prop_files.size(); ++loop){
	    quark_prop_smr[loop] = qprop[loop];
	    sink_smear2(u, quark_prop_smr[loop], 
			params.param.wvf_kind, 
			params.param.wvf_param[loop],
			params.param.wvfIntPar[loop], 
			j_decay);
	  }

	  if (Pt_src)
	     barhqlq(quark_prop_smr[1],quark_prop_smr[0], phases, 
		   t0, bc_spec, params.param.time_rev, 
		   xml_array, "Point_Shell_Wilson_Baryons");
        
	  if (Sl_src)
	     barhqlq(quark_prop_smr[1],quark_prop_smr[0], phases, 
		   t0, bc_spec, params.param.time_rev, 	
		   xml_array, "Shell_Shell_Wilson_Baryons");
        
	  if (Wl_src)
	     barhqlq(quark_prop_smr[1],quark_prop_smr[0], phases_nomom, 
		   t0, bc_spec, params.param.time_rev, 	
		   xml_array, "Wall_Shell_Wilson_Baryons");
	} // end if (Sl_snk)

	// Wall sink
	if (params.param.Wl_snk) 
	{
	  multi1d<LatticePropagator> wall_quark_prop(params.prop.prop_files.size());
	  for (int loop(0); loop < params.prop.prop_files.size(); ++loop){
	    wall_qprop(wall_quark_prop[loop], qprop[loop], phases_nomom);
	  } 

	  if (Pt_src)
	     barhqlq(wall_quark_prop[1],wall_quark_prop[0] , phases_nomom, 
		   t0, bc_spec, params.param.time_rev, 
		   xml_array, "Point_Wall_Wilson_Baryons");
        
	  if (Sl_src)
	     barhqlq(wall_quark_prop[1],wall_quark_prop[0] , phases_nomom, 
		   t0, bc_spec, params.param.time_rev, 	
		   xml_array, "Shell_Wall_Wilson_Baryons");
        
	  if (Wl_src)
	     barhqlq(wall_quark_prop[1],wall_quark_prop[0] , phases_nomom, 
		   t0, bc_spec, params.param.time_rev, 	
		   xml_array, "Wall_Wall_Wilson_Baryons");
	} // end if (Wl_snk)

      } // end if (BaryonP)

      pop(xml_array);  // array element

      

      pop(xml_array);  // Wilson_spectroscopy
      pop(xml_out);  // spectrumOct_w

      QDPIO::cout << "Wilson Octet spectroscopy ran successfully" << endl;

      END_CODE();
  } 

};
