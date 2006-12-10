// $Id: inline_spectrumQll.cc,v 1.6 2006-12-10 02:02:42 edwards Exp $
/*! \file
 * \brief Inline construction of heavy-light baryon spectrum  
 * (infinitely heavy)
 *
 * Spectrum (Qll) calculations
 */

#error "DEPRECATED - TAKEN OUT OF MAKEFILE"

#include "meas/inline/hadron/inline_spectrumQll.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "meas/smear/ape_smear.h"
#include "meas/smear/sink_smear2.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "io/param_io.h"
#include "io/qprop_io.h"
#include "meas/hadron/barQll_w.h"
#include "meas/hadron/mesQl_w.h"
#include "meas/hadron/wall_qprop_w.h"
#include "meas/glue/mesfield.h"
#include "util/gauge/taproj.h"
#include "meas/inline/make_xml_file.h"

#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{ 
  namespace InlineSpectrumQllEnv 
  { 
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineSpectrumQll(InlineSpectrumQllParams(xml_in, path));
      }
      
      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "SPECTRUM_QLL";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
    }
  }



  //! Reader for parameters
  void read(XMLReader& xml, const string& path, InlineSpectrumQllParams::Param_t& param)
  {
    XMLReader paramtop(xml, path);

    int version;
    read(paramtop, "version", version);

    param.link_smear_fact = 0;
    param.link_smear_num  = 0;

    switch (version) 
    {
    case 1:
      param.Wl_snk = false;
      break;

    case 2:
      read(paramtop, "Wl_snk", param.Wl_snk);
      break;

    case 3:
      read(paramtop, "Wl_snk", param.Wl_snk);
      read(paramtop, "link_smear_fact", param.link_smear_fact);
      read(paramtop, "link_smear_num", param.link_smear_num);
      break;

    default:
      QDPIO::cerr << "Input parameter version " << version << " unsupported." << endl;
      QDP_abort(1);
    }

    read(paramtop, "Pt_snk", param.Pt_snk);
    read(paramtop, "Sl_snk", param.Sl_snk);

    read(paramtop, "wvf_kind", param.wvf_kind);
    read(paramtop, "wvf_param", param.wvf_param);
    read(paramtop, "wvfIntPar", param.wvfIntPar);

    if (param.wvf_param.size() != param.wvfIntPar.size())
    {
      QDPIO::cerr << "wvf_param size inconsistent with wvfintpar size" << endl;
      QDP_abort(1);
    }

    read(paramtop, "Qsrc_coord", param.Qsrc_coord);

  }


  //! Writer for parameters
  void write(XMLWriter& xml, const string& path, const InlineSpectrumQllParams::Param_t& param)
  {
    push(xml, path);

    int version = 3;
    write(xml, "version", version);

    write(xml, "Pt_snk", param.Pt_snk);
    write(xml, "Sl_snk", param.Sl_snk);
    write(xml, "Wl_snk", param.Wl_snk);

    write(xml, "wvf_kind", param.wvf_kind);
    write(xml, "wvf_param", param.wvf_param);
    write(xml, "wvfIntPar", param.wvfIntPar);

    write(xml, "link_smear_fact", param.link_smear_fact);
    write(xml, "link_smear_num", param.link_smear_num);

    write(xml, "Qsrc_coord", param.Qsrc_coord);
    write(xml, "nrow", Layout::lattSize());

    pop(xml);
  }


  //! Propagator input
  void read(XMLReader& xml, const string& path, InlineSpectrumQllParams::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "gauge_id", input.gauge_id);
    read(inputtop, "prop_ids", input.prop_ids);
  }

  //! Propagator output
  void write(XMLWriter& xml, const string& path, const InlineSpectrumQllParams::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "gauge_id", input.gauge_id);
    write(xml, "prop_ids", input.prop_ids);

    pop(xml);
  }


  // Param stuff
  InlineSpectrumQllParams::InlineSpectrumQllParams() { frequency = 0; }

  InlineSpectrumQllParams::InlineSpectrumQllParams(XMLReader& xml_in, const std::string& path) 
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
      read(paramtop, "NamedObject", named_obj);

      // Possible alternate XML file pattern
      if (paramtop.count("xml_file") != 0) 
      {
	read(paramtop, "xml_file", xml_file);
      }
    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << endl;
      QDP_abort(1);
    }
  }


  void
  InlineSpectrumQllParams::write(XMLWriter& xml_out, const std::string& path) 
  {
    push(xml_out, path);
    
    Chroma::write(xml_out, "Param", param);
    Chroma::write(xml_out, "NamedObject", named_obj);
    QDP::write(xml_out, "xml_file", xml_file);

    pop(xml_out);
  }


  // Function call
  void 
  InlineSpectrumQll::operator()(unsigned long update_no,
			     XMLWriter& xml_out) 
  {
    // If xml file not empty, then use alternate
    if (params.xml_file != "")
    {
      string xml_file = makeXMLFileName(params.xml_file, update_no);

      push(xml_out, "spectrumQll_w");
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


  // Real work done here
  void 
  InlineSpectrumQll::func(unsigned long update_no,
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
      QDPIO::cerr << InlineSpectrumQllEnv::name << ": caught dynamic cast error" 
		  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << InlineSpectrumQllEnv::name << ": map call failed: " << e 
		  << endl;
      QDP_abort(1);
    }
    const multi1d<LatticeColorMatrix>& u = 
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

    push(xml_out, "spectrumQll_w");
    write(xml_out, "update_no", update_no);

    QDPIO::cout << InlineSpectrumQllEnv::name << ": Spectroscopy for Wilson-like fermions" << endl;

    /*
     * Sanity checks
     */
    if (params.param.wvf_param.size() != params.named_obj.prop_ids.size())
    {
      QDPIO::cerr << "wvf_param size inconsistent with prop_ids size" << endl;
      QDP_abort(1);
    }

    QDPIO::cout << endl << "     Gauge group: SU(" << Nc << ")" << endl;

    QDPIO::cout << "     volume: " << Layout::lattSize()[0];
    for (int i=1; i<Nd; ++i) {
      QDPIO::cout << " x " << Layout::lattSize()[i];
    }
    QDPIO::cout << endl;

    proginfo(xml_out);    // Print out basic program info

    // Write out the input
    params.write(xml_out, "Input");

    // Write out the config info
    write(xml_out, "Config_info", gauge_xml);

    push(xml_out, "Output_version");
    write(xml_out, "out_version", 2);
    pop(xml_out);


    // First calculate some gauge invariant observables just for info.
    MesPlq(xml_out, "Observables", u);

    // Keep an array of all the xml output buffers
    XMLArrayWriter xml_array(xml_out,params.named_obj.prop_ids.size());
    push(xml_array, "Wilson_hadron_measurements");


    // Now loop over the various fermion masses
    for (int loop=0; loop < params.named_obj.prop_ids.size(); ++loop)
    {
      // Read the quark propagator and extract headers
      ChromaProp_t prop_header;
      PropSourceConst_t source_header;
      QDPIO::cout << "Attempt to read propagator info" << endl;
      try
      {
	// Try the cast to see if this is a valid source
	LatticePropagator& prop_tmp = 
	  TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.prop_ids[loop]);
	
	// Snarf the source info. This is will throw if the source_id is not there
	XMLReader prop_file_xml, prop_record_xml;
	TheNamedObjMap::Instance().get(params.named_obj.prop_ids[loop]).getFileXML(prop_file_xml);
	TheNamedObjMap::Instance().get(params.named_obj.prop_ids[loop]).getRecordXML(prop_record_xml);

	// Try to invert this record XML into a ChromaProp struct
	// Also pull out the id of this source
	{
	  read(prop_record_xml, "/Propagator/ForwardProp", prop_header);
	  read(prop_record_xml, "/Propagator/PropSource", source_header);
	}
      }    
      catch (std::bad_cast)
      {
	QDPIO::cerr << InlineSpectrumQllEnv::name << ": caught dynamic cast error" 
		    << endl;
	QDP_abort(1);
      }
      catch (const string& e) 
      {
	QDPIO::cerr << InlineSpectrumQllEnv::name << ": error extracting prop_header: " << e << endl;
	QDP_abort(1);
      }

      // Should be a valid cast now
      const LatticePropagator& quark_propagator = 
	TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.prop_ids[loop]);
 
      QDPIO::cout << "Propagator successfully read and parsed" << endl;

      // Derived from input prop
      int j_decay           = source_header.j_decay;
      int t0                = source_header.t_source;
      multi1d<int> boundary = getFermActBoundary(prop_header.fermact);
      
      // Hunt around to find the mass
      // NOTE: this may be problematic in the future if actions are used with no
      // clear def. of a Mass
      QDPIO::cout << "Try action and mass" << endl;
      Real Mass = getMass(prop_header.fermact);
    
      QDPIO::cout << "FermAct = " << prop_header.fermact.id << endl;
      QDPIO::cout << "Mass = " << Mass << endl;

      // Flags
      int bc_spec = boundary[j_decay];

      // Initialize the slow Fourier transform phases with NO momenta
      SftMom phases(0, true, j_decay);

      // Next array element - name auto-written
      push(xml_array);
      write(xml_array, "loop", loop);
      write(xml_array, "Mass_mes", Mass);
      write(xml_array, "t0", t0);

      // Save prop input
      write(xml_array, "ForwardProp", prop_header);
      write(xml_array, "PropSource", source_header);

      // Sanity check - write out the norm2 of the forward prop in the j_decay direction
      // Use this for any possible verification
      {
	multi1d<Double> forward_prop_corr = sumMulti(localNorm2(quark_propagator), 
						     phases.getSet());

	push(xml_array, "Forward_prop_correlator");
	write(xml_array, "forward_prop_corr", forward_prop_corr);
	pop(xml_array);
      }

      // Determine what kind of source to use
      bool Pt_src = false;
      bool Sl_src = false;
      bool Wl_src = false;

      if (source_header.source.id == "POINT_SOURCE")
	Pt_src = true;
      else if (source_header.source.id == "SHELL_SOURCE")
	Sl_src = true;
      else if (source_header.source.id == "WALL_SOURCE")
	Wl_src = true;
      else
      {
	QDPIO::cerr << "Unsupported source type" << endl;
	QDP_abort(1);
      }


      /*
       * Smear the gauge field if needed
       */
      multi1d<LatticeColorMatrix> u_link_smr(Nd);
      u_link_smr = u;

      if (params.param.Sl_snk && params.param.link_smear_num > 0)
      {
	int BlkMax = 100;	// Maximum number of blocking/smearing iterations
	Real BlkAccu = 1.0e-5;	// Blocking/smearing accuracy

	for(int i=0; i < params.param.link_smear_num; ++i)
	{
	  multi1d<LatticeColorMatrix> u_tmp(Nd);

	  for(int mu = 0; mu < Nd; ++mu)
	  if ( mu != j_decay )
	    APE_Smear(u_link_smr, u_tmp[mu], mu, 0,
		      params.param.link_smear_fact, BlkAccu, BlkMax, 
		      j_decay);
	  else
	    u_tmp[mu] = u_link_smr[mu];
	  
	  u_link_smr = u_tmp;
	}
	QDPIO::cout << "Gauge field APE-smeared!" << endl;
      }


      // Construct {Point|Shell}-Point, if desired
      if (params.param.Pt_snk) 
	{
	  if (Pt_src){
	    Qll(u_link_smr,quark_propagator,params.param.Qsrc_coord, phases,xml_array, "Point_Point_Wilson_QllBaryons");
	    Qlbar(u_link_smr,quark_propagator,params.param.Qsrc_coord, phases,xml_array, "Point_Point_Wilson_Qlmeson");
	  }
	  if (Sl_src){
	    Qll(u_link_smr,quark_propagator,params.param.Qsrc_coord, phases,xml_array, "Shell_Point_Wilson_QllBaryons");
	    Qlbar(u_link_smr,quark_propagator,params.param.Qsrc_coord, phases,xml_array, "Shell_Point_Wilson_Qlmeson");
	  }
	  if (Wl_src){
	    Qll(u_link_smr,quark_propagator,params.param.Qsrc_coord, phases,xml_array, "Wall_Point_Wilson_QllBaryons");
	    Qlbar(u_link_smr,quark_propagator,params.param.Qsrc_coord, phases,xml_array, "Wall_Point_Wilson_Qlmeson");
	  }
	} // end if (Pt_snk)

      // Convolute the quark propagator with the sink smearing function.
      // Make a copy of the quark propagator and then overwrite it with
      // the convolution. 
      if (params.param.Sl_snk) 
	{
	  LatticePropagator quark_prop_smr;
	  quark_prop_smr = quark_propagator;
	  sink_smear2(u_link_smr, quark_prop_smr, 
		      params.param.wvf_kind, 
		      params.param.wvf_param[loop],
		      params.param.wvfIntPar[loop], 
		      j_decay);
	  {
	    if (Pt_src){
	      Qll(u_link_smr,quark_prop_smr,params.param.Qsrc_coord, phases,xml_array, "Point_Shell_Wilson_QllBaryons");
	      Qlbar(u_link_smr,quark_prop_smr,params.param.Qsrc_coord, phases,xml_array, "Point_Shell_Wilson_Qlmeson");
	    }
	    if (Sl_src){
	      Qll(u_link_smr,quark_prop_smr,params.param.Qsrc_coord, phases,xml_array, "Shell_Shell_Wilson_QllBaryons");
	      Qlbar(u_link_smr,quark_prop_smr,params.param.Qsrc_coord, phases,xml_array, "Shell_Shell_Wilson_Qlmeson");
	    }
	    if (Wl_src){
	      Qll(u_link_smr,quark_prop_smr,params.param.Qsrc_coord, phases,xml_array, "Wall_Shell_Wilson_QllBaryons");
	      Qlbar(u_link_smr,quark_prop_smr,params.param.Qsrc_coord, phases,xml_array, "Wall_Shell_Wilson_Qlmeson");
	    }
	  } 
	} // end if (Sl_snk)

	// Wall sink
      if (params.param.Wl_snk) 
	{
	  LatticePropagator wall_quark_prop;
	  wall_qprop(wall_quark_prop, quark_propagator, phases);
	  if (Pt_src){
	    Qll(u_link_smr,wall_quark_prop,params.param.Qsrc_coord, phases,xml_array, "Point_Wall_Wilson_QllBaryons");
	    Qlbar(u_link_smr,wall_quark_prop,params.param.Qsrc_coord, phases,xml_array, "Point_Wall_Wilson_Qlmeson");
	  }
	  if (Sl_src){
	    Qll(u_link_smr,wall_quark_prop,params.param.Qsrc_coord, phases,xml_array, "Shell_Wall_Wilson_QllBaryons");
	    Qlbar(u_link_smr,wall_quark_prop,params.param.Qsrc_coord, phases,xml_array, "Shell_Wall_Wilson_Qlmeson");
	  }
	  if (Wl_src){
	    Qll(u_link_smr,wall_quark_prop,params.param.Qsrc_coord, phases,xml_array, "Wall_Wall_Wilson_QllBaryons");
	    Qlbar(u_link_smr,wall_quark_prop,params.param.Qsrc_coord, phases,xml_array, "Wall_Wall_Wilson_Qlmeson");
	  }
	} // end if (Wl_snk)


      pop(xml_array);  // array element

    } // end for(loop)

    pop(xml_array);  // Wilson_spectroscopy
    pop(xml_out);  // spectrumQll_w

    snoop.stop();
    QDPIO::cout << InlineSpectrumQllEnv::name << ": total time = "
		<< snoop.getTimeInSeconds() 
		<< " secs" << endl;

    QDPIO::cout << InlineSpectrumQllEnv::name << ": ran successfully" << endl;

    END_CODE();
  }

}; // end namespace Chroma

