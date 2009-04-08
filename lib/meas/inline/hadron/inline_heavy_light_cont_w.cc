// $Id: inline_heavy_light_cont_w.cc,v 3.7 2009-04-08 18:34:11 caubin Exp $
/*! \file
 * \brief Inline construction of hadron contractions
 *  Still just does static!!
 * Contraction calculations
 */

#include "meas/inline/hadron/inline_heavy_light_cont_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "io/param_io.h"
#include "io/qprop_io.h"
#include "meas/hadron/mesQl_w.h"
#include "meas/hadron/Ql_3pt_w.h"
#include "meas/hadron/curcor2_w.h"
#include "meas/inline/make_xml_file.h"
#include "meas/inline/io/named_objmap.h"
#include "meas/smear/no_quark_displacement.h"

namespace Chroma 
{ 
  namespace InlineHeavyLightContEnv 
  { 
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineHeavyLightCont(InlineHeavyLightContParams(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "HEAVY_LIGHT_CONTRACTION";

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
  void read(XMLReader& xml, const string& path, InlineHeavyLightContParams::Param_t& param)
  {
    XMLReader paramtop(xml, path);

    int version;
    read(paramtop, "version", version);

    switch (version) 
    {
    case 1:
      break;

    default:
      QDPIO::cerr << "Input parameter version " << version << " unsupported." << endl;
      QDP_abort(1);
    }

    read(paramtop, "MesonP", param.MesonP);
    read(paramtop, "FourPt", param.FourPt);

  }


  //! Writer for parameters
  void write(XMLWriter& xml, const string& path, const InlineHeavyLightContParams::Param_t& param)
  {
    push(xml, path);

    int version = 1;
    write(xml, "version", version);

    write(xml, "MesonP", param.MesonP);
    write(xml, "FourPt", param.FourPt);
    write(xml, "nrow", Layout::lattSize());

    pop(xml);
  }


  //! Propagator input
  void read(XMLReader& xml, const string& path, InlineHeavyLightContParams::NamedObject_t::Props_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "first_id", input.first_id);
    read(inputtop, "second_id", input.second_id);
    if (inputtop.count("third_id") == 1)
      read(inputtop, "third_id", input.third_id);
    else
      input.third_id = "None";

    if (inputtop.count("heavy_id1") == 1)
      read(inputtop, "heavy_id1", input.heavy_id1);
    else
      input.heavy_id1 = "Static";
    if (inputtop.count("heavy_id2") == 1)
      read(inputtop, "heavy_id2", input.heavy_id2);
    else
      input.heavy_id2 = "Static";
  }

  //! Propagator output
  void write(XMLWriter& xml, const string& path, const InlineHeavyLightContParams::NamedObject_t::Props_t& input)
  {
    push(xml, path);

    write(xml, "first_id", input.first_id);
    write(xml, "second_id", input.second_id);
    write(xml, "third_id", input.third_id);

    write(xml, "heavy_id1", input.heavy_id1);
    write(xml, "heavy_id2", input.heavy_id2);
    pop(xml);
  }


  //! Propagator input
  void read(XMLReader& xml, const string& path, InlineHeavyLightContParams::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "gauge_id", input.gauge_id);
    read(inputtop, "sink_pairs", input.sink_pairs);
  }

  //! Propagator output
  void write(XMLWriter& xml, const string& path, const InlineHeavyLightContParams::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "gauge_id", input.gauge_id);
    write(xml, "sink_pairs", input.sink_pairs);

    pop(xml);
  }


  // Param stuff
  InlineHeavyLightContParams::InlineHeavyLightContParams()
  { 
    frequency = 0; 
  }

  InlineHeavyLightContParams::InlineHeavyLightContParams(XMLReader& xml_in, const std::string& path) 
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
  InlineHeavyLightContParams::write(XMLWriter& xml_out, const std::string& path) 
  {
    push(xml_out, path);
    
    Chroma::write(xml_out, "Param", param);
    Chroma::write(xml_out, "NamedObject", named_obj);
    QDP::write(xml_out, "xml_file", xml_file);

    pop(xml_out);
  }

  // Anonymous namespace
  namespace 
  {
    //! Useful structure holding sink props
    struct SinkPropContainer_t
    {
      ForwardProp_t prop_header;
      string quark_propagator_id;
      Real Mass;
    
      multi1d<int> bc; 
    
      string source_type;
      string source_disp_type;
      string sink_type;
      string sink_disp_type;
    };


    //! Useful structure holding sink props
    struct AllSinkProps_t
    {
      SinkPropContainer_t  sink_prop_1;
      SinkPropContainer_t  sink_prop_2;
      SinkPropContainer_t  sink_prop_3;

      SinkPropContainer_t  heavy_prop_1;
      SinkPropContainer_t  heavy_prop_2;
    };


    //! Read a sink prop
    void readSinkProp(SinkPropContainer_t& s, const std::string& id)
    {
      try
      {
	// Try a cast to see if it succeeds
	const LatticePropagator& foo = 
	  TheNamedObjMap::Instance().getData<LatticePropagator>(id);

	// Snarf the data into a copy
	s.quark_propagator_id = id;
	
	// Snarf the prop info. This is will throw if the prop_id is not there
	XMLReader prop_file_xml, prop_record_xml;
	TheNamedObjMap::Instance().get(id).getFileXML(prop_file_xml);
	TheNamedObjMap::Instance().get(id).getRecordXML(prop_record_xml);
   
	// Try to invert this record XML into a ChromaProp struct
	// Also pull out the id of this source
	{
	  string xpath;
	  read(prop_record_xml, "/SinkSmear", s.prop_header);
	  
	  read(prop_record_xml, "/SinkSmear/PropSource/Source/SourceType", s.source_type);
	  xpath = "/SinkSmear/PropSource/Source/Displacement/DisplacementType";
	  if (prop_record_xml.count(xpath) != 0)
	    read(prop_record_xml, xpath, s.source_disp_type);
	  else
	    s.source_disp_type = NoQuarkDisplacementEnv::getName();

	  read(prop_record_xml, "/SinkSmear/PropSink/Sink/SinkType", s.sink_type);
	  xpath = "/SinkSmear/PropSink/Sink/Displacement/DisplacementType";
	  if (prop_record_xml.count(xpath) != 0)
	    read(prop_record_xml, xpath, s.sink_disp_type);
	  else
	    s.sink_disp_type = NoQuarkDisplacementEnv::getName();
	}
      }
      catch( std::bad_cast ) 
      {
	QDPIO::cerr << InlineHeavyLightContEnv::name << ": caught dynamic cast error" 
		    << endl;
	QDP_abort(1);
      }
      catch (const string& e) 
      {
	QDPIO::cerr << InlineHeavyLightContEnv::name << ": error message: " << e 
		    << endl;
	QDP_abort(1);
      }

      // Derived from input prop
      // Hunt around to find the mass
      // NOTE: this may be problematic in the future if actions are used with no
      // clear def. of a Mass
      QDPIO::cout << "Try action and mass" << endl;
      s.Mass = getMass(s.prop_header.prop_header.fermact);
      s.bc = getFermActBoundary(s.prop_header.prop_header.fermact);

      QDPIO::cout << "FermAct = " << s.prop_header.prop_header.fermact.id << endl;
      QDPIO::cout << "Mass = " << s.Mass << endl;
    }


    //! Read all sinks
    void readAllSinks(AllSinkProps_t& s, 
		      InlineHeavyLightContParams::NamedObject_t::Props_t sink_pair, 
		      InlineHeavyLightContParams::Param_t params)
    {
      QDPIO::cout << "Attempt to parse forward propagator = " << sink_pair.first_id << endl;
      readSinkProp(s.sink_prop_1, sink_pair.first_id);
      QDPIO::cout << "Forward propagator successfully parsed" << endl;

      QDPIO::cout << "Attempt to parse forward propagator = " << sink_pair.second_id << endl;
      readSinkProp(s.sink_prop_2, sink_pair.second_id);
      QDPIO::cout << "Forward propagator successfully parsed" << endl;

      if (sink_pair.third_id != "None"){
	QDPIO::cout << "Attempt to parse forward propagator = " << sink_pair.third_id << endl;
	readSinkProp(s.sink_prop_3, sink_pair.third_id);
	QDPIO::cout << "Forward propagator successfully parsed" << endl;
      }
      else{
	QDPIO::cout << "Using "<<sink_pair.second_id <<" as dummy spectator." << endl;
        readSinkProp(s.sink_prop_3, sink_pair.second_id);
	QDPIO::cout << "Forward propagator successfully parsed" << endl;
      }
      if (sink_pair.heavy_id1 != "Static"){
	QDPIO::cout << "Attempt to parse forward propagator = " << sink_pair.heavy_id1 << endl;
        readSinkProp(s.heavy_prop_1, sink_pair.heavy_id1);
	QDPIO::cout << "Forward propagator successfully parsed" << endl;
      }
      if (sink_pair.heavy_id2 != "Static"){
	QDPIO::cout << "Attempt to parse forward propagator = " << sink_pair.heavy_id2 << endl;
        readSinkProp(s.heavy_prop_2, sink_pair.heavy_id2);
	QDPIO::cout << "Forward propagator successfully parsed" << endl;
      }
    }
  } // namespace anonymous



  // Function call
  void 
  InlineHeavyLightCont::operator()(unsigned long update_no,
			    XMLWriter& xml_out) 
  {
    // If xml file not empty, then use alternate
    if (params.xml_file != "")
    {
      string xml_file = makeXMLFileName(params.xml_file, update_no);

      push(xml_out, "HeavyLightCont");
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
  InlineHeavyLightCont::func(unsigned long update_no,
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
      QDPIO::cerr << InlineHeavyLightContEnv::name << ": caught dynamic cast error" 
		  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << InlineHeavyLightContEnv::name << ": map call failed: " << e 
		  << endl;
      QDP_abort(1);
    }
    const multi1d<LatticeColorMatrix>& u = 
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

    push(xml_out, "HeavyLightCont");
    write(xml_out, "update_no", update_no);

    QDPIO::cout << " HeavyLightCont: Contractions for Wilson-like fermions" ;
    QDPIO::cout << endl;
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
    write(xml_out, "out_version", 14);
    pop(xml_out);

    // First calculate some gauge invariant observables just for info.
    MesPlq(xml_out, "Observables", u);

    // Keep an array of all the xml output buffers
    push(xml_out, "Wilson_hadron_measurements");

    // Now loop over the various fermion pairs
	// Here we want to make the separate versions of the code.
	// There are three: With two static, with one static, or with
	// zero static heavy quarks. 
    for(int lpair=0; lpair < params.named_obj.sink_pairs.size(); ++lpair)
    {
      // Note that this named_obj should have the defs of the props...
      const InlineHeavyLightContParams::NamedObject_t::Props_t named_obj = params.named_obj.sink_pairs[lpair];
      const InlineHeavyLightContParams::Param_t param = params.param;
	  
      push(xml_out, "elem");
	  
      AllSinkProps_t all_sinks;
      readAllSinks(all_sinks, named_obj, param);
	  
      // Derived from input prop
      multi1d<int> t_src1
		= all_sinks.sink_prop_1.prop_header.source_header.getTSrce();
      multi1d<int> t_src2
		= all_sinks.sink_prop_2.prop_header.source_header.getTSrce();
      // Note the third prop must be here, but is ignored later if not doing four-point...
      multi1d<int> t_src3
        = all_sinks.sink_prop_3.prop_header.source_header.getTSrce();
      
      int j_decay = all_sinks.sink_prop_1.prop_header.source_header.j_decay;
      int t01      = all_sinks.sink_prop_1.prop_header.source_header.t_source;
      int t02      = all_sinks.sink_prop_2.prop_header.source_header.t_source;
      int bc_spec = all_sinks.sink_prop_1.bc[j_decay] ;
	  
      // Sanity checks
      {
	if (all_sinks.sink_prop_2.prop_header.source_header.j_decay != j_decay)
	  {
	    QDPIO::cerr << "Error!! j_decay must be the same for all propagators " << endl;
	    QDP_abort(1);
	  }
	/**
	if (named_obj.heavy_id1 != "Static"){
	if ((all_sinks.heavy_prop_1.prop_header.source_header.t_source != t01) 
	|| (all_sinks.heavy_prop_1.prop_header.source_header.t_source != t02))
	{
	QDPIO::cerr << "Error!! Heavy quark 1 must have same t_source as one of the light propagators " << endl;
	QDP_abort(1);
	}
	if (all_sinks.heavy_prop_1.prop_header.source_header.j_decay != j_decay)
	{
	QDPIO::cerr << "Error!! j_decay must be the same for all propagators " << endl;
	QDP_abort(1);
	}
	}
	if (named_obj.heavy_id2 != "Static"){
	if ((all_sinks.heavy_prop_2.prop_header.source_header.t_source != t01) 
	|| (all_sinks.heavy_prop_2.prop_header.source_header.t_source != t02))
	{
	QDPIO::cerr << "Error!! Heavy quark 2 must have same t_source as one of the light propagators " << endl;
	QDP_abort(1);
	}
	  if (all_sinks.heavy_prop_2.prop_header.source_header.j_decay != j_decay)
	  {
	  QDPIO::cerr << "Error!! j_decay must be the same for all propagators " << endl;
	  QDP_abort(1);
	  }
	  }
	
	if (all_sinks.sink_prop_2.bc[j_decay] != bc_spec)
	  {
	    QDPIO::cerr << "Error!! bc must be the same for all propagators " << endl;
	    QDP_abort(1);
	  }
	if (params.param.FourPt){
	  if (all_sinks.sink_prop_3.prop_header.source_header.j_decay != j_decay)
	    {
	      QDPIO::cerr << "Error!! j_decay must be the same for all propagators " << endl;
	      QDP_abort(1);
	    }
	  if (all_sinks.sink_prop_3.bc[j_decay] != bc_spec)
	    {
	      QDPIO::cerr << "Error!! bc must be the same for all propagators " << endl;
	      QDP_abort(1);
	    }
	}
	**/
      }
      
      /** 
	  Here we want to define the bc's for everything...in principle this 
	  should be read in from everything, but whatever.
       **/
      
      int HQBC = -1;
      
      if (named_obj.heavy_id1=="Static" && named_obj.heavy_id2=="Static"){
	
	// Masses
	write(xml_out, "Mass_1", all_sinks.sink_prop_1.Mass);
	write(xml_out, "Mass_2", all_sinks.sink_prop_2.Mass);
	if (params.param.FourPt)
	  write(xml_out, "Mass_3", all_sinks.sink_prop_3.Mass);
	write(xml_out, "t01", t01);
	write(xml_out, "t02", t02);
	
	// Initialize the slow Fourier transform phases with NO momenta
	SftMom phases(0, true, j_decay);
	
	// Save prop input
	push(xml_out, "Forward_prop_headers");
	write(xml_out, "First_forward_prop", all_sinks.sink_prop_1.prop_header);
	write(xml_out, "Second_forward_prop", all_sinks.sink_prop_2.prop_header);
	if (params.param.FourPt)
	  write(xml_out, "Third_forward_prop", all_sinks.sink_prop_3.prop_header);
	pop(xml_out);
	
	// Sanity check - write out the norm2 of the forward prop in the j_decay direction
	// Use this for any possible verification
	push(xml_out, "Forward_prop_correlator");
	{
	  const LatticePropagator& sink_prop_1 = 
	    TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks.sink_prop_1.quark_propagator_id);
	  const LatticePropagator& sink_prop_2 = 
	    TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks.sink_prop_2.quark_propagator_id);
	  const LatticePropagator& sink_prop_3 = 
	    TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks.sink_prop_3.quark_propagator_id);
	  
	  write(xml_out, "forward_prop_corr_1", sumMulti(localNorm2(sink_prop_1), phases.getSet()));
	  write(xml_out, "forward_prop_corr_2", sumMulti(localNorm2(sink_prop_2), phases.getSet()));
	  
	  write(xml_out, "forward_prop_corr_3", sumMulti(localNorm2(sink_prop_3), phases.getSet()));
	}
	pop(xml_out);
		
	
	push(xml_out, "SourceSinkType");
	{
	  QDPIO::cout << "Source_type_1 = " << all_sinks.sink_prop_1.source_type << endl;
	  QDPIO::cout << "Sink_type_1 = " << all_sinks.sink_prop_1.sink_type << endl;
	  QDPIO::cout << "Source_type_2 = " << all_sinks.sink_prop_2.source_type << endl;
	  QDPIO::cout << "Sink_type_2 = " << all_sinks.sink_prop_2.sink_type << endl;
	  if (params.param.FourPt){
	    QDPIO::cout << "Source_type_3 = " << all_sinks.sink_prop_3.source_type << endl;
	    QDPIO::cout << "Sink_type_3 = " << all_sinks.sink_prop_3.sink_type << endl;
	  }
	  write(xml_out, "source_type_1", all_sinks.sink_prop_1.source_type);
	  write(xml_out, "source_disp_type_1", all_sinks.sink_prop_1.source_disp_type);
	  write(xml_out, "sink_type_1", all_sinks.sink_prop_1.sink_type);
	  write(xml_out, "sink_disp_type_1", all_sinks.sink_prop_1.sink_disp_type);
	  
	  write(xml_out, "source_type_2", all_sinks.sink_prop_2.source_type);
	  write(xml_out, "source_disp_type_2", all_sinks.sink_prop_2.source_disp_type);
	  write(xml_out, "sink_type_2", all_sinks.sink_prop_2.sink_type);
	  write(xml_out, "sink_disp_type_2", all_sinks.sink_prop_2.sink_disp_type);
	  if (params.param.FourPt){
	    write(xml_out, "source_type_3", all_sinks.sink_prop_3.source_type);
	    write(xml_out, "source_disp_type_3", all_sinks.sink_prop_3.source_disp_type);
	    write(xml_out, "sink_type_3", all_sinks.sink_prop_3.sink_type);
	    write(xml_out, "sink_disp_type_3", all_sinks.sink_prop_3.sink_disp_type);
	  }
	}
	pop(xml_out);
	
	
	// References for use later
	const LatticePropagator& sink_prop_1 = 
	  TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks.sink_prop_1.quark_propagator_id);
	const LatticePropagator& sink_prop_2 = 
	  TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks.sink_prop_2.quark_propagator_id);
	const LatticePropagator& sink_prop_3 = 
	  TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks.sink_prop_3.quark_propagator_id);
	
	// Construct group name for output
	string src_type;
	if (all_sinks.sink_prop_1.source_type == "POINT_SOURCE")
	  src_type = "Point";
	else if (all_sinks.sink_prop_1.source_type == "SHELL_SOURCE")
	  src_type = "Shell";
	else if (all_sinks.sink_prop_1.source_type == "WALL_SOURCE")
	  src_type = "Wall";
	else
	  {
	    QDPIO::cerr << "Unsupported source type = " << all_sinks.sink_prop_1.source_type << endl;
	    QDP_abort(1);
	  }
	
	string snk_type;
	if (all_sinks.sink_prop_1.sink_type == "POINT_SINK")
	  snk_type = "Point";
	else if (all_sinks.sink_prop_1.sink_type == "SHELL_SINK")
	  snk_type = "Shell";
	else if (all_sinks.sink_prop_1.sink_type == "WALL_SINK")
	  snk_type = "Wall";
	else
	  {
	    QDPIO::cerr << "Unsupported sink type = " << all_sinks.sink_prop_1.sink_type << endl;
	    QDP_abort(1);
	  }
	
	string source_sink_type = src_type + "_" + snk_type;
	QDPIO::cout << "Source type = " << src_type << endl;
	QDPIO::cout << "Sink type = "   << snk_type << endl;
	
	if (params.param.MesonP) 
	  {
	    Qlbar(u, sink_prop_1, t_src1, phases, xml_out, source_sink_type+"_Wilson_Qlmeson1",HQBC);
	    if(all_sinks.sink_prop_1.quark_propagator_id!=all_sinks.sink_prop_2.quark_propagator_id)
	      Qlbar(u, sink_prop_2, t_src2, phases, xml_out, source_sink_type+"_Wilson_Qlmeson2",HQBC);
	    if (params.param.FourPt && all_sinks.sink_prop_3.quark_propagator_id!=all_sinks.sink_prop_2.quark_propagator_id && all_sinks.sink_prop_1.quark_propagator_id!=all_sinks.sink_prop_3.quark_propagator_id)
	      Qlbar(u, sink_prop_3, t_src3, phases, xml_out, source_sink_type + "_Wilson_Qlmeson3",HQBC);
	    // For testing backwards prop, we will calculate the Qlbar for a backwards
	    // static quark
	    QlbarBACK(u, sink_prop_1, t_src1, phases, xml_out, source_sink_type+"_Wilson_Ql_back_meson1",HQBC);
            if(all_sinks.sink_prop_1.quark_propagator_id!=all_sinks.sink_prop_2.quark_propagator_id)
              QlbarBACK(u, sink_prop_2, t_src2, phases, xml_out, source_sink_type+"_Wilson_Ql_back_meson2",HQBC);
            if (params.param.FourPt && all_sinks.sink_prop_3.quark_propagator_id!=all_sinks.sink_prop_2.quark_propagator_id && all_sinks.sink_prop_1.quark_propagator_id!=all_sinks.sink_prop_3.quark_propagator_id)
              QlbarBACK(u, sink_prop_3, t_src3, phases, xml_out, source_sink_type + "_Wilson_Ql_back_meson3",HQBC);

	    
	  } // end if (MesonP)
	
	//Now we actually will do the contractions...
		
	QlQl(u, sink_prop_1, sink_prop_2,  t_src1, t_src2, HQBC, phases, xml_out, source_sink_type +"Wilson_Ql_Ql_3pt");
	
	pop(xml_out);  // array element
      }// End of if statement for both heavies being static.
      
      else if (named_obj.heavy_id1=="Static" && named_obj.heavy_id2!="Static"){
	// Masses
	write(xml_out, "Mass_1", all_sinks.sink_prop_1.Mass);
	write(xml_out, "Mass_2", all_sinks.sink_prop_2.Mass);
	if (params.param.FourPt)
	  write(xml_out, "Mass_3", all_sinks.sink_prop_3.Mass);
	write(xml_out, "t01", t01);
	write(xml_out, "t02", t02);
	
	// Initialize the slow Fourier transform phases with NO momenta
	SftMom phases(0, true, j_decay);
	
	// Save prop input
	push(xml_out, "Forward_prop_headers");
	write(xml_out, "First_forward_prop", all_sinks.sink_prop_1.prop_header);
	write(xml_out, "Second_forward_prop", all_sinks.sink_prop_2.prop_header);
	if (params.param.FourPt)
	  write(xml_out, "Third_forward_prop", all_sinks.sink_prop_3.prop_header);
	pop(xml_out);
	
	// Sanity check - write out the norm2 of the forward prop in the j_decay direction
	// Use this for any possible verification
	push(xml_out, "Forward_prop_correlator");
	{
	  const LatticePropagator& sink_prop_1 = 
	    TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks.sink_prop_1.quark_propagator_id);
	  const LatticePropagator& sink_prop_2 = 
	    TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks.sink_prop_2.quark_propagator_id);
	  const LatticePropagator& sink_prop_3 = 
	    TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks.sink_prop_3.quark_propagator_id);
	  
	  const LatticePropagator& heavy_prop_2 = 
	    TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks.heavy_prop_2.quark_propagator_id);
	  
	  write(xml_out, "forward_prop_corr_1", sumMulti(localNorm2(sink_prop_1), phases.getSet()));
	  write(xml_out, "forward_prop_corr_2", sumMulti(localNorm2(sink_prop_2), phases.getSet()));
	  
	  write(xml_out, "forward_prop_corr_3", sumMulti(localNorm2(sink_prop_3), phases.getSet()));
	  write(xml_out, "forward_prop_corr_heavy2", sumMulti(localNorm2(heavy_prop_2), phases.getSet()));
	}
	pop(xml_out);
	
	
	push(xml_out, "SourceSinkType");
	{
	  QDPIO::cout << "Source_type_1 = " << all_sinks.sink_prop_1.source_type << endl;
	  QDPIO::cout << "Sink_type_1 = " << all_sinks.sink_prop_1.sink_type << endl;
	  QDPIO::cout << "Source_type_2 = " << all_sinks.sink_prop_2.source_type << endl;
	  QDPIO::cout << "Sink_type_2 = " << all_sinks.sink_prop_2.sink_type << endl;
	  if (params.param.FourPt){
	    QDPIO::cout << "Source_type_3 = " << all_sinks.sink_prop_3.source_type << endl;
	    QDPIO::cout << "Sink_type_3 = " << all_sinks.sink_prop_3.sink_type << endl;
	  }
	  write(xml_out, "source_type_1", all_sinks.sink_prop_1.source_type);
	  write(xml_out, "source_disp_type_1", all_sinks.sink_prop_1.source_disp_type);
	  write(xml_out, "sink_type_1", all_sinks.sink_prop_1.sink_type);
	  write(xml_out, "sink_disp_type_1", all_sinks.sink_prop_1.sink_disp_type);
	  
	  write(xml_out, "source_type_2", all_sinks.sink_prop_2.source_type);
	  write(xml_out, "source_disp_type_2", all_sinks.sink_prop_2.source_disp_type);
	  write(xml_out, "sink_type_2", all_sinks.sink_prop_2.sink_type);
	  write(xml_out, "sink_disp_type_2", all_sinks.sink_prop_2.sink_disp_type);
	  if (params.param.FourPt){
	    write(xml_out, "source_type_3", all_sinks.sink_prop_3.source_type);
	    write(xml_out, "source_disp_type_3", all_sinks.sink_prop_3.source_disp_type);
	    write(xml_out, "sink_type_3", all_sinks.sink_prop_3.sink_type);
	    write(xml_out, "sink_disp_type_3", all_sinks.sink_prop_3.sink_disp_type);
	  }
	  write(xml_out, "source_type_h2", all_sinks.heavy_prop_2.source_type);
	  write(xml_out, "source_disp_type_h2", all_sinks.heavy_prop_2.source_disp_type);
	  write(xml_out, "sink_type_h2", all_sinks.heavy_prop_2.sink_type);
	  write(xml_out, "sink_disp_type_h2", all_sinks.heavy_prop_2.sink_disp_type);
	}
	pop(xml_out);
	
	
	// References for use later
	const LatticePropagator& sink_prop_1 = 
	  TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks.sink_prop_1.quark_propagator_id);
	const LatticePropagator& sink_prop_2 = 
	  TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks.sink_prop_2.quark_propagator_id);
	const LatticePropagator& sink_prop_3 = 
	  TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks.sink_prop_3.quark_propagator_id);
	const LatticePropagator& heavy_prop_2 = 
	  TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks.heavy_prop_2.quark_propagator_id);
	
	// Construct group name for output
	string src_type;
	if (all_sinks.sink_prop_1.source_type == "POINT_SOURCE")
	  src_type = "Point";
	else if (all_sinks.sink_prop_1.source_type == "SHELL_SOURCE")
	  src_type = "Shell";
	else if (all_sinks.sink_prop_1.source_type == "WALL_SOURCE")
	  src_type = "Wall";
	else
	  {
	    QDPIO::cerr << "Unsupported source type = " << all_sinks.sink_prop_1.source_type << endl;
	    QDP_abort(1);
	  }
	
	string snk_type;
	if (all_sinks.sink_prop_1.sink_type == "POINT_SINK")
	  snk_type = "Point";
	else if (all_sinks.sink_prop_1.sink_type == "SHELL_SINK")
	  snk_type = "Shell";
	else if (all_sinks.sink_prop_1.sink_type == "WALL_SINK")
	  snk_type = "Wall";
	else
	  {
	    QDPIO::cerr << "Unsupported sink type = " << all_sinks.sink_prop_1.sink_type << endl;
	    QDP_abort(1);
	  }
	
	string source_sink_type = src_type + "_" + snk_type;
	QDPIO::cout << "Source type = " << src_type << endl;
	QDPIO::cout << "Sink type = "   << snk_type << endl;
	
	if (params.param.MesonP) 
	  {
	    Qlbar(u, sink_prop_1, t_src1, phases, xml_out, source_sink_type + "_Wilson_Qlmeson1");
	    if(all_sinks.sink_prop_1.quark_propagator_id!=all_sinks.sink_prop_2.quark_propagator_id)
	      Qlbar(u, sink_prop_2, t_src2, phases, xml_out, source_sink_type + "_Wilson_Qlmeson2");
	    if (params.param.FourPt && all_sinks.sink_prop_3.quark_propagator_id!=all_sinks.sink_prop_2.quark_propagator_id && all_sinks.sink_prop_1.quark_propagator_id!=all_sinks.sink_prop_3.quark_propagator_id)
	      Qlbar(u, sink_prop_3, t_src3, phases, xml_out, source_sink_type + "_Wilson_Qlmeson3");
	  } // end if (MesonP)
	
	//Now we actually will do the contractions...
	
	QlQl(u, sink_prop_1, sink_prop_2, heavy_prop_2,  t_src1, t_src2, t_src2, HQBC, phases, xml_out, source_sink_type +"Wilson_Ql_Ql_3pt");
	
	pop(xml_out);  // array element
      }// End of if statement for first id being static.
      
      else if (named_obj.heavy_id1!="Static" && named_obj.heavy_id2=="Static"){
	// Masses
	write(xml_out, "Mass_1", all_sinks.sink_prop_1.Mass);
	write(xml_out, "Mass_2", all_sinks.sink_prop_2.Mass);
	if (params.param.FourPt)
	  write(xml_out, "Mass_3", all_sinks.sink_prop_3.Mass);
	write(xml_out, "t01", t01);
	write(xml_out, "t02", t02);
	
	// Initialize the slow Fourier transform phases with NO momenta
	SftMom phases(0, true, j_decay);
	
	// Save prop input
	push(xml_out, "Forward_prop_headers");
	write(xml_out, "First_forward_prop", all_sinks.sink_prop_1.prop_header);
	write(xml_out, "Second_forward_prop", all_sinks.sink_prop_2.prop_header);
	if (params.param.FourPt)
	  write(xml_out, "Third_forward_prop", all_sinks.sink_prop_3.prop_header);
	pop(xml_out);
	
	// Sanity check - write out the norm2 of the forward prop in the j_decay direction
	// Use this for any possible verification
	push(xml_out, "Forward_prop_correlator");
	{
	  const LatticePropagator& sink_prop_1 = 
	    TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks.sink_prop_1.quark_propagator_id);
	  const LatticePropagator& sink_prop_2 = 
	    TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks.sink_prop_2.quark_propagator_id);
	  const LatticePropagator& sink_prop_3 = 
	    TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks.sink_prop_3.quark_propagator_id);
	  
	  const LatticePropagator& heavy_prop_1 = 
	    TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks.heavy_prop_1.quark_propagator_id);
	  
	  write(xml_out, "forward_prop_corr_1", sumMulti(localNorm2(sink_prop_1), phases.getSet()));
	  write(xml_out, "forward_prop_corr_2", sumMulti(localNorm2(sink_prop_2), phases.getSet()));
	  
	  write(xml_out, "forward_prop_corr_3", sumMulti(localNorm2(sink_prop_3), phases.getSet()));
	  write(xml_out, "forward_prop_corr_h1", sumMulti(localNorm2(heavy_prop_1), phases.getSet()));
	}
	pop(xml_out);
	    
	push(xml_out, "SourceSinkType");
	{
	  QDPIO::cout << "Source_type_1 = " << all_sinks.sink_prop_1.source_type << endl;
	  QDPIO::cout << "Sink_type_1 = " << all_sinks.sink_prop_1.sink_type << endl;
	  QDPIO::cout << "Source_type_2 = " << all_sinks.sink_prop_2.source_type << endl;
	  QDPIO::cout << "Sink_type_2 = " << all_sinks.sink_prop_2.sink_type << endl;
	  if (params.param.FourPt){
	    QDPIO::cout << "Source_type_3 = " << all_sinks.sink_prop_3.source_type << endl;
	    QDPIO::cout << "Sink_type_3 = " << all_sinks.sink_prop_3.sink_type << endl;
	  }
	  write(xml_out, "source_type_1", all_sinks.sink_prop_1.source_type);
	  write(xml_out, "source_disp_type_1", all_sinks.sink_prop_1.source_disp_type);
	  write(xml_out, "sink_type_1", all_sinks.sink_prop_1.sink_type);
	  write(xml_out, "sink_disp_type_1", all_sinks.sink_prop_1.sink_disp_type);
	  
	  write(xml_out, "source_type_2", all_sinks.sink_prop_2.source_type);
	  write(xml_out, "source_disp_type_2", all_sinks.sink_prop_2.source_disp_type);
	  write(xml_out, "sink_type_2", all_sinks.sink_prop_2.sink_type);
	  write(xml_out, "sink_disp_type_2", all_sinks.sink_prop_2.sink_disp_type);
	  if (params.param.FourPt){
	    write(xml_out, "source_type_3", all_sinks.sink_prop_3.source_type);
	    write(xml_out, "source_disp_type_3", all_sinks.sink_prop_3.source_disp_type);
	    write(xml_out, "sink_type_3", all_sinks.sink_prop_3.sink_type);
	    write(xml_out, "sink_disp_type_3", all_sinks.sink_prop_3.sink_disp_type);
	  }
	  write(xml_out, "source_type_h1", all_sinks.heavy_prop_1.source_type);
	  write(xml_out, "source_disp_type_h1", all_sinks.heavy_prop_1.source_disp_type);
	  write(xml_out, "sink_type_h1", all_sinks.heavy_prop_1.sink_type);
	  write(xml_out, "sink_disp_type_h1", all_sinks.heavy_prop_1.sink_disp_type);
	}
	pop(xml_out);
	
	// References for use later
	const LatticePropagator& sink_prop_1 = 
	  TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks.sink_prop_1.quark_propagator_id);
	const LatticePropagator& sink_prop_2 = 
	  TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks.sink_prop_2.quark_propagator_id);
	const LatticePropagator& sink_prop_3 = 
	  TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks.sink_prop_3.quark_propagator_id);
	
	const LatticePropagator& heavy_prop_1 = 
	  TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks.heavy_prop_1.quark_propagator_id);
	
	// Construct group name for output
	string src_type;
	if (all_sinks.sink_prop_1.source_type == "POINT_SOURCE")
	  src_type = "Point";
	else if (all_sinks.sink_prop_1.source_type == "SHELL_SOURCE")
	  src_type = "Shell";
	else if (all_sinks.sink_prop_1.source_type == "WALL_SOURCE")
	  src_type = "Wall";
	else
	  {
	    QDPIO::cerr << "Unsupported source type = " << all_sinks.sink_prop_1.source_type << endl;
	    QDP_abort(1);
	  }
	
	string snk_type;
	if (all_sinks.sink_prop_1.sink_type == "POINT_SINK")
	  snk_type = "Point";
	else if (all_sinks.sink_prop_1.sink_type == "SHELL_SINK")
	  snk_type = "Shell";
	else if (all_sinks.sink_prop_1.sink_type == "WALL_SINK")
	  snk_type = "Wall";
	else
	  {
	    QDPIO::cerr << "Unsupported sink type = " << all_sinks.sink_prop_1.sink_type << endl;
	    QDP_abort(1);
	  }
	
	string source_sink_type = src_type + "_" + snk_type;
	QDPIO::cout << "Source type = " << src_type << endl;
	QDPIO::cout << "Sink type = "   << snk_type << endl;
	
	if (params.param.MesonP) 
	  {
	    Qlbar(u, sink_prop_1, t_src1, phases, xml_out, source_sink_type + "_Wilson_Qlmeson1");
	    if(all_sinks.sink_prop_1.quark_propagator_id!=all_sinks.sink_prop_2.quark_propagator_id)
	      Qlbar(u, sink_prop_2, t_src2, phases, xml_out, source_sink_type + "_Wilson_Qlmeson2");
	    if (params.param.FourPt && all_sinks.sink_prop_3.quark_propagator_id!=all_sinks.sink_prop_2.quark_propagator_id && all_sinks.sink_prop_1.quark_propagator_id!=all_sinks.sink_prop_3.quark_propagator_id)
	      Qlbar(u, sink_prop_3, t_src3, phases, xml_out, source_sink_type + "_Wilson_Qlmeson3");
	  } // end if (MesonP)
	
	//Now we actually will do the contractions...
	
	QlQl(u, sink_prop_1, sink_prop_2, heavy_prop_1,  t_src1, t_src2, t_src1, HQBC, phases, xml_out, source_sink_type +"Wilson_Ql_Ql_3pt");
	
	pop(xml_out);  // array element
      }// End of if statement for second id being static.
      
      else if (named_obj.heavy_id1!="Static" && named_obj.heavy_id2!="Static"){
	// Masses
	write(xml_out, "Mass_1", all_sinks.sink_prop_1.Mass);
	write(xml_out, "Mass_2", all_sinks.sink_prop_2.Mass);
	if (params.param.FourPt)
	  write(xml_out, "Mass_3", all_sinks.sink_prop_3.Mass);
	write(xml_out, "t01", t01);
	write(xml_out, "t02", t02);
	
	// Initialize the slow Fourier transform phases with NO momenta
	SftMom phases(0, true, j_decay);
	
	// Save prop input
	push(xml_out, "Forward_prop_headers");
	write(xml_out, "First_forward_prop", all_sinks.sink_prop_1.prop_header);
	write(xml_out, "Second_forward_prop", all_sinks.sink_prop_2.prop_header);
	if (params.param.FourPt)
	  write(xml_out, "Third_forward_prop", all_sinks.sink_prop_3.prop_header);
	pop(xml_out);
	
	// Sanity check - write out the norm2 of the forward prop in the j_decay direction
	// Use this for any possible verification
	push(xml_out, "Forward_prop_correlator");
	{
	  const LatticePropagator& sink_prop_1 = 
	    TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks.sink_prop_1.quark_propagator_id);
	  const LatticePropagator& sink_prop_2 = 
	    TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks.sink_prop_2.quark_propagator_id);
	  const LatticePropagator& sink_prop_3 = 
	    TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks.sink_prop_3.quark_propagator_id);
	  
	  const LatticePropagator& heavy_prop_1 = 
	    TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks.heavy_prop_1.quark_propagator_id);
	  const LatticePropagator& heavy_prop_2 = 
	    TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks.heavy_prop_2.quark_propagator_id);
	  
	  write(xml_out, "forward_prop_corr_1", sumMulti(localNorm2(sink_prop_1), phases.getSet()));
	  write(xml_out, "forward_prop_corr_2", sumMulti(localNorm2(sink_prop_2), phases.getSet()));
	  
	  write(xml_out, "forward_prop_corr_3", sumMulti(localNorm2(sink_prop_3), phases.getSet()));
	  write(xml_out, "forward_prop_corr_h1", sumMulti(localNorm2(heavy_prop_1), phases.getSet()));
	  write(xml_out, "forward_prop_corr_h2", sumMulti(localNorm2(heavy_prop_2), phases.getSet()));
	  
	}
	pop(xml_out);
	
	
	push(xml_out, "SourceSinkType");
	{
	  QDPIO::cout << "Source_type_1 = " << all_sinks.sink_prop_1.source_type << endl;
	  QDPIO::cout << "Sink_type_1 = " << all_sinks.sink_prop_1.sink_type << endl;
	  QDPIO::cout << "Source_type_2 = " << all_sinks.sink_prop_2.source_type << endl;
	  QDPIO::cout << "Sink_type_2 = " << all_sinks.sink_prop_2.sink_type << endl;
	  if (params.param.FourPt){
	    QDPIO::cout << "Source_type_3 = " << all_sinks.sink_prop_3.source_type << endl;
	    QDPIO::cout << "Sink_type_3 = " << all_sinks.sink_prop_3.sink_type << endl;
	  }
	  write(xml_out, "source_type_1", all_sinks.sink_prop_1.source_type);
	  write(xml_out, "source_disp_type_1", all_sinks.sink_prop_1.source_disp_type);
	  write(xml_out, "sink_type_1", all_sinks.sink_prop_1.sink_type);
	  write(xml_out, "sink_disp_type_1", all_sinks.sink_prop_1.sink_disp_type);
	  
	  write(xml_out, "source_type_2", all_sinks.sink_prop_2.source_type);
	  write(xml_out, "source_disp_type_2", all_sinks.sink_prop_2.source_disp_type);
	  write(xml_out, "sink_type_2", all_sinks.sink_prop_2.sink_type);
	  write(xml_out, "sink_disp_type_2", all_sinks.sink_prop_2.sink_disp_type);
	  if (params.param.FourPt){
	    write(xml_out, "source_type_3", all_sinks.sink_prop_3.source_type);
	    write(xml_out, "source_disp_type_3", all_sinks.sink_prop_3.source_disp_type);
	    write(xml_out, "sink_type_3", all_sinks.sink_prop_3.sink_type);
	    write(xml_out, "sink_disp_type_3", all_sinks.sink_prop_3.sink_disp_type);
	  }
	  write(xml_out, "source_type_h1", all_sinks.heavy_prop_1.source_type);
	  write(xml_out, "source_disp_type_h1", all_sinks.heavy_prop_1.source_disp_type);
	  write(xml_out, "sink_type_h1", all_sinks.heavy_prop_1.sink_type);
	  write(xml_out, "sink_disp_type_h1", all_sinks.heavy_prop_1.sink_disp_type);
	  
	  write(xml_out, "source_type_h2", all_sinks.heavy_prop_2.source_type);
	  write(xml_out, "source_disp_type_h2", all_sinks.heavy_prop_2.source_disp_type);
	  write(xml_out, "sink_type_h2", all_sinks.heavy_prop_2.sink_type);
	  write(xml_out, "sink_disp_type_h2", all_sinks.heavy_prop_2.sink_disp_type);
	  
	}
	pop(xml_out);
	
	
	// References for use later
	const LatticePropagator& sink_prop_1 = 
	  TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks.sink_prop_1.quark_propagator_id);
	const LatticePropagator& sink_prop_2 = 
	  TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks.sink_prop_2.quark_propagator_id);
	const LatticePropagator& sink_prop_3 = 
	  TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks.sink_prop_3.quark_propagator_id);
	
	const LatticePropagator& heavy_prop_1 = 
	  TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks.heavy_prop_1.quark_propagator_id);
	const LatticePropagator& heavy_prop_2 = 
	  TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks.heavy_prop_2.quark_propagator_id);
	// Construct group name for output
	string src_type;
	if (all_sinks.sink_prop_1.source_type == "POINT_SOURCE")
	  src_type = "Point";
	else if (all_sinks.sink_prop_1.source_type == "SHELL_SOURCE")
	  src_type = "Shell";
	else if (all_sinks.sink_prop_1.source_type == "WALL_SOURCE")
	  src_type = "Wall";
	else
	  {
	    QDPIO::cerr << "Unsupported source type = " << all_sinks.sink_prop_1.source_type << endl;
	    QDP_abort(1);
	  }
	
	string snk_type;
	if (all_sinks.sink_prop_1.sink_type == "POINT_SINK")
	  snk_type = "Point";
	else if (all_sinks.sink_prop_1.sink_type == "SHELL_SINK")
	  snk_type = "Shell";
	else if (all_sinks.sink_prop_1.sink_type == "WALL_SINK")
	  snk_type = "Wall";
	else
	  {
	    QDPIO::cerr << "Unsupported sink type = " << all_sinks.sink_prop_1.sink_type << endl;
	    QDP_abort(1);
	  }
	
	string source_sink_type = src_type + "_" + snk_type;
	QDPIO::cout << "Source type = " << src_type << endl;
	QDPIO::cout << "Sink type = "   << snk_type << endl;
	
	QDPIO::cout << "Note that we are not calculating Heavy-Light Spectrum,"<<endl;
	QDPIO::cout << "because both heavy quarks are not static."<<endl;
	
	//Now we actually will do the contractions...
	
	QlQl(u, sink_prop_1, sink_prop_2, heavy_prop_1, heavy_prop_2,  t_src1, t_src2, phases, xml_out, source_sink_type +"Wilson_Ql_Ql_3pt");
	
	pop(xml_out);  // array element
      }// End of if statement for no statics.
      
    }// End of loop over pairs...
    pop(xml_out);  // Wilson_spectroscopy
    pop(xml_out);  // HeavyLightCont
    
    snoop.stop();
    QDPIO::cout << InlineHeavyLightContEnv::name << ": total time = "
		<< snoop.getTimeInSeconds() 
		<< " secs" << endl;

    QDPIO::cout << InlineHeavyLightContEnv::name << ": ran successfully" << endl;

    END_CODE();
  } 

};
