// $Id: inline_npr_vertex_w.cc,v 1.1 2006-10-30 22:31:43 edwards Exp $
/*! \file
 * \brief Inline construction of NPR vertices
 *
 * NPR vertices on  props
 */

#include "meas/inline/hadron/inline_npr_vertex_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "util/ft/sftmom.h"
#include "meas/hadron/BuildingBlocks_w.h"
#include "util/info/proginfo.h"
#include "meas/inline/make_xml_file.h"

#include "meas/inline/io/named_objmap.h"

#include "actions/ferm/fermstates/ferm_createstate_factory_w.h"
#include "actions/ferm/fermstates/ferm_createstate_aggregate_w.h"

namespace Chroma 
{ 
  namespace InlineNprVertexEnv 
  { 
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineNprVertex(InlineNprVertexParams(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "BUILDING_BLOCKS";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= CreateFermStateEnv::registerAll();
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
    }
  }


  //! Param input
  void read(XMLReader& xml, const string& path, InlineNprVertexParams::Param_t& input)
  {
    XMLReader paramtop(xml, path);

    int version;
    read(paramtop, "version", version);

    if (paramtop.count("FermState") != 0)
      input.cfs = readXMLGroup(paramtop, "FermState", "Name");
    else
      input.cfs = CreateFermStateEnv::nullXMLGroup();

    switch (version) 
    {
    case 1:
      break;

    default :
      QDPIO::cerr << InlineNprVertexEnv::name << ": input parameter version " 
		  << version << " unsupported." << endl;
      QDP_abort(1);
    }
    
    read(paramtop, "links_max", input.links_max);
  }


  //! Param write
  void write(XMLWriter& xml, const string& path, const InlineNprVertexParams::Param_t& input)
  {
    push(xml, path);

    int version = 1;
    write(xml, "version", version);
    write(xml, "links_max", input.links_max);
    xml << input.cfs.xml;

    pop(xml);
  }

  //! Propagator input
  void read(XMLReader& xml, const string& path, InlineNprVertexParams::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "gauge_id", input.gauge_id);
    read(inputtop, "prop_id", input.prop_id);
    read(inputtop, "bb_file_name", input.bb_file_name);
  }

  //! Propagator output
  void write(XMLWriter& xml, const string& path, const InlineNprVertexParams::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "gauge_id", input.gauge_id);
    write(xml, "prop_id", input.prop_id);
    write(xml, "bb_file_name", input.bb_file_name);

    pop(xml);
  }


  // Param stuff
  InlineNprVertexParams::InlineNprVertexParams() {frequency = 0;}

  InlineNprVertexParams::InlineNprVertexParams(XMLReader& xml_in, const std::string& path) 
  {
    try 
    {
      XMLReader paramtop(xml_in, path);

      if (paramtop.count("Frequency") == 1)
	read(paramtop, "Frequency", frequency);
      else
	frequency = 1;

      // Read program parameters
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
  InlineNprVertexParams::write(XMLWriter& xml_out, const std::string& path) 
  {
    push(xml_out, path);
    
    Chroma::write(xml_out, "Param", param); 
    Chroma::write(xml_out, "NamedObject", named_obj);
    QDP::write(xml_out, "xml_file", xml_file);

    pop(xml_out);
  }


  //###################################################################################//
  // Accept All Link Patterns                                                          //
  //###################################################################################//

  void AllLinkPatterns( bool &                          DoThisPattern,
			bool &                          DoFurtherPatterns,
			multi1d< unsigned short int > & LinkPattern )
  {
    DoThisPattern     = true;
    DoFurtherPatterns = true;
    
    return;
  }


  // Function call
  void 
  InlineNprVertex::operator()(unsigned long update_no,
				   XMLWriter& xml_out) 
  {
    // If xml file not empty, then use alternate
    if (params.xml_file != "")
    {
      string xml_file = makeXMLFileName(params.xml_file, update_no);

      push(xml_out, "NprVertex");
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


  // Function call
  void 
  InlineNprVertex::func(unsigned long update_no,
			     XMLWriter& XmlOut) 
  {
    START_CODE();

    StopWatch snoop;
    snoop.reset();
    snoop.start();

    push(XmlOut, "NprVertex");
    write(XmlOut, "update_no", update_no);

    QDPIO::cout << " ExampleNprVertex" << endl;
    QDPIO::cout << "     volume: " << QDP::Layout::lattSize()[0];
    for (int i=1; i<Nd; ++i) {
      QDPIO::cout << " x " << QDP::Layout::lattSize()[i];
    }
    QDPIO::cout << endl;

    //#################################################################################//
    // XML output
    //#################################################################################//

    proginfo(XmlOut);    // Print out basic program info

    push(XmlOut, "Output_version");
    write(XmlOut, "out_version", 2);
    pop(XmlOut);

    //###############################################################################//
    // Read Gauge Field                                                              //
    //###############################################################################//

    Out << "Attempt to initialize the gauge field" << "\n";  Out.flush();

    // Grab the gauge field
    multi1d<LatticeColorMatrix> U;
    XMLBufferWriter gauge_xml;

    try
    {
      U = TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.bb.GaugeId);
      TheNamedObjMap::Instance().get(params.bb.GaugeId).getRecordXML(gauge_xml);

      // Set the construct state and modify the fields
      {
QDPIO::cout << "cfs=XX" << params.param.cfs.xml << "XX" << endl;
	std::istringstream  xml_s(params.param.cfs.xml);
	XMLReader  fermtop(xml_s);

	Handle<CreateFermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > 
	  cfs(TheCreateFermStateFactory::Instance().createObject(params.param.cfs.id,
								 fermtop,
								 params.param.cfs.path));

	Handle<FermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > 
	  state((*cfs)(U));

	// Pull the u fields back out from the state since they might have been
	// munged with fermBC's
	U = state->getLinks();
      }
    
    }
    catch( std::bad_cast ) 
    {
      QDPIO::cerr << InlineNprVertexEnv::name << ": caught dynamic cast error" 
		  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << InlineNprVertexEnv::name << ": map call failed: " << e 
		  << endl;
      QDP_abort(1);
    }
    catch( ... )
    {
      QDPIO::cerr << InlineNprVertexEnv::name << ": caught generic exception "
		  << endl;
      QDP_abort(1);
    }

    // Write out the input
    params.write(XmlOut, "Input");

    // Write out the config info
    write(XmlOut, "Config_info", gauge_xml);

    // check that the gauge field seems normal
    Double ave_plaq, ave_spacelike_plaq, ave_timelike_plaq, ave_link_trace;
    MesPlq( U, ave_plaq, ave_spacelike_plaq, ave_timelike_plaq, ave_link_trace );
    Out << "basic gauge field observables"                         << "\n";
    Out << "average plaquette            = " << ave_plaq           << "\n";
    Out << "average space-like plaquette = " << ave_spacelike_plaq << "\n";
    Out << "average time-like plaquette  = " << ave_timelike_plaq  << "\n";
    Out << "average link trace           = " << ave_link_trace     << "\n";

    push(XmlOut, "Observables");
    write(XmlOut, "ave_plaq", ave_plaq);
    write(XmlOut, "ave_spacelike_plaq", ave_spacelike_plaq);
    write(XmlOut, "ave_timelike_plaq", ave_timelike_plaq);
    write(XmlOut, "ave_link_trace", ave_link_trace);
    pop(XmlOut);

    //#################################################################################//
    // Read Forward Propagator                                                         //
    //#################################################################################//

    SftMom phases_nomom( 0, true, Nd-1 );  // used to check props. Fix to Nd-1 direction.

    LatticePropagator F;
    ChromaProp_t prop_header;
    PropSourceConst_t source_header;
    QDPIO::cout << "Attempt to parse forward propagator" << endl;
    Out << "parsing forward propagator " << params.bb.FrwdPropId << " ... " << "\n";  Out.flush();

    try
    {
      // Snarf a copy
      F = TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.prop_id);
	
      // Snarf the frwd prop info. This is will throw if the frwd prop id is not there
      XMLReader PropXML;
      TheNamedObjMap::Instance().get(params.bb.FrwdPropId).getFileXML(PropXML);
      TheNamedObjMap::Instance().get(params.bb.FrwdPropId).getRecordXML(PropRecordXML);

      // Try to invert this record XML into a ChromaProp struct
      {
	read(PropRecordXML, "/Propagator/ForwardProp", prop_header);
	read(PropRecordXML, "/Propagator/PropSource", source_header);
      }

      // Sanity check - write out the norm2 of the forward prop in the j_decay direction
      // Use this for any possible verification
      {
	multi1d<Double> FrwdPropCheck = 
	  sumMulti( localNorm2( F ), phases_nomom.getSet() );

	Out << "forward propagator check = " << PropCheck[0] << "\n";  Out.flush();

	// Write out the forward propagator header
	push(XmlOut, "ForwardProp");
	write(XmlOut, "FrwdPropId", params.bb.FrwdPropId);
	write(XmlOut, "FrwdPropXML", FrwdPropXML);
	write(XmlOut, "FrwdPropRecordXML", FrwdPropRecordXML);
	write(XmlOut, "FrwdPropCheck", FrwdPropCheck);
	pop(XmlOut);
      }
    }
    catch( std::bad_cast ) 
    {
      QDPIO::cerr << InlineNprVertexEnv::name << ": caught dynamic cast error" 
		  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << InlineNprVertexEnv::name << ": forward prop: error message: " << e 
		  << endl;
      QDP_abort(1);
    }

    QDPIO::cout << "Forward propagator successfully parsed" << endl;
    Out << "finished parsing forward propagator " << params.bb.FrwdPropId << "\n";  Out.flush();


    // Derived from input prop
    int  j_decay = source_header.j_decay;
    multi1d<int> t_srce = source_header.getTSrce() ;


    //#################################################################################//
    // Construct Building Blocks                                                       //
    //#################################################################################//
    
    swatch.reset();
    Out << "calculating building blocks" << "\n";  Out.flush();
    QDPIO::cout << "calculating building blocks" << endl;

    const signed short int T1 = 0;
    const signed short int T2 = QDP::Layout::lattSize()[j_decay] - 1;
    const signed short int DecayDir = j_decay;
    const signed short int Tsrc = source_header.t_source;
    const signed short int Tsnk = seqsource_header.t_sink;

    swatch.start();
    NprVertex(B, F, U, 
	      GammaInsertions, Flavors,
	      params.param.links_max, AllLinkPatterns, 
	      Phases, PhasesCanonical,
	      Files, T1, T2,
	      Tsrc, Tsnk,
	      seqsource_header.seqsrc.id, seqsource_header.sink_mom, DecayDir,
	      params.param.time_reverse,
	      params.param.translate);
    swatch.stop();
      
    QDPIO::cout << "finished calculating building blocks for loop = " << loop 
		<< "  time= "
		<< swatch.getTimeInSeconds() 
		<< " secs" << endl;

    pop(XmlOut);   // NprVertex

    snoop.stop();
    QDPIO::cout << InlineNprVertexEnv::name << ": total time = "
		<< snoop.getTimeInSeconds() 
		<< " secs" << endl;

    QDPIO::cout << InlineNprVertexEnv::name << ": ran successfully" << endl;

    END_CODE();
  } 

}
