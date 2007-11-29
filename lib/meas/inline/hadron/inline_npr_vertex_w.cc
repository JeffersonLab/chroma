// $Id: inline_npr_vertex_w.cc,v 1.7 2007-11-29 04:13:07 kostas Exp $
/*! \file
 * \brief Inline construction of NPR vertices
 *
 * NPR vertices on  props
 */

#include "meas/inline/hadron/inline_npr_vertex_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "util/ft/sftmom.h"
// #include "meas/hadron/BuildingBlocks_w.h"
#include "util/info/proginfo.h"
#include "meas/inline/make_xml_file.h"

#include "meas/hadron/npr_vertex_w.h"

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

    const std::string name = "NPR_VERTEX";

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

  DPropagator FTpropagator(const LatticePropagator& prop,
			   const multi1d<int> mom,
			   const multi1d<int> t_src){
    // Initialize the slow Fourier transform phases
    // Makes no sence to do momemtum average here...
    SftMom phases(0, t_src, mom);
    multi1d<int> mm = phases.numToMom(0);
    QDPIO::cout<<"Sink momentum: " ;
    QDPIO::cout<<mm[0]<<mm[1]<<mm[2]<<mm[3]<<endl;
    return sum(prop*phases[0])/Double(Layout::vol());
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
    read(paramtop, "file_name", input.file_name);
  }


  //! Param write
  void write(XMLWriter& xml, const string& path, const InlineNprVertexParams::Param_t& input)
  {
    push(xml, path);

    int version = 1;
    write(xml, "version", version);
    write(xml, "links_max", input.links_max);
    write(xml, "file_name", input.file_name);    
    xml << input.cfs.xml;

    pop(xml);
  }

  //! Propagator input
  void read(XMLReader& xml, const string& path, InlineNprVertexParams::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "gauge_id", input.gauge_id);
    read(inputtop, "prop_id", input.prop_id);
  }

  //! Propagator output
  void write(XMLWriter& xml, const string& path, const InlineNprVertexParams::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "gauge_id", input.gauge_id);
    write(xml, "prop_id", input.prop_id);

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
			multi1d< int > & LinkPattern )
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

    QDPIO::cout << "Attempt to initialize the gauge field" << endl << flush ;

    // Grab the gauge field
    multi1d<LatticeColorMatrix> U;
    XMLBufferWriter gauge_xml;

    try
    {
      U = TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);
      TheNamedObjMap::Instance().get(params.named_obj.gauge_id).getRecordXML(gauge_xml);

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

    // Calculate some gauge invariant observables just for info.
    MesPlq(XmlOut, "Observables", U);

    //#################################################################################//
    // Read Forward Propagator                                                         //
    //#################################################################################//

    SftMom phases_nomom( 0, true, Nd-1 );  // used to check props. Fix to Nd-1 direction.

    LatticePropagator F;
    ChromaProp_t prop_header;
    PropSourceConst_t source_header;
    QDPIO::cout << "Attempt to parse forward propagator" << endl;
    QDPIO::cout << "parsing forward propagator " << params.named_obj.prop_id << " ... " << endl << flush;

    try
    {
      // Snarf a copy
      F = TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.prop_id);
	
      // Snarf the frwd prop info. This is will throw if the frwd prop id is not there
      XMLReader PropXML, PropRecordXML;
      TheNamedObjMap::Instance().get(params.named_obj.prop_id).getFileXML(PropXML);
      TheNamedObjMap::Instance().get(params.named_obj.prop_id).getRecordXML(PropRecordXML);

      // Try to invert this record XML into a ChromaProp struct
      {
	read(PropRecordXML, "/Propagator/ForwardProp", prop_header);
	read(PropRecordXML, "/Propagator/PropSource", source_header);
      }

      // Sanity check - write out the norm2 of the forward prop in the j_decay direction
      // Use this for any possible verification
      {
	multi1d<Double> PropCheck = 
	  sumMulti( localNorm2( F ), phases_nomom.getSet() );

	QDPIO::cout << "forward propagator check = " << PropCheck[0] << endl;

	// Write out the forward propagator header
	push(XmlOut, "ForwardProp");
	write(XmlOut, "PropXML", PropXML);
	write(XmlOut, "PropRecordXML", PropRecordXML);
	write(XmlOut, "PropCheck", PropCheck);
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


    // Get the momentum from the header

    multi1d<int> mom  ;
    multi1d<int> t_src ;
    int       t_dir = source_header.j_decay ;

    try{
      mom = source_header.getMom() ;
      t_src = source_header.getTSrce() ;
    }
    catch (const string& e){
      QDPIO::cerr << InlineNprVertexEnv::name << ": propagator does not have a momentum source or t_src not present: error message: " << e << endl;
      QDP_abort(1);
    }


    //#################################################################################//
    // Construct Building Blocks                                                       //
    //#################################################################################//
    QDP::StopWatch swatch;
    swatch.reset();
    
    XMLBufferWriter file_xml;
    push(file_xml, "NprVertex");
    write(file_xml, "Param", params.param);
    write(file_xml, "ForwardProp", prop_header);
    write(file_xml, "PropSource", source_header);
    write(file_xml, "Config", gauge_xml);
    pop(file_xml);

    QDPFileWriter qio_file(file_xml, params.param.file_name,QDPIO_SINGLEFILE, 
			   QDPIO_SERIAL, QDPIO_OPEN); 


    //Fourier transform the propagator 
    QDPIO::cout << "Fourier Transforming propagator" << endl;
    swatch.start();
    multi1d<int> neg_mom(mom.size());
    //need the oposit momentum on the sink
    for(int i(0);i<mom.size();i++)
      neg_mom[i] = -mom[i] ;
    DPropagator FTprop(FTpropagator(F,neg_mom,t_src));
    swatch.stop();
    XMLBufferWriter prop_xml;
    push(prop_xml,"QuarkPropagator");
    write(prop_xml,"mom",mom);
    write(prop_xml,"origin",t_src);
    write(prop_xml,"t_dir",t_dir);
    pop(prop_xml) ;
    write(qio_file, prop_xml, FTprop);

    QDPIO::cout << "finished Fourier Transforming propagator"
                << "  time= "
                << swatch.getTimeInSeconds()
                << " secs" << endl;

    QDPIO::cout << "Calculating building blocks" << endl;
    swatch.reset();
    swatch.start();
    NprVertex(F, U, params.param.links_max, AllLinkPatterns, qio_file);
    swatch.stop();
      
    close(qio_file);

    QDPIO::cout << "finished calculating NprVertex"
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
