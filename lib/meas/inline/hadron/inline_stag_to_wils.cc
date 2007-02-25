// $Id: inline_stag_to_wils.cc,v 3.5 2007-02-25 22:39:29 edwards Exp $
/*! \file
 * \brief Inline construction of propagator
 *
 * Propagator calculations
 */

#include "fermact.h"
#include "util/ferm/transf.h"
#include "actions/ferm/fermacts/asqtad_fermact_s.h"
#include "meas/inline/hadron/inline_stag_to_wils.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "util/info/unique_id.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"
#include "meas/inline/make_xml_file.h"

#include "meas/inline/io/named_objmap.h"


namespace Chroma 
{ 
  namespace InlineStagToWilsEnv 
  { 
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineStagToWils(InlineStagToWilsParams(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "STAG_TO_WILS";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= WilsonTypeFermActsEnv::registerAll();
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
    }
  } // end namespace


  //! StagToWils input
  void read(XMLReader& xml, const string& path, InlineStagToWilsParams::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "wils_id", input.wils_id);
    read(inputtop, "stag_id", input.stag_id);
  }

  //! StagToWils output
  void write(XMLWriter& xml, const string& path, const InlineStagToWilsParams::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "stag_id", input.stag_id);
    write(xml, "wils_id", input.wils_id);

    pop(xml);
  }


  // Param stuff
  InlineStagToWilsParams::InlineStagToWilsParams() { frequency = 0; }

  InlineStagToWilsParams::InlineStagToWilsParams(XMLReader& xml_in, const std::string& path) 
  {
    try 
    {
      XMLReader paramtop(xml_in, path);

      if (paramtop.count("Frequency") == 1)
	read(paramtop, "Frequency", frequency);
      else
	frequency = 1;

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
  InlineStagToWilsParams::writeXML(XMLWriter& xml_out, const std::string& path) 
  {
    push(xml_out, path);
    
    write(xml_out, "NamedObject", named_obj);

    pop(xml_out);
  }


  // Function call
  void 
  InlineStagToWils::operator()(unsigned long update_no,
			       XMLWriter& xml_out) 
  {
    // If xml file not empty, then use alternate
    if (params.xml_file != "")
    {
      string xml_file = makeXMLFileName(params.xml_file, update_no);

      push(xml_out, "stag_to_wils");
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
  InlineStagToWils::func(unsigned long update_no,
			 XMLWriter& xml_out) 
  {
    START_CODE();

    StopWatch snoop;
    snoop.reset();
    snoop.start();

    push(xml_out, "stag_to_wils");
    write(xml_out, "update_no", update_no);

    QDPIO::cout << InlineStagToWilsEnv::name << ": propagator conversion" << endl;

    proginfo(xml_out);    // Print out basic program info

    // Write out the input
    params.writeXML(xml_out, "Input");

    push(xml_out, "Output_version");
    write(xml_out, "out_version", 1);
    pop(xml_out);

    //
    // Read in the source along with relevant information.
    // 
    XMLReader stag_file_xml, stag_record_xml;

    multi1d<int> t_srce  ;
    int t0;
    int j_decay;
    bool make_sourceP = false;
    bool SmearedSink = false;

    QDPIO::cout << "Snarf the staggered propagator from a named buffer" << endl;
    try
    {
      // Try the cast to see if this is a valid propagator
      LatticeStaggeredPropagator& stag_tmp = 
	TheNamedObjMap::Instance().getData<LatticeStaggeredPropagator>(params.named_obj.stag_id);

      // Snarf the propagator info. 
      // This is will throw if the stag_id is not there
      TheNamedObjMap::Instance().get(params.named_obj.stag_id).getFileXML(stag_file_xml);
      TheNamedObjMap::Instance().get(params.named_obj.stag_id).getRecordXML(stag_record_xml);

      // Try to invert this record XML into a source struct
      // First identify what kind of source might be here
      if (stag_record_xml.count("/Propagator") != 0)
      {
	PropSourceConst_t source_header;

	read(stag_record_xml, "/Propagator/PropSource", source_header);
	j_decay = source_header.j_decay;
	t0 = source_header.t_source;
	t_srce = source_header.getTSrce() ;
	//Need to figure out what you do for Wall noise, or Corner sources
	//plane wall or noise source should not be supported.
	// corner sources do not exist...
	// All this should be fixed in the future
	make_sourceP = true;
      }
      // When smeared sink propagators are found 
      // I need to do something different
      else if (stag_record_xml.count("/SinkSmear") != 0)
      {
	SmearedSink=true ;
	PropSourceConst_t source_header;

	read(stag_record_xml, "/SinkSmear/PropSource", source_header);
	j_decay = source_header.j_decay;
	t0 = source_header.t_source;
	t_srce = source_header.getTSrce() ;
	//Need to figure out what you do for Wall noise, or Corner sources
	//plane wall or noise source should not be supported.
	// corner sources do not exist...
	// All this should be fixed in the future
	make_sourceP = true;
      }
      
      else
      {
	throw std::string("No appropriate header found");
      }

      // Write out the source header
      write(xml_out, "StaggeredProp_file_info", stag_file_xml);
      write(xml_out, "StaggeredProp_record_info", stag_record_xml);
    }    
    catch (std::bad_cast)
    {
      QDPIO::cerr << InlineStagToWilsEnv::name << ": caught dynamic cast error" 
		  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << InlineStagToWilsEnv::name << ": error extracting source_header: " << e << endl;
      QDP_abort(1);
    }

    // Should be a valid cast now
    const LatticeStaggeredPropagator& stag_prop = 
      TheNamedObjMap::Instance().getData<LatticeStaggeredPropagator>(params.named_obj.stag_id);
 
    QDPIO::cout << "Staggered propagator successfully read and parsed" << endl;

    // Sanity check - write out the norm2 of the source in the Nd-1 direction
    // Use this for any possible verification
    {
      // Initialize the slow Fourier transform phases
      SftMom phases(0, true, Nd-1);

      multi1d<Double> stag_corr = sumMulti(localNorm2(stag_prop), 
					     phases.getSet());

      push(xml_out, "StaggeredProp_correlator");
      write(xml_out, "stag_prop_corr", stag_corr);
      pop(xml_out);
    }

    //
    // Get the Wilson propagator
    //
    try
    {
      TheNamedObjMap::Instance().create<LatticePropagator>(params.named_obj.wils_id);
    }
    catch (std::bad_cast)
    {
      QDPIO::cerr << InlineStagToWilsEnv::name << ": caught dynamic cast error" 
		  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << InlineStagToWilsEnv::name << ": error creating wils_prop: " << e << endl;
      QDP_abort(1);
    }

    // Cast should be valid now
    LatticePropagator& wils_prop = 
      TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.wils_id);

    StopWatch swatch;
    

    if (!make_sourceP){
      QDPIO::cerr<<"Error: Only forward props supported" << endl;
      QDP_abort(1);
    }
	
    LatticeColorMatrix xx = peekSpin(stag_prop,0,0) ;
    wils_prop = zero ;
    for(int s(0);s<Ns;s++){
      pokeSpin(wils_prop,xx,s,s);
    }//diagonal in spin

    swatch.start() ;
    //Multiply the source by Gamma(1)^x Gamma(2)^y Gamma(4)^z Gamma(8)^t 
    QDPIO::cout<<"SOURCE: " ;
    for(int d(0);d<Nd;d++)
      QDPIO::cout<<t_srce[d]<<" ";
    QDPIO::cout<<endl;
    SpinMatrix Omega = 1.0 ;
    int g(0);
    for(int mu(0);mu<Nd;mu++)
      g = g | ((t_srce[mu]%2)<<mu) ;
    QDPIO::cout<<"SOURCE gamma: "<<g<<endl ;
    Omega = Gamma(g)*Omega ;
    
    //now the sink location
    LatticeSpinMatrix lOmega=1.0 ;
    LatticeInteger lg = zero;
    for(int mu(0);mu<Nd;mu++){
      lg = lg | ((Layout::latticeCoordinate(mu)%2)<<mu) ;
    }
    for(int g(0);g<(Ns*Ns-1);g++)//loop over all possible gamma matrices
      lOmega = where((lg==g),Gamma(g)*lOmega,lOmega) ;

    wils_prop = lOmega*wils_prop*adj(Omega) ;

    swatch.stop();
    QDPIO::cout << "Staggered Propagator coverted to Wilson: time= " 
		<< swatch.getTimeInSeconds() 
		<< " secs" << endl;

 
    // Sanity check - write out the propagator (pion) correlator in the Nd-1 direction
    {
      // Initialize the slow Fourier transform phases
      SftMom phases(0, true, Nd-1);

      multi1d<Double> prop_corr=sumMulti(localNorm2(wils_prop),phases.getSet());

      push(xml_out, "WilsonProp_correlator");
      write(xml_out, "wils_corr", prop_corr);
      pop(xml_out);
    }

    // Save the propagator info
    try
    {
      QDPIO::cout << "Start writing propagator info" << endl;

      XMLBufferWriter file_xml;
      push(file_xml, "propagator");
      write(file_xml, "id", uniqueId());  // NOTE: new ID form
      pop(file_xml);

      XMLBufferWriter record_xml;
      if (make_sourceP)
      {
	if(SmearedSink){
	  XMLReader xml_tmp(stag_record_xml, "/SinkSmear");
	  
	  push(record_xml, "SinkSmear");
	  //write(record_xml, "ForwardProp", params.param);
	  //record_xml << params.stateInfo;  // write out the stateinfo - might be empty
	  record_xml << xml_tmp;  // write out all the stuff under MakeSource
	  pop(record_xml);
	}
	else{
	  XMLReader xml_tmp(stag_record_xml, "/Propagator");
	  
	  push(record_xml, "Propagator");
	  //write(record_xml, "ForwardProp", params.param);
	  //record_xml << params.stateInfo;  // write out the stateinfo - might be empty
	  record_xml << xml_tmp;  // write out all the stuff under MakeSource
	  pop(record_xml);
	}
      } 

      // Write the propagator xml info
      TheNamedObjMap::Instance().get(params.named_obj.wils_id).setFileXML(file_xml);
      TheNamedObjMap::Instance().get(params.named_obj.wils_id).setRecordXML(record_xml);

      QDPIO::cout << "Propagator successfully updated" << endl;
    }
    catch (std::bad_cast)
    {
      QDPIO::cerr << InlineStagToWilsEnv::name << ": caught dynamic cast error" 
		  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << InlineStagToWilsEnv::name << ": error extracting prop_header: " << e << endl;
      QDP_abort(1);
    }

    pop(xml_out);  // propagator

    snoop.stop();
    QDPIO::cout << InlineStagToWilsEnv::name << ": total time = "
		<< snoop.getTimeInSeconds() 
		<< " secs" << endl;

    QDPIO::cout << InlineStagToWilsEnv::name << ": ran successfully" << endl;

    END_CODE();
  } 

}
