// $Id: inline_npr_w.cc,v 1.6 2007-11-16 22:27:33 kostas Exp $
/*! \file
 * \brief Inline construction of NPR propagator
 *
 * Propagator calculations
 */

#include "fermact.h"
#include "meas/inline/hadron/inline_npr_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"
#include "meas/inline/make_xml_file.h"

#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{ 
  namespace InlineNprEnv 
  { 
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineNpr(InlineNprParams(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "NPR";

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


  //! Npr input
  void read(XMLReader& xml, const string& path, InlineNprParams::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "gauge_id", input.gauge_id);
    //read(inputtop, "source_id", input.source_id);
    read(inputtop, "prop_id", input.prop_id);
  }

  //! Npr output
  void write(XMLWriter& xml, const string& path, const InlineNprParams::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "gauge_id", input.gauge_id);
    //write(xml, "source_id", input.source_id);
    write(xml, "prop_id", input.prop_id);

    pop(xml);
  }


  // Param stuff
  InlineNprParams::InlineNprParams() { frequency = 0; }

  InlineNprParams::InlineNprParams(XMLReader& xml_in, const std::string& path) 
  {
    try 
    {
      XMLReader paramtop(xml_in, path);

      if (paramtop.count("Frequency") == 1)
	read(paramtop, "Frequency", frequency);
      else
	frequency = 1;

       // Read in the output npr/source configuration info
      read(paramtop, "max_mom2", max_mom2);
      output_type = "LIME";
      if(paramtop.count("output_type")==1) 
	read(paramtop, "output_type", output_type);

      read(paramtop, "filename", filename);

      // Read in the output npr/source configuration info
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
  InlineNprParams::write(XMLWriter& xml_out, const std::string& path) 
  {
    push(xml_out, path);
    
    QDP::write(xml_out, "filename", filename);
    QDP::write(xml_out, "max_mom2", max_mom2);
    QDP::write(xml_out, "output_type", output_type);
    Chroma::write(xml_out, "NamedObject", named_obj);

    pop(xml_out);
  }

  /**
  void InlineNpr::make_source(LatticePropagator& src,
			      const Handle<const ConnectState>& state,
			      const mult1d<ind>& t_source,
			      int mu){
    multi1d<LatticeColorMatrix>& u = state.getLinks() ;
    for(int color_source = 0; color_source < Nc; ++color_source){
      QDPIO::cout << "color = " << color_source << endl; 
      LatticeColorVector cvec = zero;
      // Make a point source at coordinates t_source
      srcfil(src_color_vec, t_source, color_source);
      if((mu>=0)&&(mu<Nd)){
	LatticeColorVector tt = cvec ;
	cvec=0.5*(u[mu]*shift(tt,FORWARD,mu) - shift(adj(u[mu])*tt,BACKWARD,mu));
      }
      for(int spin_source = 0; spin_source < Ns; ++spin_source){
	QDPIO::cout << "spin = " << spin_source << endl; 
	// Insert a ColorVector into spin index spin_source
	// This only overwrites sections, so need to initialize first
	LatticeFermion chi = zero;
	CvToFerm(cvec, chi, spin_source);
	FermToProp(chi, src, color_source, spin_source);
      }
    }
  }
  **/

  // Function call
  void 
  InlineNpr::operator()(unsigned long update_no,
			       XMLWriter& xml_out) 
  {
    // If xml file not empty, then use alternate
    if (params.xml_file != "")
    {
      string xml_file = makeXMLFileName(params.xml_file, update_no);

      push(xml_out, "npr");
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
  InlineNpr::func(unsigned long update_no, XMLWriter& xml_out) 
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
	QDPIO::cerr << InlineNprEnv::name << ": caught dynamic cast error" 
		    << endl;
	QDP_abort(1);
      }
    catch (const string& e) 
      {
	QDPIO::cerr << InlineNprEnv::name << ": map call failed: " << e 
		    << endl;
	QDP_abort(1);
      }
    const multi1d<LatticeColorMatrix>& u = 
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);
    
    XMLReader prop_file_xml,prop_xml ;
    TheNamedObjMap::Instance().get(params.named_obj.prop_id).getFileXML(prop_file_xml);
    TheNamedObjMap::Instance().get(params.named_obj.prop_id).getRecordXML(prop_xml);
    LatticePropagator quark_propagator = TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.prop_id);

    push(xml_out, "npr");
    write(xml_out, "update_no", update_no);
    XMLBufferWriter file_xml;
    push(file_xml, "NPR_W");
    push(file_xml, "Output_version");
    write(file_xml, "out_version", 1);
    pop(file_xml);
    proginfo(file_xml);    // Print out basic program info
    params.write(file_xml, "Input");

    QDPIO::cout << InlineNprEnv::name << ": npr calculation" << endl;
    
    proginfo(xml_out);    // Print out basic program info
    
    // Write out the input
    params.write(xml_out, "Input");
    
    // Write out the config header
    write(xml_out, "Config_info", gauge_xml);
    
    push(xml_out, "Output_version");
    write(xml_out, "out_version", 1);
    pop(xml_out);
    write(xml_out, "Config_info", gauge_xml);
    write(xml_out, "Propagator_info", prop_xml);

    write(file_xml, "Config_info", gauge_xml);
    write(file_xml, "Propagator_info", prop_xml);
    pop(file_xml) ; //NPR_W

    // Calculate some gauge invariant observables just for info.
    MesPlq(xml_out, "Observables", u);
    
    // Make sure that the source location is irrelevant
    SftMom phases(params.max_mom2, false, Nd); // 4D fourier transform
        	      
    //now need to Fourier transform 
    multi1d<DPropagator> FF(phases.numMom());//= sumMulti(quark_propagator,phases.getSet());
    for (int m(0); m < FF.size(); m++){
      multi1d<DPropagator> tt ;
      tt = sumMulti(phases[m]*quark_propagator, phases.getSet()) ;
      FF[m] =  tt[0] ;
    }
    
    // Either write out prop as text or binary. The options are to set the option
    if(params.output_type=="XML"){
      push(xml_out,"FT_prop") ;
      for(int p(0) ; p<phases.numMom();p++){
	push(xml_out, "prop_desc");//write out the momemtum of each bit
	write(xml_out, "mom", phases.numToMom(p));
	pop(xml_out);//prop_desc
	write(xml_out,"prop",FF[p]);
      }
      pop(xml_out);//FT_prop
    }
    else{//Default output is binary
      QDPFileWriter FFfile(file_xml, params.filename, QDPIO_SINGLEFILE, QDPIO_SERIAL, QDPIO_OPEN);
      for(int p(0) ; p<phases.numMom();p++){
	XMLBufferWriter record_xml;
	push(record_xml, "prop_desc");//write out the momemtum of each bit
	write(record_xml, "mom", phases.numToMom(p));
	pop(record_xml);
	write(FFfile,record_xml,FF[p]);
      }
    }

    // Sanity check - write out the propagator (pion) correlator in the Nd-1 direction
    {
      // Initialize the slow Fourier transform phases
      SftMom ph(0, true, Nd-1);
      
      multi1d<Double> prop_corr = sumMulti(localNorm2(quark_propagator), 
					   ph.getSet());
      
      push(xml_out, "Prop_correlator");
      write(xml_out, "prop_corr", prop_corr);
      pop(xml_out);
    }
              
    pop(xml_out);  // npr

    snoop.stop();
    QDPIO::cout << InlineNprEnv::name << ": total time = "
		<< snoop.getTimeInSeconds() 
		<< " secs" << endl;
    
    QDPIO::cout << InlineNprEnv::name << ": ran successfully" << endl;
    
    END_CODE();
  } 
  
} //name space Chroma
