// $Id: inline_propagator_w.cc,v 3.11 2007-10-13 20:46:29 edwards Exp $
/*! \file
 * \brief Inline construction of propagator
 *
 * Propagator calculations
 */

#include "fermact.h"
#include "meas/inline/hadron/inline_propagator_w.h"
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
  namespace InlinePropagatorEnv 
  { 
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlinePropagator(InlinePropagatorParams(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "PROPAGATOR";

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


  //! Propagator input
  void read(XMLReader& xml, const string& path, InlinePropagatorParams::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "gauge_id", input.gauge_id);
    read(inputtop, "source_id", input.source_id);
    read(inputtop, "prop_id", input.prop_id);
  }

  //! Propagator output
  void write(XMLWriter& xml, const string& path, const InlinePropagatorParams::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "gauge_id", input.gauge_id);
    write(xml, "source_id", input.source_id);
    write(xml, "prop_id", input.prop_id);

    pop(xml);
  }


  // Param stuff
  InlinePropagatorParams::InlinePropagatorParams() { frequency = 0; }

  InlinePropagatorParams::InlinePropagatorParams(XMLReader& xml_in, const std::string& path) 
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
  InlinePropagatorParams::writeXML(XMLWriter& xml_out, const std::string& path) 
  {
    push(xml_out, path);
    
    write(xml_out, "Param", param);
    write(xml_out, "NamedObject", named_obj);

    pop(xml_out);
  }


  // Function call
  void 
  InlinePropagator::operator()(unsigned long update_no,
			       XMLWriter& xml_out) 
  {
    // If xml file not empty, then use alternate
    if (params.xml_file != "")
    {
      string xml_file = makeXMLFileName(params.xml_file, update_no);

      push(xml_out, "propagator");
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
  InlinePropagator::func(unsigned long update_no,
			 XMLWriter& xml_out) 
  {
    START_CODE();

    QDPIO::cout << InlinePropagatorEnv::name << ": propagator calculation" << endl;

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
      QDPIO::cerr << InlinePropagatorEnv::name << ": caught dynamic cast error" 
		  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << InlinePropagatorEnv::name << ": map call failed: " << e 
		  << endl;
      QDP_abort(1);
    }
    const multi1d<LatticeColorMatrix>& u = 
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

    push(xml_out, "propagator");
    write(xml_out, "update_no", update_no);

    proginfo(xml_out);    // Print out basic program info

    // Write out the input
    params.writeXML(xml_out, "Input");

    // Write out the config header
    write(xml_out, "Config_info", gauge_xml);

    push(xml_out, "Output_version");
    write(xml_out, "out_version", 1);
    pop(xml_out);

    // Calculate some gauge invariant observables just for info.
    MesPlq(xml_out, "Observables", u);

    //
    // Read in the source along with relevant information.
    // 
    XMLReader source_file_xml, source_record_xml;

    // These pesky variables are needed in the quarkprop call - only chiral dudes
    // need this stuff, but it must be there for the bleeping virtual function
    // to live at the base class
    int t0;
    int j_decay;

    // Record the type of header
    bool make_sourceP = false;
    bool seqsourceP = false;

    QDPIO::cout << "Snarf the source from a named buffer" << endl;
    try
    {
      // Try the cast to see if this is a valid source
      LatticePropagator& source_tmp = 
	TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.source_id);

      // Snarf the source info. This is will throw if the source_id is not there
      TheNamedObjMap::Instance().get(params.named_obj.source_id).getFileXML(source_file_xml);
      TheNamedObjMap::Instance().get(params.named_obj.source_id).getRecordXML(source_record_xml);

      // Try to invert this record XML into a source struct
      // First identify what kind of source might be here
      if (source_record_xml.count("/MakeSource") != 0)
      {
	make_sourceP = true;
	MakeSourceProp_t  orig_header;
	read(source_record_xml, "/MakeSource", orig_header);

	j_decay = orig_header.source_header.j_decay;
	t0      = orig_header.source_header.t_source;
      }
      else if (source_record_xml.count("/SequentialSource") != 0)
      {
	seqsourceP = true;
	SequentialSource_t   orig_header;
	read(source_record_xml, "/SequentialSource", orig_header);

	j_decay = orig_header.seqsource_header.j_decay;
	t0      = orig_header.seqsource_header.t_sink;   // funky, but probably not really needed
      }
      else
      {
	throw std::string("No appropriate header found");
      }

      // Write out the source header
      write(xml_out, "Source_file_info", source_file_xml);
      write(xml_out, "Source_record_info", source_record_xml);
    }    
    catch (std::bad_cast)
    {
      QDPIO::cerr << InlinePropagatorEnv::name << ": caught dynamic cast error" 
		  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << InlinePropagatorEnv::name << ": error extracting source_header: " << e << endl;
      QDP_abort(1);
    }

    // Should be a valid cast now
    const LatticePropagator& quark_prop_source = 
      TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.source_id);
 
    QDPIO::cout << "Source successfully read and parsed" << endl;

    // Sanity check - write out the norm2 of the source in the Nd-1 direction
    // Use this for any possible verification
    {
      // Initialize the slow Fourier transform phases
      SftMom phases(0, true, Nd-1);

      multi1d<Double> source_corr = sumMulti(localNorm2(quark_prop_source), 
					     phases.getSet());

      push(xml_out, "Source_correlator");
      write(xml_out, "source_corr", source_corr);
      pop(xml_out);
    }

    //
    // Loop over the source color and spin, creating the source
    // and calling the relevant propagator routines. The QDP
    // terminology is that a propagator is a matrix in color
    // and spin space
    //
    try
    {
      TheNamedObjMap::Instance().create<LatticePropagator>(params.named_obj.prop_id);
    }
    catch (std::bad_cast)
    {
      QDPIO::cerr << InlinePropagatorEnv::name << ": caught dynamic cast error" 
		  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << InlinePropagatorEnv::name << ": error creating prop: " << e << endl;
      QDP_abort(1);
    }

    // Cast should be valid now
    LatticePropagator& quark_propagator = 
      TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.prop_id);
    int ncg_had = 0;

    //
    // Initialize fermion action
    //
    std::istringstream  xml_s(params.param.fermact.xml);
    XMLReader  fermacttop(xml_s);
    QDPIO::cout << "FermAct = " << params.param.fermact.id << endl;


    //
    // Try the factories
    //
    bool success = false;

    if (! success)
    {
      try
      {
	StopWatch swatch;
	swatch.reset();
	QDPIO::cout << "Try the various factories" << endl;

	// Typedefs to save typing
	typedef LatticeFermion               T;
	typedef multi1d<LatticeColorMatrix>  P;
	typedef multi1d<LatticeColorMatrix>  Q;

	// Generic Wilson-Type stuff
	Handle< FermionAction<T,P,Q> >
	  S_f(TheFermionActionFactory::Instance().createObject(params.param.fermact.id,
							       fermacttop,
							       params.param.fermact.path));

	Handle< FermState<T,P,Q> > state(S_f->createState(u));

	QDPIO::cout << "Suitable factory found: compute the quark prop" << endl;
	swatch.start();
	S_f->quarkProp(quark_propagator, 
		       xml_out, 
		       quark_prop_source,
		       t0, j_decay,
		       state, 
		       params.param.invParam, 
		       params.param.quarkSpinType,
		       params.param.obsvP,
		       ncg_had);
	swatch.stop();
	QDPIO::cout << "Propagator computed: time= " 
		    << swatch.getTimeInSeconds() 
		    << " secs" << endl;
      
	success = true;
      }
      catch (const std::string& e) 
      {
	QDPIO::cout << InlinePropagatorEnv::name << ": caught exception around quarkprop: " << e << endl;
      }
    }


    if (! success)
    {
      QDPIO::cerr << "Error: no fermact found" << endl;
      QDP_abort(1);
    }


    push(xml_out,"Relaxation_Iterations");
    write(xml_out, "ncg_had", ncg_had);
    pop(xml_out);

    // Sanity check - write out the propagator (pion) correlator in the Nd-1 direction
    {
      // Initialize the slow Fourier transform phases
      SftMom phases(0, true, Nd-1);

      multi1d<Double> prop_corr = sumMulti(localNorm2(quark_propagator), 
					   phases.getSet());

      push(xml_out, "Prop_correlator");
      write(xml_out, "prop_corr", prop_corr);
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
	MakeSourceProp_t  orig_header;
	read(source_record_xml, "/MakeSource", orig_header);

	Propagator_t  new_header;   // note, abandoning state_info
	new_header.prop_header   = params.param;
	new_header.source_header = orig_header.source_header;
	new_header.gauge_header  = orig_header.gauge_header;
	write(record_xml, "Propagator", new_header);  
      } 
      else if (seqsourceP)
      {
	SequentialSource_t  orig_header;
	read(source_record_xml, "/SequentialSource", orig_header);

	SequentialProp_t  new_header;   // note, abandoning state_info
	new_header.seqprop_header   = params.param;
	new_header.sink_header      = orig_header.sink_header;
	new_header.seqsource_header = orig_header.seqsource_header;
	new_header.forward_props    = orig_header.forward_props;
	new_header.gauge_header     = orig_header.gauge_header;
	write(record_xml, "SequentialProp", new_header);  
      }

      // Write the propagator xml info
      TheNamedObjMap::Instance().get(params.named_obj.prop_id).setFileXML(file_xml);
      TheNamedObjMap::Instance().get(params.named_obj.prop_id).setRecordXML(record_xml);

      QDPIO::cout << "Propagator successfully updated" << endl;
    }
    catch (std::bad_cast)
    {
      QDPIO::cerr << InlinePropagatorEnv::name << ": caught dynamic cast error" 
		  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << InlinePropagatorEnv::name << ": error extracting prop_header: " << e << endl;
      QDP_abort(1);
    }

    pop(xml_out);  // propagator

    snoop.stop();
    QDPIO::cout << InlinePropagatorEnv::name << ": total time = "
		<< snoop.getTimeInSeconds() 
		<< " secs" << endl;

    QDPIO::cout << InlinePropagatorEnv::name << ": ran successfully" << endl;

    END_CODE();
  } 

}
