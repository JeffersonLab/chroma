// $Id: inline_propagator_w.cc,v 1.7 2005-08-24 03:12:53 edwards Exp $
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
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"
#include "meas/inline/make_xml_file.h"

namespace Chroma 
{ 
  namespace InlinePropagatorEnv 
  { 
    AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					    const std::string& path) 
    {
      return new InlinePropagator(InlinePropagatorParams(xml_in, path));
    }

    bool registerAll()
    {
      bool foo = true;
      foo &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
      foo &= WilsonTypeFermActsEnv::registered;
      return foo;
    }

    const std::string name = "PROPAGATOR";
    const bool registered = registerAll();
  };


  //! Propagator input
  void read(XMLReader& xml, const string& path, InlinePropagatorParams::Prop_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "source_file", input.source_file);
    read(inputtop, "prop_file", input.prop_file);
    read(inputtop, "prop_volfmt", input.prop_volfmt);  // singlefile or multifile
  }

  //! Propagator output
  void write(XMLWriter& xml, const string& path, const InlinePropagatorParams::Prop_t& input)
  {
    push(xml, path);

    write(xml, "source_file", input.source_file);
    write(xml, "prop_file", input.prop_file);
    write(xml, "prop_volfmt", input.prop_volfmt);

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

      // Read any auxiliary state information
      if( paramtop.count("StateInfo") == 1 ) {
	XMLReader xml_state_info(paramtop, "StateInfo");
	std::ostringstream os;
	xml_state_info.print(os);
	stateInfo = os.str();
      }
      else { 
	XMLBufferWriter s_i_xml;
	push(s_i_xml, "StateInfo");
	pop(s_i_xml);
	stateInfo = s_i_xml.printCurrentContext();
      }

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
  InlinePropagatorParams::write(XMLWriter& xml_out, const std::string& path) 
  {
    push(xml_out, path);
    
    Chroma::write(xml_out, "Param", param);
    {
      //QDP::write(xml_out, "StateInfo", stateInfo);
      istringstream header_is(stateInfo);
      XMLReader xml_header(header_is);
      xml_out << xml_header;
    }
    Chroma::write(xml_out, "Prop", prop);

    pop(xml_out);
  }


  // Function call
  void 
  InlinePropagator::operator()(const multi1d<LatticeColorMatrix>& u,
			       XMLBufferWriter& gauge_xml,
			       unsigned long update_no,
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
      func(u, gauge_xml, update_no, xml);
    }
    else
    {
      func(u, gauge_xml, update_no, xml_out);
    }
  }


  // Real work done here
  void 
  InlinePropagator::func(const multi1d<LatticeColorMatrix>& u,
			 XMLBufferWriter& gauge_xml,
			 unsigned long update_no,
			 XMLWriter& xml_out) 
  {
    START_CODE();

    push(xml_out, "propagator");
    write(xml_out, "update_no", update_no);

    QDPIO::cout << "PROPAGATOR: propagator calculation" << endl;

    //
    // Read in the source along with relevant information.
    // 
    StopWatch swatch;

    LatticePropagator quark_prop_source;
    XMLReader source_file_xml, source_record_xml;
    int t0;
    int j_decay;
    bool make_sourceP = false;
    bool seqsourceP = false;
    {
      swatch.reset();

      // ONLY SciDAC mode is supported for propagators!!
      QDPIO::cout << "Attempt to read source" << endl;
      swatch.start();
      readQprop(source_file_xml, 
		source_record_xml, quark_prop_source,
		params.prop.source_file, QDPIO_SERIAL);
      swatch.stop();
      QDPIO::cout << "Source successfully read: time= " 
		  << swatch.getTimeInSeconds() 
		  << " secs" << endl;

      // Try to invert this record XML into a source struct
      try
      {
	// First identify what kind of source might be here
	if (source_record_xml.count("/MakeSource") != 0)
	{
	  PropSource_t source_header;

	  read(source_record_xml, "/MakeSource/PropSource", source_header);
	  j_decay = source_header.j_decay;
	  t0 = source_header.t_source[j_decay];
	  make_sourceP = true;
	}
	else if (source_record_xml.count("/SequentialSource") != 0)
	{
	  ChromaProp_t prop_header;
	  PropSource_t source_header;
	  SeqSource_t seqsource_header;

	  read(source_record_xml, "/SequentialSource/SeqSource", seqsource_header);
	  // Any source header will do for j_decay
	  read(source_record_xml, "/SequentialSource/ForwardProps/elem[1]/ForwardProp", 
	       prop_header);
	  read(source_record_xml, "/SequentialSource/ForwardProps/elem[1]/PropSource", 
	       source_header);
	  j_decay = source_header.j_decay;
	  t0 = seqsource_header.t_sink;
	  seqsourceP = true;
	}
	else
	  throw std::string("No appropriate header found");
      }
      catch (const string& e) 
      {
	QDPIO::cerr << "Error extracting source_header: " << e << endl;
	QDP_abort(1);
      }
    }    
    QDPIO::cout << "Source successfully read and parsed" << endl;
 
    proginfo(xml_out);    // Print out basic program info

    // Write out the input
    params.write(xml_out, "Input");

    // Write out the config header
    write(xml_out, "Config_info", gauge_xml);

    // Write out the source header
    write(xml_out, "Source_file_info", source_file_xml);

    write(xml_out, "Source_record_info", source_record_xml);

    push(xml_out, "Output_version");
    write(xml_out, "out_version", 1);
    pop(xml_out);

    // Calculate some gauge invariant observables just for info.
    MesPlq(xml_out, "Observables", u);

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
    LatticePropagator quark_propagator;
    int ncg_had = 0;

    //
    // Initialize fermion action
    //
    std::istringstream  xml_s(params.param.fermact);
    XMLReader  fermacttop(xml_s);
    const string fermact_path = "/FermionAction";
    string fermact;

    try
    {
      read(fermacttop, fermact_path + "/FermAct", fermact);
    }
    catch (const std::string& e) 
    {
      QDPIO::cerr << "Error reading fermact: " << e << endl;
      QDP_abort(1);
    }

    QDPIO::cout << "FermAct = " << fermact << endl;


    // Deal with auxiliary (and polymorphic) state information
    // eigenvectors, eigenvalues etc. The XML for this should be
    // stored as a string called "stateInfo" in the param struct.

    // Make a reader for the stateInfo
    std::istringstream state_info_is(params.stateInfo);
    XMLReader state_info_xml(state_info_is);
    string state_info_path="/StateInfo";
    //
    // Try the factories
    //
    bool success = false;

    if (! success)
    {
      try
      {
	QDPIO::cout << "Try the various factories" << endl;

	// Generic Wilson-Type stuff
	Handle< FermionAction<LatticeFermion> >
	  S_f(TheFermionActionFactory::Instance().createObject(fermact,
							       fermacttop,
							       fermact_path));


	Handle<const ConnectState> state(S_f->createState(u,
							  state_info_xml,
							  state_info_path));  // uses phase-multiplied u-fields

	QDPIO::cout << "Suitable factory found: compute the quark prop" << endl;

	S_f->quarkProp(quark_propagator, xml_out, quark_prop_source,
		       t0, j_decay,
		       state, 
		       params.param.invParam, 
		       params.param.nonRelProp,
		       params.param.obsvP,
		       ncg_had);
      
	success = true;
      }
      catch (const std::string& e) 
      {
	QDPIO::cout << "4D: " << e << endl;
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


    // Save the propagator
    // ONLY SciDAC output format is supported!
    {
      swatch.reset();
      QDPIO::cout << "Start writing propagator" << endl;

      XMLBufferWriter file_xml;
      push(file_xml, "propagator");
      int id = 0;    // NEED TO FIX THIS - SOMETHING NON-TRIVIAL NEEDED
      write(file_xml, "id", id);
      pop(file_xml);

      XMLBufferWriter record_xml;
      if (make_sourceP)
      {
	XMLReader xml_tmp(source_record_xml, "/MakeSource");

	push(record_xml, "Propagator");
	write(record_xml, "ForwardProp", params.param);
	record_xml << params.stateInfo;  // write out the stateinfo - might be empty
	record_xml << xml_tmp;  // write out all the stuff under MakeSource
	pop(record_xml);
      } 
      else if (seqsourceP)
      {
	XMLReader xml_tmp(source_record_xml, "/SequentialSource");

	push(record_xml, "SequentialProp");
	write(record_xml, "SeqProp", params.param);
	record_xml << xml_tmp;  // write out all the stuff under SequentialSource
	pop(record_xml);
      }

      // Write the source
      swatch.start();
      writeQprop(file_xml, record_xml, quark_propagator,
		 params.prop.prop_file, params.prop.prop_volfmt, QDPIO_SERIAL);
      swatch.stop();

      QDPIO::cout << "Propagator successfully written: time= " 
		  << swatch.getTimeInSeconds() 
		  << " secs" << endl;
    }

    pop(xml_out);  // propagator

    QDPIO::cout << "Propagator ran successfully" << endl;

    END_CODE();
  } 

};
