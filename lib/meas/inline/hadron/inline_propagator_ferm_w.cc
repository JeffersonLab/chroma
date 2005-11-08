// $Id: inline_propagator_ferm_w.cc,v 2.1 2005-11-08 21:16:23 edwards Exp $
/*! \file
 * \brief Inline construction of propagator returning only a single lattice fermion
 *
 * Propagator calculations
 */

#include "fermact.h"
#include "meas/inline/hadron/inline_propagator_ferm_w.h"
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
  namespace InlinePropagatorFermEnv 
  { 
    AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					    const std::string& path) 
    {
      return new InlinePropagatorFerm(InlinePropagatorFermParams(xml_in, path));
    }

    bool registerAll()
    {
      bool foo = true;
      foo &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
      foo &= WilsonTypeFermActsEnv::registered;
      return foo;
    }

    const std::string name = "PROPAGATOR_FERM";
    const bool registered = registerAll();
  };


  //! Propagator input
  void read(XMLReader& xml, const string& path, InlinePropagatorFermParams::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "source_id", input.source_id);
    read(inputtop, "prop_id", input.prop_id);
  }

  //! Propagator output
  void write(XMLWriter& xml, const string& path, const InlinePropagatorFermParams::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "source_id", input.source_id);
    write(xml, "prop_id", input.prop_id);

    pop(xml);
  }


  // Param stuff
  InlinePropagatorFermParams::InlinePropagatorFermParams() { frequency = 0; }

  InlinePropagatorFermParams::InlinePropagatorFermParams(XMLReader& xml_in, const std::string& path) 
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


  //! Propagator input
  void read(XMLReader& xml, const string& path, InlinePropagatorFermParams& input)
  {
    InlinePropagatorFermParams tmp(xml, path);
    input = tmp;
  }

  //! Propagator output
  void write(XMLWriter& xml, const string& path, const InlinePropagatorFermParams& input)
  {
    push(xml, path);
    
    Chroma::write(xml, "Param", input.param);
    {
      //QDP::write(xml, "StateInfo", stateInfo);
      istringstream header_is(input.stateInfo);
      XMLReader xml_header(header_is);
      xml << xml_header;
    }
    Chroma::write(xml, "NamedObject", input.named_obj);

    pop(xml);
  }



  // Function call
  void 
  InlinePropagatorFerm::operator()(const multi1d<LatticeColorMatrix>& u,
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
  InlinePropagatorFerm::func(const multi1d<LatticeColorMatrix>& u,
			     XMLBufferWriter& gauge_xml,
			     unsigned long update_no,
			     XMLWriter& xml_out) 
  {
    START_CODE();

    StopWatch snoop;
    snoop.reset();
    snoop.start();

    push(xml_out, "propagator");
    write(xml_out, "update_no", update_no);

    QDPIO::cout << InlinePropagatorFermEnv::name << ": propagator calculation" << endl;

    proginfo(xml_out);    // Print out basic program info

    // Write out the input
    write(xml_out, "Input", params);

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
    LatticeFermion quark_source;

    bool make_sourceP = false;
    bool seqsourceP = false;
    QDPIO::cout << "Snarf the source from a named buffer" << endl;
    try
    {
      quark_source = 
	TheNamedObjMap::Instance().getData<LatticeFermion>(params.named_obj.source_id);

      // Snarf the source info. This is will throw if the source_id is not there
      TheNamedObjMap::Instance().get(params.named_obj.source_id).getFileXML(source_file_xml);
      TheNamedObjMap::Instance().get(params.named_obj.source_id).getRecordXML(source_record_xml);

      // First identify what kind of source might be here
      if (source_record_xml.count("/MakeSource") != 0)
      {
	make_sourceP = true;
      }
      else if (source_record_xml.count("/SequentialSource") != 0)
      {
	seqsourceP = true;
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
      QDPIO::cerr << InlinePropagatorFermEnv::name << ": caught dynamic cast error" 
		  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << InlinePropagatorFermEnv::name << ": error extracting source_header: " << e << endl;
      QDP_abort(1);
    }

    QDPIO::cout << "Source successfully read and parsed" << endl;

    // Sanity check - write out the norm2 of the source in the Nd-1 direction
    // Use this for any possible verification
    {
      // Initialize the slow Fourier transform phases
      SftMom phases(0, true, Nd-1);

      multi1d<Double> source_corr = sumMulti(localNorm2(quark_source), 
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
    LatticeFermion quark_soln = zero;
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
    try
    {
      StopWatch swatch;
      swatch.reset();
      QDPIO::cout << "Try the various factories" << endl;

      // Generic Wilson-Type stuff
      Handle< FermionAction<LatticeFermion> >
	S_f(TheFermionActionFactory::Instance().createObject(fermact,
							     fermacttop,
							     fermact_path));

      
      Handle<const ConnectState> state(S_f->createState(u,
							state_info_xml,
							state_info_path));  // uses phase-multiplied u-fields

      Handle< const SystemSolver<LatticeFermion> > PP = S_f->qprop(state,
								   params.param.invParam);
      
      QDPIO::cout << "Suitable factory found: compute the quark prop" << endl;
      swatch.start();

      int ncg_had = (*PP)(quark_soln, quark_source);

      swatch.stop();
      QDPIO::cout << "Propagator computed: time= " 
		  << swatch.getTimeInSeconds() 
		  << " secs" << endl;
    }
    catch (const std::string& e) 
    {
      QDPIO::cout << InlinePropagatorFermEnv::name 
		  << ": caught exception around qprop: " << e << endl;
      QDP_abort(1);
    }


    push(xml_out,"Relaxation_Iterations");
    write(xml_out, "ncg_had", ncg_had);
    pop(xml_out);

    // Sanity check - write out the propagator (pion) correlator in the Nd-1 direction
    {
      // Initialize the slow Fourier transform phases
      SftMom phases(0, true, Nd-1);

      multi1d<Double> prop_corr = sumMulti(localNorm2(quark_soln), 
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

      // Write the propagator and info
      TheNamedObjMap::Instance().create<LatticeFermion>(params.named_obj.prop_id);
      TheNamedObjMap::Instance().getData<LatticeFermion>(params.named_obj.prop_id) = quark_soln;
      TheNamedObjMap::Instance().get(params.named_obj.prop_id).setFileXML(file_xml);
      TheNamedObjMap::Instance().get(params.named_obj.prop_id).setRecordXML(record_xml);

      QDPIO::cout << "Propagator successfully updated" << endl;
    }
    catch (std::bad_cast)
    {
      QDPIO::cerr << InlinePropagatorFermEnv::name << ": caught dynamic cast error" 
		  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << InlinePropagatorFermEnv::name << ": error extracting prop_header: " << e << endl;
      QDP_abort(1);
    }

    pop(xml_out);  // propagator

    snoop.stop();
    QDPIO::cout << InlinePropagatorFermEnv::name << ": total time = "
		<< snoop.getTimeInSeconds() 
		<< " secs" << endl;

    QDPIO::cout << InlinePropagatorFermEnv::name << ": ran successfully" << endl;

    END_CODE();
  } 

};
