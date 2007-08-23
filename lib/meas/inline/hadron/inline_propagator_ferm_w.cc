// $Id: inline_propagator_ferm_w.cc,v 3.6 2007-08-23 19:02:44 edwards Exp $
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
#include "util/info/unique_id.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"
#include "meas/inline/make_xml_file.h"

#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{ 
  namespace InlinePropagatorFermEnv 
  {
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlinePropagatorFerm(InlinePropagatorFermParams(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "PROPAGATOR_FERM";

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
  }


  //! Propagator input
  void read(XMLReader& xml, const string& path, InlinePropagatorFermParams::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "gauge_id", input.gauge_id);
    read(inputtop, "source_id", input.source_id);
    read(inputtop, "prop_id", input.prop_id);
  }

  //! Propagator output
  void write(XMLWriter& xml, const string& path, const InlinePropagatorFermParams::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "gauge_id", input.gauge_id);
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
    Chroma::write(xml, "NamedObject", input.named_obj);

    pop(xml);
  }



  // Function call
  void 
  InlinePropagatorFerm::operator()(unsigned long update_no,
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
  InlinePropagatorFerm::func(unsigned long update_no,
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
      QDPIO::cerr << InlinePropagatorFermEnv::name << ": caught dynamic cast error" 
		  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << InlinePropagatorFermEnv::name << ": map call failed: " << e 
		  << endl;
      QDP_abort(1);
    }
    const multi1d<LatticeColorMatrix>& u = 
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

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
    std::istringstream  xml_s(params.param.fermact.xml);
    XMLReader  fermacttop(xml_s);
    QDPIO::cout << "FermAct = " << params.param.fermact.id << endl;


    //
    // Try the factories
    //
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

      Handle< SystemSolver<LatticeFermion> > PP = S_f->qprop(state,
							     params.param.invParam);
      
      QDPIO::cout << "Suitable factory found: compute the quark prop" << endl;
      swatch.start();

      SystemSolverResults_t res = (*PP)(quark_soln, quark_source);
      ncg_had = res.n_count;

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
      write(file_xml, "id", uniqueId());  // NOTE: new ID form
      pop(file_xml);

      XMLBufferWriter record_xml;
      if (make_sourceP)
      {
	XMLReader xml_tmp(source_record_xml, "/MakeSource");

	push(record_xml, "Propagator");
	write(record_xml, "ForwardProp", params.param);
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
