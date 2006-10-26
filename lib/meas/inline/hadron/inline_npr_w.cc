// $Id: inline_npr_w.cc,v 1.1 2006-10-26 21:20:50 kostas Exp $
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

      // Parameters for source construction
      read(paramtop, "Param", param);

      // Read any auxiliary state information
      if( paramtop.count("Param/StateInfo") == 1 ) {
	XMLReader xml_state_info(paramtop, "Param/StateInfo");
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
  InlineNprParams::writeXML(XMLWriter& xml_out, const std::string& path) 
  {
    push(xml_out, path);
    
    write(xml_out, "Param", param);
    {
      //QDP::write(xml_out, "StateInfo", stateInfo);
      istringstream header_is(stateInfo);
      XMLReader xml_header(header_is);
      xml_out << xml_header;
    }
    write(xml_out, "NamedObject", named_obj);

    pop(xml_out);
  }


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
  InlineNpr::func(unsigned long update_no,
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

    push(xml_out, "npr");
    write(xml_out, "update_no", update_no);

    QDPIO::cout << InlineNprEnv::name << ": npr calculation" << endl;

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

    //Landau gauge fix
    // Landau Gauge fix
    coulGauge(u, nrl_gf, Nd+1, GFAccu, GFMax, OrlxDo, OrPara);
    //                   ^^^^  Makes it Landau gauge
    // gauge field is now gauge-fixed
    // Calculate some gauge invariant observables just for info.
    MesPlq(xml_out, "GaugeFixedObservables", u);

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
      QDPIO::cerr << InlineNprEnv::name << ": caught dynamic cast error" 
		  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << InlineNprEnv::name << ": error creating prop: " << e << endl;
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
//				 state_info_xml,
//				 state_info_path));  // uses phase-multiplied u-fields

	LatticePropagator quark_prop_source ;
	multi1d<int> t_source(Nd);
	
	// Set the source in the midle of the lattice
	for(int mu(0);mu<Nd;mu++)
	  t_source[mu] = Layout::lattSize()[mu]/2 ; 
	
	j_decay = Nd ; // no jdecay
	//first the point source
	make_source(quark_prop_source,state,t_source,Nd);

	QDPIO::cout << "Suitable factory found: compute the quark prop" << endl;
	swatch.start();
	S_f->quarkProp(quark_propagator, 
		       xml_out, 
		       quark_prop_source,
		       t_source[Nd-1], j_decay,
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
	QDPIO::cout << InlineNprEnv::name << ": caught exception around quarkprop: " << e << endl;
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
      push(file_xml, "npr");
      int id = 0;    // NEED TO FIX THIS - SOMETHING NON-TRIVIAL NEEDED
      write(file_xml, "id", id);
      pop(file_xml);

      XMLBufferWriter record_xml;
      if (make_sourceP)
      {
	XMLReader xml_tmp(source_record_xml, "/MakeSource");

	push(record_xml, "Npr");
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

      // Write the npr xml info
      TheNamedObjMap::Instance().get(params.named_obj.prop_id).setFileXML(file_xml);
      TheNamedObjMap::Instance().get(params.named_obj.prop_id).setRecordXML(record_xml);

      QDPIO::cout << "Propagator successfully updated" << endl;
    }
    catch (std::bad_cast)
    {
      QDPIO::cerr << InlineNprEnv::name << ": caught dynamic cast error" 
		  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << InlineNprEnv::name << ": error extracting prop_header: " << e << endl;
      QDP_abort(1);
    }

    pop(xml_out);  // npr

    snoop.stop();
    QDPIO::cout << InlineNprEnv::name << ": total time = "
		<< snoop.getTimeInSeconds() 
		<< " secs" << endl;

    QDPIO::cout << InlineNprEnv::name << ": ran successfully" << endl;

    END_CODE();
  } 

}
