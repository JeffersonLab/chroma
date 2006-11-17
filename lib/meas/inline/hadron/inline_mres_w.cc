// $Id: inline_mres_w.cc,v 3.4 2006-11-17 02:17:31 edwards Exp $
/*! \file
 * \brief Inline construction of mres
 *
 * Mres calculations
 */

#include "fermact.h"
#include "meas/inline/hadron/inline_mres_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "util/ft/sftmom.h"
#include "util/ferm/transf.h"
#include "util/info/proginfo.h"
#include "io/qprop_io.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"
#include "actions/ferm/fermacts/overlap_fermact_base_w.h"
#include "meas/inline/make_xml_file.h"
#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{ 
  namespace InlineMresEnv 
  { 
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineMres(InlineMresParams(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "MRES";

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


  // Reader
  void read(XMLReader& xml, const string& path, InlineMresParams::Param_t& input)
  {
    XMLReader inputtop(xml, path);

    if (inputtop.count("FermionAction") != 0)
      input.fermact = readXMLGroup(inputtop, "FermionAction", "FermAct");
  }


  // Writer
  void write(XMLWriter& xml, const string& path, const InlineMresParams::Param_t& input)
  {
    push(xml, path);

    QDPIO::cout << "write mresparams: fermact=XX" << input.fermact.xml << "XX\n";

    if (input.fermact.xml != "")
      xml << input.fermact.xml;
      
    pop(xml);
  }


  // Reader
  void read(XMLReader& xml, const string& path, InlineMresParams::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "gauge_id", input.gauge_id);
    read(inputtop, "prop_id", input.prop_id);
  }


  // Writer
  void write(XMLWriter& xml, const string& path, const InlineMresParams::NamedObject_t& input)
  {
    push(xml, path);
    write(xml, "gauge_id", input.gauge_id);
    write(xml, "prop_id", input.prop_id);
    pop(xml);
  }


  // Param stuff
  InlineMresParams::InlineMresParams() { frequency = 0; }

  InlineMresParams::InlineMresParams(XMLReader& xml, const std::string& path) 
  {
    XMLReader inputtop(xml, path);

    // Read the input
    try
    {
      // The parameters holds the version number
      read(inputtop, "Param", param);

      // Read any auxiliary state information
      if( inputtop.count("StateInfo") == 1 ) {
	XMLReader xml_state_info(inputtop, "StateInfo");
	std::ostringstream os;
	xml_state_info.print(os);
	stateInfo = os.str();
      }
      else { 
	XMLBufferWriter s_i_xml;
	push(s_i_xml, "StateInfo");
	pop(s_i_xml);
	stateInfo = s_i_xml.str();
      }

      // Read in the source/propagator info
      read(inputtop, "NamedObject", named_obj);

      // Possible alternate XML file pattern
      if (inputtop.count("xml_file") != 0) 
      {
	read(inputtop, "xml_file", xml_file);
      }
    }
    catch (const string& e) 
    {
      QDPIO::cerr << "Error reading data: " << e << endl;
      throw;
    }
  }

  void
  InlineMresParams::write(XMLWriter& xml, const std::string& path) 
  {
    push(xml, path);
    
    Chroma::write(xml, "Param", param);
    {
      //QDP::write(xml, "StateInfo", stateInfo);
      istringstream header_is(stateInfo);
      XMLReader xml_header(header_is);
      xml << xml_header;
    }
    Chroma::write(xml, "NamedObject", named_obj);
    QDP::write(xml, "xml_file", xml_file);

    pop(xml);
  }


  // Function call
  void 
  InlineMres::operator()(unsigned long update_no,
			 XMLWriter& xml_out) 
  {
    // If xml file not empty, then use alternate
    if (params.xml_file != "")
    {
      string xml_file = makeXMLFileName(params.xml_file, update_no);

      push(xml_out, "mres");
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
  InlineMres::func(unsigned long update_no,
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
      QDPIO::cerr << InlineMresEnv::name << ": caught dynamic cast error" 
		  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << InlineMresEnv::name << ": map call failed: " << e 
		  << endl;
      QDP_abort(1);
    }
    const multi1d<LatticeColorMatrix>& u = 
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

    push(xml_out, "mres");
    write(xml_out, "update_no", update_no);

    QDPIO::cout << " MRES" << endl;

    //
    // Read the prop
    //
    LatticePropagator quark_propagator;
    ChromaProp_t prop_header;
    PropSourceConst_t source_header;
    string       stateInfo;
    XMLReader prop_file_xml, prop_record_xml;

    try
    {
      // Snarf the data into a copy
      quark_propagator =
	TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.prop_id);

      // Snarf the source info. This is will throw if the source_id is not there
      XMLReader prop_file_xml, prop_record_xml;
      TheNamedObjMap::Instance().get(params.named_obj.prop_id).getFileXML(prop_file_xml);
      TheNamedObjMap::Instance().get(params.named_obj.prop_id).getRecordXML(prop_record_xml);
   
      // Try to invert this record XML into a ChromaProp struct
      // Also pull out the id of this source
      {
	read(prop_record_xml, "/Propagator/ForwardProp", prop_header);
	read(prop_record_xml, "/Propagator/PropSource", source_header);

	// Read any auxiliary state information
	if( prop_record_xml.count("/Propagator/StateInfo") == 1 ) 
	{
	  XMLReader xml_state_info(prop_record_xml, "/Propagator/StateInfo");
	  std::ostringstream os;
	  xml_state_info.print(os);
	  stateInfo = os.str();
	}
	else 
	{ 
	  // The user better have supplied something
	  stateInfo = params.stateInfo;
	}
      }
    }
    catch( std::bad_cast ) 
    {
      QDPIO::cerr << InlineMresEnv::name << ": caught dynamic cast error" 
		  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << InlineMresEnv::name << ": map call failed: " << e 
		  << endl;
      QDP_abort(1);
    }
      
    proginfo(xml_out);    // Print out basic program info

    // Write out the input
    params.write(xml_out, "Input");

    // Write out the config header
    write(xml_out, "Config_info", gauge_xml);

    // Write out the propagator header
    write(xml_out, "Propagator_file_info", prop_file_xml);
    write(xml_out, "Propagator_record_info", prop_record_xml);

    push(xml_out, "Output_version");
    write(xml_out, "out_version", 1);
    pop(xml_out);

    // Calculate some gauge invariant observables just for info.
    MesPlq(xml_out, "Observables", u);

    // Basic parameters
    int j_decay = source_header.j_decay;
    int t0      = source_header.t_source;

    // Initialize the slow Fourier transform phases
    SftMom phases(0, true, j_decay);
    {
      multi1d<Double> forward_prop_corr = sumMulti(localNorm2(quark_propagator),
						   phases.getSet());
                                                                                
      push(xml_out, "Forward_prop_correlator");
      write(xml_out, "forward_prop_corr", forward_prop_corr);
      pop(xml_out);
    }


    //
    // Initialize fermion action
    //
    string ferm_act_str;
    if (params.param.fermact.xml == "")
      ferm_act_str = prop_header.fermact.xml;
    else
      ferm_act_str = params.param.fermact.xml;
      
    std::istringstream  xml_s(ferm_act_str);
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
      throw;
    }

    QDPIO::cout << "FermAct = " << fermact << endl;
 
    // Make a reader for the stateInfo
    std::istringstream state_info_is(stateInfo);
    XMLReader state_info_xml(state_info_is);
    string state_info_path="/StateInfo";


    // 
    LatticePropagator delta_prop;
    LatticePropagator deltaSq_prop;


    try 
    { 
      QDPIO::cout << "Try generic WilsonTypeFermAct actions" << endl;

      // Typedefs to save typing
      typedef LatticeFermion               T;
      typedef multi1d<LatticeColorMatrix>  P;
      typedef multi1d<LatticeColorMatrix>  Q;

      // Generic Wilson-Type Array stuff
      FermionAction<T,P,Q>* S_f =
	TheFermionActionFactory::Instance().createObject(fermact,
							 fermacttop,
							 fermact_path);
    
      Handle< FermState<T,P,Q> > state(S_f->createState(u,
							state_info_xml,
							state_info_path)); 
      
      LinearOperator<T>* DelLs;

      // Possible actions
      WilsonTypeFermAct5D<T,P,Q>* S_dwf = dynamic_cast< WilsonTypeFermAct5D<T,P,Q>*>(S_f);
      OverlapFermActBase* S_ov = dynamic_cast< OverlapFermActBase*>(S_f);

      if (S_dwf != 0)
      {
	DelLs = S_dwf->DeltaLs(state,prop_header.invParam);
      }
      else if (S_ov != 0)
      {
	DelLs = S_ov->DeltaLs(state);
      }
      else
      {
	throw string("No suitable cast found");
      }
	
      for(int col = 0; col < Nc; col++) 
      {
	for(int spin = 0; spin < Ns; spin++) 
	{
	  LatticeFermion tmp1, tmp2;

	  // Move component from prop to ferm
	  PropToFerm(quark_propagator, tmp1, col, spin);

	  // Apply DeltaLs -> tmp2 = Delta_Ls tmp1
	  (*DelLs)(tmp2, tmp1, PLUS);
	
	  FermToProp(tmp2, delta_prop, col, spin);
	}
      }
    
      delete DelLs;
      delete S_f;
    }
    catch(const string& e) { 
      QDPIO::cerr << InlineMresEnv::name << ": Error:" << e << endl;
      QDP_abort(1);
    }
    catch(bad_cast) { 
      QDPIO::cout << InlineMresEnv::name << ": Action entered is not suitable to be cast to DWF " << endl;
      QDP_abort(1);
    }


    multi1d<Double> pseudo_prop_corr = sumMulti(localNorm2(quark_propagator),
						phases.getSet());

    multi1d<DComplex> delta_prop_corr = sumMulti(trace(adj(quark_propagator)*delta_prop),
						 phases.getSet());
  
    multi1d<DComplex> deltaSq_prop_corr = sumMulti(trace(adj(delta_prop)*delta_prop), 
						   phases.getSet());
   
    int length = pseudo_prop_corr.size();
    multi1d<Real> shifted_pseudo(length);
    multi1d<Real> shifted_delta(length);
    multi1d<Real> shifted_deltaSq(length);
  
    for(int t=0; t<length; t++)
    {
      int t_eff( (t - t0 + length) % length ) ;
      shifted_pseudo[t_eff] = real(pseudo_prop_corr[t]);
      shifted_delta[t_eff]  = real(delta_prop_corr[t]);
      shifted_deltaSq[t_eff]= real(deltaSq_prop_corr[t]);
    }
  
    push(xml_out, "DeltaProp_correlator");
    write(xml_out, "delta_prop_corr", shifted_delta);
    pop(xml_out);
  
    push(xml_out, "DeltaSqProp_correlator");
    write(xml_out, "delta_sq_prop_corr", shifted_deltaSq);
    pop(xml_out);
  
    push(xml_out, "PsuedoPseudo_correlator");
    write(xml_out, "pseudo_prop_corr", shifted_pseudo);
    pop(xml_out);
  
    pop(xml_out); // mres

    snoop.stop();
    QDPIO::cout << InlineMresEnv::name << ": total time = "
		<< snoop.getTimeInSeconds() 
		<< " secs" << endl;

    QDPIO::cout << InlineMresEnv::name << ": ran successfully" << endl;

    END_CODE();
  } 

};
