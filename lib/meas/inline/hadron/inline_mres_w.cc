// $Id: inline_mres_w.cc,v 1.3 2005-04-19 17:11:07 edwards Exp $
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

namespace Chroma 
{ 
  namespace InlineMresEnv 
  { 
    AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					    const std::string& path) 
    {
      return new InlineMres(InlineMresParams(xml_in, path));
    }

    bool registerAll()
    {
      bool foo = true;
      foo &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
      foo &= WilsonTypeFermActsEnv::registered;
      return foo;
    }

    const std::string name = "MRES";
    const bool registered = registerAll();
  };


  // Reader
  void read(XMLReader& xml, const string& path, InlineMresParams::Param_t& input)
  {
    XMLReader inputtop(xml, path);
    read(inputtop, "nrow", input.nrow);
  }


  // Writer
  void write(XMLWriter& xml, const string& path, const InlineMresParams::Param_t& input)
  {
    push(xml, path);
    write(xml, "nrow", input.nrow);
    pop(xml);
  }


  // Reader
  void read(XMLReader& xml, const string& path, InlineMresParams::Prop_t& input)
  {
    XMLReader inputtop(xml, path);
    read(inputtop, "prop_file", input.prop_file);
  }


  // Writer
  void write(XMLWriter& xml, const string& path, const InlineMresParams::Prop_t& input)
  {
    push(xml, path);
    write(xml, "prop_file", input.prop_file);
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
      read(inputtop, "Prop", prop);
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
    QDP::write(xml, "StateInfo", stateInfo);
    Chroma::write(xml, "Prop", prop);
    QDP::write(xml, "xml_file", xml_file);

    pop(xml);
  }


  // Function call
  void 
  InlineMres::operator()(const multi1d<LatticeColorMatrix>& u,
			 XMLBufferWriter& gauge_xml,
			 unsigned long update_no,
			 XMLWriter& xml_out) 
  {
    // If xml file not empty, then use alternate
    if (params.xml_file != "")
    {
      Handle<XMLFileWriter> xml(makeXMLFileWriter(params.xml_file, update_no));
      func(u, gauge_xml, update_no, *xml);
    }
    else
    {
      func(u, gauge_xml, update_no, xml_out);
    }
  }


  // Function call
  void 
  InlineMres::func(const multi1d<LatticeColorMatrix>& u,
		   XMLBufferWriter& gauge_xml,
		   unsigned long update_no,
		   XMLWriter& xml_out) 
  {
    push(xml_out, "mres");
    write(xml_out, "update_no", update_no);

    QDPIO::cout << " MRES" << endl;

    //
    // Read the prop
    //
    LatticePropagator quark_propagator;
    ChromaProp_t prop_header;
    PropSource_t source_header;
    string       stateInfo;
    XMLReader prop_file_xml, prop_record_xml;

    {
      readQprop(prop_file_xml,
		prop_record_xml,
		quark_propagator,
		params.prop.prop_file,
		QDPIO_SERIAL);

      // Try to invert this record XML into a ChromaProp struct
      // Also pull out the id of this source
      try {
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
      catch (const string& e) {
	QDPIO::cerr << "Error extracting forward_prop header: " << e << endl;
	throw;
      }
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
    multi1d<int> t_source = source_header.t_source;
    // Flags
    int t0      = t_source[j_decay];

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
    std::istringstream  xml_s(prop_header.fermact);
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

      // Generic Wilson-Type Array stuff
      FermionAction<LatticeFermion>* S_f =
	TheFermionActionFactory::Instance().createObject(fermact,
							 fermacttop,
							 fermact_path);
    
      Handle<const ConnectState> state(S_f->createState(u,
							state_info_xml,
							state_info_path)); 
    
      LinearOperator<LatticeFermion>* DelLs;

      // Possible actions
      const WilsonTypeFermAct5D< LatticeFermion, multi1d<LatticeColorMatrix> >* S_dwf = 
	dynamic_cast<const WilsonTypeFermAct5D< LatticeFermion, multi1d<LatticeColorMatrix> >*>(S_f);

      const OverlapFermActBase* S_ov = dynamic_cast<const OverlapFermActBase*>(S_f);

      if (S_dwf != 0)
      {
	DelLs = const_cast<LinearOperator<LatticeFermion>*>(S_dwf->DeltaLs(state,prop_header.invParam));
      }
      else if (S_ov != 0)
      {
	DelLs = const_cast<LinearOperator<LatticeFermion>*>(S_ov->DeltaLs(state));
      }
      else
      {
	throw string("no suitable cast found");
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
      QDPIO::cout << "Error: " << e << endl;
    }
    catch(bad_cast) { 
      QDPIO::cout << "Action entered is not suitable to be cast to DWF " << endl;
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
    int t_src = t_source[j_decay];
  
    for(int t=0; t<length; t++){
      int t_eff( (t - t_src + length) % length ) ;
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

    QDPIO::cout << "Mres ran successfully" << endl;

    END_CODE();
  } 

};
