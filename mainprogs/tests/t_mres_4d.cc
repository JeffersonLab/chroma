// $Id: t_mres_4d.cc,v 3.1 2006-09-20 20:28:06 edwards Exp $

#include "chroma.h"

using namespace Chroma;

struct Prop_t
{
  string          prop_file;
};

struct AppInput_t 
{
  multi1d<int> nrow;
  std::string  fermact;
  Cfg_t cfg;
  std::string stateInfo;
  Prop_t prop;
};

//! To insure linking of code, place the registered code flags here
/*! This is the bit of code that dictates what fermacts are in use */
bool linkage_hack()
{
  bool foo = true;

  // All actions
  foo &= WilsonTypeFermActsEnv::registerAll();

  return foo;
}


void read(XMLReader& xml, const string& path, Prop_t& input)
{
  XMLReader inputtop(xml, path);
  read(inputtop, "prop_file", input.prop_file);
}


void read(XMLReader& xml, const string& path, AppInput_t& input)
{
  XMLReader inputtop(xml, path);

  // Read the input
  try
  {
    // The parameters holds the version number
    read(inputtop, "nrow", input.nrow);

    //
    XMLReader xml_tmp(inputtop, "FermionAction");
    std::ostringstream os;
    xml_tmp.print(os);
    input.fermact = os.str();

    // Read in the gauge configuration info
    read(inputtop, "Cfg", input.cfg);

    // Read any auxiliary state information
    if( inputtop.count("StateInfo") == 1 ) {
      XMLReader xml_state_info(inputtop, "StateInfo");
      std::ostringstream os;
      xml_state_info.print(os);
      input.stateInfo = os.str();
    }
    else { 
      XMLBufferWriter s_i_xml;
      push(s_i_xml, "StateInfo");
      pop(s_i_xml);
      input.stateInfo = s_i_xml.str();
    }

    // Read in the source/propagator info
    read(inputtop, "Prop", input.prop);
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading data: " << e << endl;
    throw;
  }
}



int main(int argc, char **argv)
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);



  AppInput_t input;
  XMLReader xml_in(Chroma::getXMLInputFileName());

  try {
    read(xml_in, "/mres4D", input);
  }
  catch( const string& e) { 
    QDPIO::cerr << "Caught Exception : " << e << endl;
    QDP_abort(1);
  }


  // Setup the lattice
  Layout::setLattSize(input.nrow);
  Layout::create();

  multi1d<LatticeColorMatrix> u(Nd);
  XMLReader gauge_file_xml, gauge_xml;
  gaugeStartup(gauge_file_xml, gauge_xml, u, input.cfg);

  XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();
  push(xml_out,"t_mres4D");

  proginfo(xml_out);
  write(xml_out, "Input", xml_in);

  // Measure the plaquette on the gauge
  MesPlq(xml_out, "Observables", u);
  xml_out.flush();

  // Read the prop

  LatticePropagator quark_propagator;
  ChromaProp_t prop_header;
  PropSourceConst_t source_header;

  {
    XMLReader prop_file_xml, prop_record_xml;
    readQprop(prop_file_xml,
	      prop_record_xml,
	      quark_propagator,
	      input.prop.prop_file,
	      QDPIO_SERIAL);

    // Try to invert this record XML into a ChromaProp struct
    // Also pull out the id of this source
    try {
      read(prop_record_xml, "/Propagator/ForwardProp", prop_header);
      read(prop_record_xml, "/Propagator/PropSource", source_header);
    }
    catch (const string& e) {
      QDPIO::cerr << "Error extracting forward_prop header: " << e << endl;
      throw;
    }
  }

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
  std::istringstream  xml_s(input.fermact);
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
  std::istringstream state_info_is(input.stateInfo);
  XMLReader state_info_xml(state_info_is);
  string state_info_path="/StateInfo";


  // 
  LatticePropagator delta_prop;
  LatticePropagator deltaSq_prop;

  try 
  {
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
    QDPIO::cout << "Wilson Factory Error: " << e << endl;
    QDP_abort(1);
  }
  catch(bad_cast) { 
    QDPIO::cout << "Action entered is not suitable to be cast to OverlapFermActBase " << endl;
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
  

  pop(xml_out);
  xml_out.close();
  
  Chroma::finalize();
  exit(0);
}
