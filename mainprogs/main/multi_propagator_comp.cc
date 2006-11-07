// $Id: multi_propagator_comp.cc,v 3.1 2006-11-07 21:51:25 edwards Exp $
/*! \file
 *  \brief Main code for propagator generation
 */

#include <iostream>
#include <cstdio>
#include <iomanip>
#include "chroma.h"

using namespace Chroma;


// define MRES_CALCULATION in order to run the code computing the residual mass
// and the pseudoscalar to conserved axial current correlator
#define MRES_CALCULATION

/*
 * Input 
 */
struct Prop_t
{
  string          source_file;
  string          prop_file;
  QDP_volfmt_t    prop_volfmt;
};


struct Component_t { 
  int color;
  int spin;
};

struct PropagatorComponent_input_t
{
  ChromaMultiProp_t     param;
  Cfg_t            cfg;
  Prop_t           prop;
  multi1d<Component_t> components;
};


void read(XMLReader& xml, const string& path, Component_t &comp)
{
  XMLReader top(xml,path);
  try {
    read(top, "color", comp.color);
    read(top, "spin",  comp.spin);
  }
  catch( const string& e ) {
    QDPIO::cerr << "Caught Exception : " << e << endl;
    QDP_abort(1);
  }  
  if( comp.color < 0 || comp.color >= Nc ) { 
    QDPIO::cerr << "Component color >= Nc. color = " << comp.color << endl;
    QDP_abort(1);
  }

  if( comp.spin < 0 || comp.spin >= Ns ) { 
    QDPIO::cerr << "Component spin >= Ns.  spin = " << comp.spin << endl;
    QDP_abort(1);
  }
}

void write(XMLWriter& xml, const string& path, const Component_t &comp)
{
  
  push( xml, path );

  write(xml, "color", comp.color);
  write(xml, "spin",  comp.spin);
  
  pop( xml );
}


// Propagator parameters
void read(XMLReader& xml, const string& path, Prop_t& input)
{
  XMLReader inputtop(xml, path);

  read(inputtop, "source_file", input.source_file);
  read(inputtop, "prop_file", input.prop_file);
  read(inputtop, "prop_volfmt", input.prop_volfmt);  // singlefile or multifile
}


// Reader for input parameters
void read(XMLReader& xml, const string& path, PropagatorComponent_input_t& input)
{
  XMLReader inputtop(xml, path);

  // Read the input
  try
  {
    // The parameters holds the version number
    read(inputtop, "Param", input.param);

    // Read in the gauge configuration info
    read(inputtop, "Cfg", input.cfg);

    // Read in the source/propagator info
    read(inputtop, "Prop", input.prop);

    read(inputtop, "Components", input.components);
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading data: " << e << endl;
    throw;
  }
}

// Forward declaration
void saveComponents(const ChromaMultiProp_t& param, 
		    const Prop_t& prop,
		    XMLReader& source_record_xml,
		    const Component_t& component,
		    XMLReader& gauge_xml,
		    XMLWriter& xml_out,
		    const multi1d<LatticeFermion>& psi);

//! Propagator generation
/*! \defgroup multi_propagator_comp Propagator generation
 *  \ingroup main
 *
 * Main program for propagator generation. 
 */

int main(int argc, char **argv)
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);

  START_CODE();

  // Input parameter structure
  PropagatorComponent_input_t  input;

  // Instantiate xml reader for DATA
  XMLReader xml_in(Chroma::getXMLInputFileName());

  // Read data
  read(xml_in, "/multiPropagatorComp", input);

  // Specify lattice size, shape, etc.
  Layout::setLattSize(input.param.nrow);
  Layout::create();

  QDPIO::cout << "multiPropagatorComp" << endl;

  // Read in the configuration along with relevant information.
  multi1d<LatticeColorMatrix> u(Nd);
  XMLReader gauge_file_xml, gauge_xml;

  gaugeStartup(gauge_file_xml, gauge_xml, u, input.cfg);

  // Read in the source along with relevant information.
  LatticePropagator quark_prop_source;
  XMLReader source_file_xml, source_record_xml;

  // ONLY SciDAC mode is supported for propagators!!
  readQprop(source_file_xml, 
	    source_record_xml, quark_prop_source,
	    input.prop.source_file, QDPIO_SERIAL);

  // Try to invert this record XML into a source struct
  // Also pull out the id of this source
  PropSourceConst_t source_header;

  try
  {
    read(source_record_xml, "/MakeSource/PropSource", source_header);
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error extracting source_header: " << e << endl;
    throw;
  }


  // Instantiate XML writer for XMLDAT
  XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();
  push(xml_out, "multiPropagatorComp");

  proginfo(xml_out);    // Print out basic program info

  // Write out the input
  write(xml_out, "Input", xml_in);

  // Write out the config header
  write(xml_out, "Config_info", gauge_xml);

  // Write out the source header
  write(xml_out, "Source_file_info", source_file_xml);
  write(xml_out, "Source_record_info", source_record_xml);

  push(xml_out, "Output_version");
  write(xml_out, "out_version", 1);
  pop(xml_out);

  xml_out.flush();


  // Check if the gauge field configuration is unitarized
  unitarityCheck(u);

  // Calculate some gauge invariant observables just for info.
  MesPlq(xml_out, "Observables", u);
  xml_out.flush();

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

  xml_out.flush();

  //
  // Initialize fermion action
  //
  FermionAction<LatticeFermion>* S_f_ptr = 0;
  FermionAction< multi1d<LatticeFermion> >* S_f_a_ptr = 0;

  //
  // Loop over the source color and spin, creating the source
  // and calling the relevant propagator routines. The QDP
  // terminology is that a propagator is a matrix in color
  // and spin space
  //
  int num_mass = input.param.MultiMasses.size();

  multi1d<LatticeFermion> psi(num_mass);
  int ncg_had = 0;

  //
  // Initialize fermion action
  //
  switch (input.param.FermActHandle->getFermActType()) {


    /***********************************************************************/
    /* 4D ZOlotarev                                                        */
    /***********************************************************************/
  case FERM_ACT_ZOLOTAREV_4D:
    {
      QDPIO::cout << "FERM_ACT_ZOLOTAREV_4D" << endl;
      const Zolotarev4DFermActParams& zolo4d = dynamic_cast<const Zolotarev4DFermActParams& > (*(input.param.FermActHandle));
      
      // Construct Fermact -- now uses constructor from the zolo4d params
      // struct
      S_f_ptr = new Zolotarev4DFermAct(zolo4d, xml_out);
    }
    break;
  default:
    QDPIO::cerr << "Unsupported fermion action" << endl;
    QDP_abort(1);
  }

    // Create a useable handle on the action
  // The handle now owns the pointer
  Handle< FermionAction<LatticeFermion> > S_f(S_f_ptr);
  Handle< FermionAction< multi1d<LatticeFermion> > > S_f_a(S_f_a_ptr);
      
  // FIrst we have to set up the state -- this is fermact dependent
  const ConnectState *state_ptr;

  switch(input.param.FermActHandle->getFermActType()) {
  case FERM_ACT_ZOLOTAREV_4D:
    {    
      const Zolotarev4DFermActParams& zolo4d = dynamic_cast<const Zolotarev4DFermActParams& > (*(input.param.FermActHandle));
      const Zolotarev4DFermAct& S_zolo4 = dynamic_cast<const Zolotarev4DFermAct&>(*S_f);
      
      state_ptr = S_zolo4.createState(u, zolo4d.StateInfo, xml_out,zolo4d.AuxFermActHandle->getMass());
      
    }
    break;
  default:
    QDPIO::cerr << "Unsupported fermion action (state creation)" << endl;
    QDP_abort(1);
  }
  
  // Now do the quarkprop here again depending on which of the
  // action pointers is not null
  Handle<const ConnectState> state(state_ptr);  // inserts any BC


  for(int comp = 0; comp < input.components.size(); comp++) { 
	
    LatticeFermion chi;
    PropToFerm(quark_prop_source, chi, input.components[comp].color,
	       input.components[comp].spin);
    

    // Zero out initial guess
    for(int i=0; i < num_mass; i++) {
      psi[i] = zero;
    }

    // Normalise source
    Real fact = Real(1) / sqrt(norm2(chi));
    chi *= fact;

    int n_count = 0;

    // Do the appropriate inversion.
    // Unfortunately Fermion action does not have multiQprop
    // In it so its switch and Dynamic Cast again
    switch(input.param.FermActHandle->getFermActType()) {
    case FERM_ACT_ZOLOTAREV_4D:
      {
	const Zolotarev4DFermAct& S=dynamic_cast<const Zolotarev4DFermAct&>(*S_f);
	S.multiQprop(psi, 
		     input.param.MultiMasses, state, chi, 
		     input.param.invParam.invType, 
		     input.param.invParam.RsdCG, 
		     1, 
		     input.param.invParam.MaxCG, 
		     n_count);
      }
      break;
    default:
      QDPIO::cerr << "Fermion action unsupported " << endl;
      break;
    }
    
    // Remove source normalisation
    fact = Real(1) / fact;
    for(int i=0; i < num_mass; i++) { 
      psi[i] *= fact;
      
    }
    
    push(xml_out,"Relaxation_Iterations");
    write(xml_out, "ncg_had", n_count);
    pop(xml_out);
    
    saveComponents( input.param,
		    input.prop,
		    source_record_xml,
		    input.components[comp],
		    gauge_xml,
		    xml_out,
		    psi );
    
  }

  pop(xml_out);  // propagator
  
  xml_out.close();
  xml_in.close();
  
  END_CODE();

  // Time to bolt
  QDP_finalize();
  
  exit(0);
}

void saveComponents(const ChromaMultiProp_t& param, 
		    const Prop_t& prop,
		    XMLReader& source_record_xml,
		    const Component_t& component,
		    XMLReader& gauge_xml,
		    XMLWriter& xml_out,
		    const multi1d<LatticeFermion>& psi)
{
  START_CODE();

  // Initialize the slow Fourier transform phases
  int num_mass = param.MultiMasses.size();

  SftMom phases(0, true, Nd-1);
  for(int m=0; m < num_mass; m++) { 

    multi1d<Double> prop_corr = sumMulti(localNorm2(psi[m]), 
					 phases.getSet());
	    
    push(xml_out, "PropComp_correlator");
    write(xml_out, "Number", m);
    write(xml_out, "Mass", param.MultiMasses[m]);
    write(xml_out, "prop_corr", prop_corr);
    pop(xml_out);

    XMLBufferWriter file_xml;
    push(file_xml, "propagatorComponent");
    int id = 0;    // NEED TO FIX THIS - SOMETHING NON-TRIVIAL NEEDED
    write(file_xml, "id", id);
    pop(file_xml);

    

    XMLBufferWriter record_xml;
    push(record_xml, "PropagatorComponent");
    
    // Jiggery pokery. Substitute the ChromaMultiProp_t with a 
    // ChromaProp. This is a pisser because of the FermActParams
    // THIS IS NOT TOTALLY KOSHER AS IT CHANGES THE MASS IN INPUT
    // PARAM as well. However, at this stage we have no further need
    // for input param.
    // I Will eventually write Copy Constructors.
    
    ChromaProp_t out_param(param, m);
 
    write(record_xml, "ForwardProp", out_param);
    XMLReader xml_tmp(source_record_xml, "/MakeSource");
    record_xml << xml_tmp;
    write(record_xml, "Component", component);
   
    pop(record_xml);
    
    ostringstream outfile;
    outfile << prop.prop_file << "_component_s" << component.spin
	    << "_c" << component.color << "_" 
	    << setw(3) << setfill('0') << m;
	  
    QDPIO::cout << "Attempting to write " << outfile.str() << endl;
	  
    // Write the source
    writeFermion(file_xml, record_xml, psi[m],
		 outfile.str(), prop.prop_volfmt, QDPIO_SERIAL);
  }

  Chroma::finalize();

  END_CODE();
}
