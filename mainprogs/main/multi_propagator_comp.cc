// $Id: multi_propagator_comp.cc,v 3.0 2006-04-03 04:59:13 edwards Exp $
// $Log: multi_propagator_comp.cc,v $
// Revision 3.0  2006-04-03 04:59:13  edwards
// Major overhaul of fermion and gauge action interface. Basically,
// all fermacts and gaugeacts now carry out  <T,P,Q>  template parameters. These are
// the fermion type, the "P" - conjugate momenta, and "Q" - generalized coordinates
// in the sense of Hamilton's equations. The fermbc's have been rationalized to never
// be over multi1d<T>. The "createState" within the FermionAction is now fixed meaning
// the "u" fields are now from the coordinate type. There are now "ConnectState" that
// derive into FermState<T,P,Q> and GaugeState<P,Q>.
//
// Revision 2.1  2005/11/30 04:47:27  edwards
// Changed PropSource_t to PropSourceConst_t and added a new PropSourceSmear_t.
// Renamed PropSink_t to PropSinkSmear_t .
//
// Revision 2.0  2005/09/25 21:04:45  edwards
// Moved to version 2.0
//
// Revision 1.10  2005/04/11 02:01:30  edwards
// Fixed up doxygen comments.
//
// Revision 1.9  2005/03/02 00:44:18  edwards
// Changed to new Chroma initialize/finalize format. Changed
// all XMLReader("DATA") to use a command-line param arg.
// Changed all XMLFileWriter(XMLDAT) to use the singleton instance.
//
// Revision 1.8  2005/02/28 03:34:46  edwards
// Collapsed code surrounding MesPlq call to a single sub call.
//
// Revision 1.7  2005/01/14 20:13:08  edwards
// Removed all using namespace QDP/Chroma from lib files. The library
// should now be 100% in the Chroma namespace. All mainprogs need a
// using namespace Chroma.
//
// Revision 1.6  2004/12/24 04:19:22  edwards
// Removed explict FermBC args to FermAct factory functions.
//
// Revision 1.5  2004/07/28 03:08:04  edwards
// Added START/END_CODE to all routines. Changed some to not pass an
// argument.
//
// Revision 1.4  2004/04/28 16:37:53  bjoo
// Adheres to new propagator structure
//
// Revision 1.3  2004/04/26 11:19:13  bjoo
// Added support for Reading EV-s in SZIN format. Must provide e-values in XML tho.
//
// Revision 1.2  2004/04/22 16:49:23  bjoo
// Added overlap_state_info Param struct and gauge startup convenience function. Tidyed up propagator zolo4d case
//
// Revision 1.1  2004/04/20 13:08:12  bjoo
// Added multi mass component based solves and propagator collection program
//
// Revision 1.3  2004/04/16 22:03:59  bjoo
// Added Zolo 4D test files and tidyed up
//
// Revision 1.2  2004/04/16 20:18:03  bjoo
// Zolo seems to work now
//
// Revision 1.1  2004/04/16 17:04:49  bjoo
// Added multi_propagator for Zolo4D multi massery. Seems to work even
//
// Revision 1.49  2004/04/15 14:43:25  bjoo
// Added generalised qprop_io FermAct reading
//
// Revision 1.48  2004/04/06 04:20:33  edwards
// Added SZINQIO support.
//
// Revision 1.47  2004/04/01 18:10:22  edwards
// Added support for non-relativistic quark props.
//
// Revision 1.46  2004/02/23 03:13:58  edwards
// Major overhaul of input/output model! Now using EXCLUSIVELY
// SciDAC propagator format for propagators. Now, Param part of input
// files directly matches source/sink/propagator/seqprop headers
// of propagators. All ``known'' input of a propagator is derived
// from its header(s) and used for subsequent calculations.
//
// Revision 1.45  2004/02/11 12:51:35  bjoo
// Stripped out Read() and Write()
//
// Revision 1.44  2004/02/07 04:51:58  edwards
// Removed SSE hack - pushed it down into the SSE code where it belongs.
//
// Revision 1.43  2004/02/06 22:31:00  edwards
// Put in sse hack for the short term.
//
// Revision 1.42  2004/02/06 17:39:05  edwards
// Added a flush to xml_out.
//
// Revision 1.41  2004/02/05 04:18:56  edwards
// Changed call of quarkProp4 to write to xml_out instead of xml buffer.
//
// Revision 1.40  2004/02/04 17:01:55  edwards
// Changed getSubset() to the now correct getSet().
//
// Revision 1.39  2004/01/31 23:22:01  edwards
// Added proginfo call.
//
// Revision 1.38  2004/01/30 21:35:22  kostas
// added calls to calculate mres for dwf
// 
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
