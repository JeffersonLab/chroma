// $Id: wallformfac.cc,v 1.8 2004-04-04 03:54:40 edwards Exp $
/*! \file
 * \brief Main program for computing 3pt functions with a wall sink
 *
 * Main program for computing 3pt functions with a wall sink
 */

#include "chroma.h"

using namespace QDP;


/*
 * Input 
 */
struct Prop_t
{
  string       forwprop_file;
  string       backprop_file;
};


// Parameters which must be determined from the XML input
// and written to the XML output
struct Param_t
{
  int mom2_max;            // (mom)^2 <= mom2_max. mom2_max=7 in szin.

  multi1d<int> formfac_type;
  multi1d<int> nrow;
};

struct WallFormFac_input_t
{
  Param_t          param;
  Cfg_t            cfg;
  Prop_t           prop;
};


//! Propagator filenames
void read(XMLReader& xml, const string& path, Prop_t& input)
{
  XMLReader inputtop(xml, path);

  read(inputtop, "forwprop_file",input.forwprop_file);
  read(inputtop, "backprop_file",input.backprop_file);
}


//! Parameter input
void read(XMLReader& xml, const string& path, Param_t& input)
{
  XMLReader paramtop(xml, path);

  int version;
  read(paramtop, "version", version);

  switch (version) 
  {
    /**************************************************************************/
  case 1:
    break;

  default:
    /**************************************************************************/
    QDPIO::cerr << "Input parameter version " << version 
		<< " unsupported." << endl;
    QDP_abort(1);
  }

  read(paramtop, "mom2_max", param.mom2_max);
  read(paramtop, "formfac_type", param.formfac_type);
  read(paramtop, "nrow", param.nrow);
}



// Reader for input parameters
void read(XMLReader& xml, const string& path, WallFormFac_input_t& input)
{
  XMLReader inputtop(xml, path);

  // Read the input
  try
  {
    // Parameters for source construction
    read(inputtop, "Param", input.param);

    // Read in the gauge configuration info
    read(inputtop, "Cfg", input.cfg);

    // Read in the output propagator/source configuration info
    read(inputtop, "Prop", input.prop);
  }
  catch(const string& e)
  {
    QDP_error_exit("Error reading in wallformfac: %s", e.c_str());
  }
}




//! Main program for computing 3pt functions
/*! Main program */
int
main(int argc, char *argv[])
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  // Input parameter structure
  WallFormFac_input_t  input;

  // Instantiate xml reader for DATA
  XMLReader xml_in("DATA");

  // Read data
  read(xml_in, "/WallFormFac", input);

  // Specify lattice size, shape, etc.
  Layout::setLattSize(input.param.nrow);
  Layout::create();

  QDPIO::cout << " WALLFORMFAC: Form factors for Wilson-like fermions" << endl;
  QDPIO::cout << endl << "     Gauge group: SU(" << Nc << ")" << endl;
  QDPIO::cout << "     volume: " << input.param.nrow[0];
  for (int i=1; i<Nd; ++i) {
    QDPIO::cout << " x " << input.param.nrow[i];
  }
  QDPIO::cout << endl;

  // Read in the configuration along with relevant information.
  QDPIO::cout << "Attempt to initialize the gauge field" << endl;

  multi1d<LatticeColorMatrix> u(Nd);
  XMLReader gauge_xml;

  switch (input.param.cfg_type) 
  {
  case CFG_TYPE_SZIN :
    readSzin(gauge_xml, u, input.cfg.cfg_file);
    break;
  default :
    QDPIO::cerr << "Configuration type is unsupported." << endl;
    QDP_abort(1);
  }

  // Next check the gauge field configuration by reunitarizing.
  unitarityCheck(u);

  QDPIO::cout << "Gauge field successfully initialized" << endl;


  // Instantiate XML writer for XMLDAT
  XMLFileWriter xml_out("XMLDAT");
  push(xml_out, "wallFormFac");

  proginfo(xml_out);    // Print out basic program info

  // Write out the input
  write(xml_out, "Input", xml_in);

  // Write out the config info
  write(xml_out, "Config_info", gauge_xml);

  push(xml_out, "Output_version");
  write(xml_out, "out_version", 1);
  pop(xml_out);

  // First calculate some gauge invariant observables just for info.
  // This is really cheap.
  Double w_plaq, s_plaq, t_plaq, link;
  MesPlq(u, w_plaq, s_plaq, t_plaq, link);

  push(xml_out, "Observables");
  write(xml_out, "w_plaq", w_plaq);
  write(xml_out, "s_plaq", s_plaq);
  write(xml_out, "t_plaq", t_plaq);
  write(xml_out, "link", link);
  pop(xml_out);

  xml_out.flush();


  //
  // Read the quark propagator and extract headers
  //
  XMLReader forwprop_file_xml, forwprop_record_xml;
  LatticePropagator forward_quark_prop;
  ChromaProp_t forward_prop_header;
  PropSource_t forward_source_header;
  {
    QDPIO::cout << "Attempt to read forward propagator" << endl;
    readQprop(forwprop_file_xml, 
	      forwprop_record_xml, forward_quark_prop,
	      input.prop.forwprop_file, QDPIO_SERIAL);
   
    // Try to invert this record XML into a ChromaProp struct
    // Also pull out the id of this source
    try
    {
      read(forwprop_record_xml, "/Propagator/ForwardProp", forward_prop_header);
      read(forwprop_record_xml, "/Propagator/PropSource", forward_source_header);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << "Error extracting forward_prop header: " << e << endl;
      throw;
    }
  }
  QDPIO::cout << "Forward propagator successfully read" << endl;

  // Derived from input prop
  int  j_decay = forward_source_header.j_decay;
  multi1d<int> t_source = forward_source_header.t_source;

  // Sanity check - write out the norm2 of the forward prop in the j_decay direction
  // Use this for any possible verification
  {
    // Initialize the slow Fourier transform phases
    SftMom phases(0, true, j_decay);

    multi1d<Double> forward_prop_corr = sumMulti(localNorm2(forward_quark_prop), 
						 phases.getSet());

    push(xml_out, "Forward_prop_correlator");
    write(xml_out, "forward_prop_corr", forward_prop_corr);
    pop(xml_out);
  }


  // Read the backward propagator
  XMLReader backprop_file_xml, backprop_record_xml;
  LatticePropagator backward_quark_prop;
  ChromaProp_t backward_prop_header;
  PropSource_t backward_source_header;
  {
    QDPIO::cout << "Attempt to read backward propagator" << endl;
    readQprop(backprop_file_xml, 
	      backprop_record_xml, backward_quark_prop,
	      input.prop.backprop_file, QDPIO_SERIAL);
   
    // Try to invert this record XML into a ChromaProp struct
    // Also pull out the id of this source
    try
    {
      read(backprop_record_xml, "/Propagator/ForwardProp", backward_prop_header);
      read(backprop_record_xml, "/Propagator/PropSource", backward_source_header);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << "Error extracting backward_prop header: " << e << endl;
      throw;
    }
  }
  QDPIO::cout << "Backward propagator successfully read" << endl;
   
  // Sanity check - write out the norm2 of the backward prop in the j_decay direction
  // Use this for any possible verification
  {
    // Initialize the slow Fourier transform phases
    SftMom phases(0, true, j_decay);

    multi1d<Double> backward_prop_corr = sumMulti(localNorm2(backward_quark_prop), 
						 phases.getSet());

    push(xml_out, "Backward_prop_correlator");
    write(xml_out, "backward_prop_corr", backward_prop_corr);
    pop(xml_out);
  }

  
  //
  // Now the 3pt contractions
  //
  SftMom phases(input.param.mom2_max, false, input.param.j_decay);

#if 1
  wallPionFormFac(xml_out,
		  u, forward_quark_prop, backward_quark_prop, 
		  phases, 
		  input.param.t_srce[input.param.j_decay],
		  input.param.t_sink);
#else
  wallNucleonFormFac(xml_out,
		     u, forward_quark_prop, backward_quark_prop, phases, 
		     input.param.t_srce[input.param.j_decay]);
#endif

  
  // Close the namelist output file XMLDAT
  pop(xml_out);     // wallFormFac

  xml_in.close();
  xml_out.close();

  // Time to bolt
  QDP_finalize();

  exit(0);
}
