// $Id: wallformfac.cc,v 1.3 2004-01-13 04:32:19 edwards Exp $
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
// Parameters which must be determined from the XML input
// and written to the XML output
struct Param_t
{
  FermType FermTypeP;

//  multi1d<Real> Mass;      // array of mass values

  CfgType cfg_type;        // storage order for stored gauge configuration
  int j_decay;             // direction to measure propagation

  SourceType source_type;  // source type (POINT_SOURCE, SHELL_SOURCE)

  int t_sink;

  int mom2_max;            // (mom)^2 <= mom2_max. mom2_max=7 in szin.

  SmearingParam_t  smearParam;

  multi1d<int> formfac_type;

  multi1d<int> nrow;
  multi1d<int> t_srce;
};

struct WallFormFac_input_t
{
  IO_version_t     io_version;
  Param_t          param;
  Cfg_t            cfg;
};


// Reader for input parameters
void read(XMLReader& xml, const string& path, WallFormFac_input_t& input)
{
  XMLReader inputtop(xml, path);


  // First, read the input parameter version.  Then, if this version
  // includes 'Nc' and 'Nd', verify they agree with values compiled
  // into QDP++

  // Read in the IO_version
  try
  {
    read(inputtop, "IO_version/version", input.io_version.version);
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading data: " << e << endl;
    throw;
  }


  // Currently, in the supported IO versions, there is only a small difference
  // in the inputs. So, to make code simpler, extract the common bits 

  // Read the uncommon bits first
  try
  {
    XMLReader paramtop(inputtop, "param"); // push into 'param' group

    switch (input.io_version.version) 
    {
      /**************************************************************************/
    case 1:
      /**************************************************************************/
      break;

    default:
      /**************************************************************************/
      QDPIO::cerr << "Input parameter version " << input.io_version.version 
		  << " unsupported." << endl;
      QDP_abort(1);
    }
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading data: " << e << endl;
    throw;
  }

  QDPIO::cout << " FORMFAC: Baryon form factors for Wilson fermions" << endl;

  // Read the common bits
  try 
  {
    XMLReader paramtop(inputtop, "param"); // push into 'param' group

    read(paramtop, "FermTypeP", input.param.FermTypeP);
//    read(paramtop, "Mass", input.param.Mass);

#if 0
    for (int i=0; i < input.param.Mass.size(); ++i) {
      if (toBool(input.param.Mass[i] < 0.0)) {
	QDPIO::cerr << "Unreasonable value for Mass." << endl;
	QDPIO::cerr << "  Mass[" << i << "] = " << input.param.Mass[i] << endl;
	QDP_abort(1);
      } else {
	QDPIO::cout << " Spectroscopy Mass: " << input.param.Mass[i] << endl;
      }
    }
#endif

    read(paramtop, "cfg_type", input.param.cfg_type);
    read(paramtop, "j_decay", input.param.j_decay);
    if (input.param.j_decay < 0 || input.param.j_decay >= Nd) {
      QDPIO::cerr << "Bad value: j_decay = " << input.param.j_decay << endl;
      QDP_abort(1);
    }

    read(paramtop, "source_type", input.param.source_type);

    read(paramtop, "mom2_max", input.param.mom2_max);

    if (input.param.source_type == SRC_TYPE_SHELL_SOURCE)
      read(paramtop, "SmearingParam", input.param.smearParam);

    read(paramtop, "nrow", input.param.nrow);
    read(paramtop, "t_srce", input.param.t_srce);

    // Now we read in the information associated with the sequential sources
    read(paramtop, "t_sink", input.param.t_sink);
    read(paramtop, "formfac_type", input.param.formfac_type);

#if 1
    for (int seq_src_ctr=0; seq_src_ctr<input.param.formfac_type.size(); ++seq_src_ctr) 
    {
      QDPIO::cout << "Computing formfactors for "
		  << input.param.formfac_type[seq_src_ctr] << endl;
    }
#endif
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading data: " << e << endl;
    throw;
  }


  // Read in the gauge configuration file name
  try
  {
    read(inputtop, "Cfg/cfg_file",input.cfg.cfg_file);
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading data: " << e << endl;
    throw;
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
  read(xml_in, "/wallFormFac", input);

  // Specify lattice size, shape, etc.
  Layout::setLattSize(input.param.nrow);
  Layout::create();

  // Sanity checks
  for (int i=0; i<Nd; ++i) {
    if (input.param.t_srce[i] < 0 || input.param.t_srce[i] >= input.param.nrow[i]) {
      QDPIO::cerr << "Quark propagator source coordinate incorrect." << endl;
      QDPIO::cerr << "t_srce[" << i << "] = " << input.param.t_srce[i] << endl;
      QDP_abort(1);
    }
  }

  if (input.param.t_sink < 0 || input.param.t_sink >= input.param.nrow[input.param.j_decay]) {
    QDPIO::cerr << "Sink time coordinate incorrect." << endl;
    QDPIO::cerr << "t_sink = " << input.param.t_sink << endl;
    QDP_abort(1);
  }

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

  // Write out the input
  write(xml_out, "Input", xml_in);

  // Write out the config info
  write(xml_out, "Config_info", gauge_xml);

  push(xml_out, "IO_version");
  write(xml_out, "version", input.io_version.version);
  pop(xml_out);

  push(xml_out, "Output_version");
  write(xml_out, "out_version", 1);
  pop(xml_out);

  // First calculate some gauge invariant observables just for info.
  // This is really cheap.
  Double w_plaq, s_plaq, t_plaq, link;
  MesPlq(u, w_plaq, s_plaq, t_plaq, link);

  push(xml_out, "Observables");
  Write(xml_out, w_plaq);
  Write(xml_out, s_plaq);
  Write(xml_out, t_plaq);
  Write(xml_out, link);
  pop(xml_out);

  xml_out.flush();


  // Read the quark propagator
  QDPIO::cout << "Attempt to read forward propagator" << endl;
  
  LatticePropagator forward_quark_prop;
  XMLReader forwprop_xml;
  {
    stringstream prop_file;
    prop_file << "propagator_0";
    readSzinQprop(forwprop_xml, forward_quark_prop, prop_file.str());
    
    write(xml_out, "Forward_prop_info", forwprop_xml);
  }

  QDPIO::cout << "Forward propagator successfully read" << endl;
   

  // Read the backward propagator
  QDPIO::cout << "Attempt to read backward propagator" << endl;

  LatticePropagator backward_quark_prop;
  XMLReader backprop_xml;
  {
    stringstream prop_file;
    prop_file << "backprop_0";
    readSzinQprop(backprop_xml, backward_quark_prop, prop_file.str());
    
    write(xml_out, "Backward_prop_info", backprop_xml);
  }

  QDPIO::cout << "Backward propagator successfully read" << endl;
   
  
  // Now the 3pt contractions
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
