// $Id: t_propagator_s.cc,v 1.1 2003-12-10 12:38:14 bjoo Exp $
/*! \file
 *  \brief Main code for propagator generation
 */

#include <iostream>
#include <cstdio>

#define MAIN

#include "chroma.h"
#include "mesphas_follana_s.h"
#include "improvement_terms_s.h"
#include "qprop.h"
#include "asqtad_fermact_s.h"
#include "util/ferm/transf_w.h"

/*
 *  Here we have various temporary definitions
 */
enum CfgType {
  CFG_TYPE_MILC = 0,
  CFG_TYPE_NERSC,
  CFG_TYPE_SCIDAC,
  CFG_TYPE_SZIN,
  CFG_TYPE_UNKNOWN
} ;

enum PropType {
  PROP_TYPE_SCIDAC = 2,
  PROP_TYPE_SZIN,
  PROP_TYPE_UNKNOWN
} ;

enum FermType {
  FERM_TYPE_STAGGERED,
  FERM_TYPE_UNKNOWN
};



using namespace QDP;


/*
 * Input 
 */
struct IO_version_t
{
  int version;
};

// Parameters which must be determined from the XML input
// and written to the XML output
struct Param_t
{
  FermType     FermTypeP;
  Real         Mass;      // Staggered mass
  Real         u0;        // Tadpole Factor
 
  CfgType  cfg_type;       // storage order for stored gauge configuration
  PropType prop_type;      // storage order for stored propagator

//  enum InvType  invType;            // Inverter type
  Real RsdCG;
  int MaxCG;		   // Iteration parameters

  multi1d<int> nrow;
  multi1d<int> boundary;
  multi1d<int> t_srce;
};

struct Cfg_t
{
  string       cfg_file;
};

struct Prop_t
{
  string       source_file;
  string       prop_file;
};

struct Propagator_input_t
{
  IO_version_t     io_version;
  Param_t          param;
  Cfg_t            cfg;
  Prop_t           prop;
};


//
void read(XMLReader& xml, const string& path, Cfg_t& input)
{
  XMLReader inputtop(xml, path);

  read(inputtop, "cfg_file", input.cfg_file);
}


//
void read(XMLReader& xml, const string& path, Prop_t& input)
{
  XMLReader inputtop(xml, path);

//  read(inputtop, "source_file", input.source_file);
  read(inputtop, "prop_file", input.prop_file);
}



// Reader for input parameters
void read(XMLReader& xml, const string& path, Propagator_input_t& input)
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
    case 1 :
      /**************************************************************************/
      break;

    default :
      /**************************************************************************/

      QDPIO::cerr << "Input parameter version " << input.io_version.version << " unsupported." << endl;
      QDP_abort(1);
    }
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading data: " << e << endl;
    throw;
  }


  // Read the common bits
  try 
  {
    XMLReader paramtop(inputtop, "param"); // push into 'param' group

    {
      string ferm_type_str;
      read(paramtop, "FermTypeP", ferm_type_str);
      if (ferm_type_str == "STAGGERED") {
	input.param.FermTypeP = FERM_TYPE_STAGGERED;
      } else {
	input.param.FermTypeP = FERM_TYPE_UNKNOWN;
      }
    }

    // GTF NOTE: I'm going to switch on FermTypeP here because I want
    // to leave open the option of treating masses differently.
    switch (input.param.FermTypeP) {
    case FERM_TYPE_STAGGERED :

      QDPIO::cout << " PROPAGATOR: Propagator for Staggered fermions" << endl;

//      read(paramtop, "numKappa", input.param.numKappa);
      read(paramtop, "Mass", input.param.Mass);
      read(paramtop, "u0" , input.param.u0);

#if 0
      for (int i=0; i < input.param.numKappa; ++i) {
	if (toBool(input.param.Kappa[i] < 0.0)) {
	  QDPIO::cerr << "Unreasonable value for Kappa." << endl;
	  QDPIO::cerr << "  Kappa[" << i << "] = " << input.param.Kappa[i] << endl;
	  QDP_abort(1);
	} else {
	  QDPIO::cout << " Spectroscopy Kappa: " << input.param.Kappa[i] << endl;
	}
      }
#endif

      break;

    default :
      QDPIO::cerr << "Fermion type not supported." << endl;
      if (input.param.FermTypeP == FERM_TYPE_UNKNOWN) {
	QDPIO::cerr << "  FermTypeP = UNKNOWN" << endl;
      }
      QDP_abort(1);
    }

    {
      string cfg_type_str;
      read(paramtop, "cfg_type", cfg_type_str);
      if (cfg_type_str == "NERSC") {
	input.param.cfg_type = CFG_TYPE_NERSC;
      } else {
	input.param.cfg_type = CFG_TYPE_UNKNOWN;
      }
    }

    {
      string prop_type_str;
      read(paramtop, "prop_type", prop_type_str);
      if (prop_type_str == "SZIN") {
	input.param.prop_type = PROP_TYPE_SZIN;
      } else {
	input.param.prop_type = PROP_TYPE_UNKNOWN;
      }
    }

//    read(paramtop, "invType", input.param.invType);
//    input.param.invType = CG_INVERTER;   //need to fix this
    read(paramtop, "RsdCG", input.param.RsdCG);
    read(paramtop, "MaxCG", input.param.MaxCG);

    read(paramtop, "nrow", input.param.nrow);
    read(paramtop, "boundary", input.param.boundary);
    read(paramtop, "t_srce", input.param.t_srce);
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading data: " << e << endl;
    throw;
  }


  // Read in the gauge configuration file name
  try
  {
    read(inputtop, "Cfg", input.cfg);
    read(inputtop, "Prop", input.prop);
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading data: " << e << endl;
    throw;
  }
}



//! Propagator generation
/*! \defgroup propagator Propagator generation
 *  \ingroup main
 *
 * Main program for propagator generation. 
 */

int main(int argc, char **argv)
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  // Input parameter structure
  Propagator_input_t  input;

  // Phases
  multi1d<LatticeInteger> alpha(Nd); // KS Phases
  multi1d<LatticeInteger> beta(Nd);  // Auxiliary phases for this work (not needed here)

  // Instantiate xml reader for DATA
  XMLReader xml_in("DATA");

  // Read data
  read(xml_in, "/propagator", input);

  // Specify lattice size, shape, etc.
  Layout::setLattSize(input.param.nrow);
  Layout::create();

  // Read in the configuration along with relevant information.
  multi1d<LatticeColorMatrix> u(Nd);
  multi1d<LatticeColorMatrix> u_fat(Nd);
  multi1d<LatticeColorMatrix> u_triple(Nd);
  XMLReader gauge_xml;

  switch (input.param.cfg_type) 
  {
  case CFG_TYPE_NERSC :
    readArchiv(gauge_xml, u, input.cfg.cfg_file);
    break;
  default :
    QDP_error_exit("Configuration type is unsupported.");
  }


  // Read in the source along with relevant information.
  LatticePropagator quark_prop_source;
  XMLReader source_xml;

  switch (input.param.prop_type) 
  {
  case PROP_TYPE_SZIN :
//    readSzinQprop(source_xml, quark_prop_source, input.prop.source_file);
    quark_prop_source = 1;
    break;
  default :
    QDP_error_exit("Propagator type is unsupported.");
  }


  // Instantiate XML writer for XMLDAT
  XMLFileWriter xml_out("XMLDAT");
  push(xml_out, "propagator");

  // Write out the input
  write(xml_out, "Input", xml_in);

  // Write out the config header
  write(xml_out, "Config_info", gauge_xml);

  // Write out the source header
  write(xml_out, "Source_info", source_xml);

  push(xml_out, "Output_version");
  write(xml_out, "out_version", 1);
  pop(xml_out);

  xml_out.flush();


  // Check if the gauge field configuration is unitarized
//  unitarityCheck(u);   worry about this some other time


// Test to see if the phases are working!
// Turn them on and then off again and calculate plaq
//
//  mesPhasFollana(alpha, beta);
//  for(int mu=0; mu < Nd; ++mu){
//   u[mu] *= alpha[mu];
//  }
//
//  for(int mu=0; mu < Nd; ++mu){
//   u[mu] *= alpha[mu];  
//  }
//  WORKS FINE!


  // Calculate some gauge invariant observables just for info.
  Double w_plaq, s_plaq, t_plaq, link;
  MesPlq(u, w_plaq, s_plaq, t_plaq, link);

  push(xml_out, "Observables");
  Write(xml_out, w_plaq);
  Write(xml_out, s_plaq);
  Write(xml_out, t_plaq);
  Write(xml_out, link);
  pop(xml_out);

  push(xml_out, "Gauge_Field");
  Write(xml_out, u);
  pop(xml_out);

  xml_out.flush();


// Turn on KS phases
// first get them!
  mesPhasFollana(alpha, beta);
  for(int mu=0; mu<Nd; mu++) {
    u[mu] *= alpha[mu];
  }

// Make fat7 and triple links

  Fat7_Links(u, u_fat, input.param.u0);
  Triple_Links(u, u_triple, input.param.u0);

  push(xml_out, "U_FAT");
  Write(xml_out, u_fat);
  pop(xml_out);

  push(xml_out, "U_TRIPLE");
  Write(xml_out, u_triple);
  pop(xml_out);
     
  xml_out.flush();

  //
  // Initialize fermion action
  //
  AsqtadFermAct S_f(input.param.Mass);

//  FermAct = ASQTAD;  // global
//  input.param.invType = CG_INVERTER;  // enum


  //
  // Loop over the source color, creating the source
  // and calling the relevant propagator routines. The QDP
  // terminology is that a staggered propagator is a matrix in color space
  // 
  //
  LatticePropagator quark_propagator;
//  XMLBufferWriter xml_buf;
  int ncg_had = 0;

//  quarkprop4 is generic to Wilson not staggered, so explicitly
//  call qprop from here as in the past.

//  quarkProp4(quark_propagator, xml_buf, quark_prop_source,
//	     S_f, u, input.param.invType, input.param.RsdCG, input.param.MaxCG, ncg_had);

  LatticeFermion psi = zero;   // note this is ``zero'' and not 0

  for(int color_source = 0; color_source < Nc; ++color_source)
  {
    int spin_source = 0;
    LatticeFermion chi;

    // Extract a fermion source
    PropToFerm(quark_prop_source, chi, color_source, spin_source);

    push(xml_out, "SOURCE");
    Write(xml_out, psi);
    Write(xml_out, chi);
    pop(xml_out);

    // Use the last initial guess as the current guess

    // Compute the propagator for given source color/spin 
    int n_count;

    S_f.qprop(psi, u_fat, u_triple, chi, CG_INVERTER, 
              input.param.RsdCG, input.param.MaxCG, input.param.Mass, n_count);
    ncg_had += n_count;
      
    push(xml_out,"Qprop");
    write(xml_out, "Mass" , input.param.Mass);
    write(xml_out, "RsdCG", input.param.RsdCG);
    Write(xml_out, n_count);
    pop(xml_out);

    /*
     * Move the solution to the appropriate components
     * of quark propagator.
     */
    FermToProp(psi, quark_propagator, color_source, spin_source);
  }


  // Instantiate XML buffer to make the propagator header
  XMLBufferWriter prop_xml;
  push(prop_xml, "propagator");

  // Write out the input
  write(prop_xml, "Input", xml_in);

  // Write out the config header
  write(prop_xml, "Config_info", gauge_xml);

  // Write out the source header
  write(prop_xml, "Source_info", source_xml);

  pop(prop_xml);


  // Save the propagator
  switch (input.param.prop_type) 
  {
  case PROP_TYPE_SZIN:
    writeSzinQprop(quark_propagator, input.prop.prop_file, input.param.Mass);
    break;

//  case PROP_TYPE_SCIDAC:
//    writeQprop(prop_xml, quark_propagator, input.prop.prop_file);
//    break;

  default :
    QDP_error_exit("Propagator type is unsupported.");
  }


  xml_out.close();
  xml_in.close();

  // Time to bolt
  QDP_finalize();

  exit(0);
}
