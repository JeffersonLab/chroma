// $Id: propagator.cc,v 1.26 2004-01-02 03:19:41 edwards Exp $
/*! \file
 *  \brief Main code for propagator generation
 */

#include <iostream>
#include <cstdio>

#define MAIN

#include "chroma.h"


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
  FERM_TYPE_WILSON,
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
  Real         Kappa;      // Wilson mass
 
  CfgType  cfg_type;       // storage order for stored gauge configuration
  PropType prop_type;      // storage order for stored propagator

  enum InvType  invType;   // Inverter type
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

  read(inputtop, "source_file", input.source_file);
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
      if (ferm_type_str == "WILSON") {
	input.param.FermTypeP = FERM_TYPE_WILSON;
      } else {
	input.param.FermTypeP = FERM_TYPE_UNKNOWN;
      }
    }

    // GTF NOTE: I'm going to switch on FermTypeP here because I want
    // to leave open the option of treating masses differently.
    switch (input.param.FermTypeP) {
    case FERM_TYPE_WILSON :

      QDPIO::cout << " PROPAGATOR: Propagator for Wilson fermions" << endl;

//      read(paramtop, "numKappa", input.param.numKappa);
      read(paramtop, "Kappa", input.param.Kappa);

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
      if (cfg_type_str == "SZIN") {
	input.param.cfg_type = CFG_TYPE_SZIN;
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

  // Instantiate xml reader for DATA
  XMLReader xml_in("DATA");

  // Read data
  read(xml_in, "/propagator", input);

  // Specify lattice size, shape, etc.
  Layout::setLattSize(input.param.nrow);
  Layout::create();

  // Read in the configuration along with relevant information.
  multi1d<LatticeColorMatrix> u(Nd);
  XMLReader gauge_xml;

  switch (input.param.cfg_type) 
  {
  case CFG_TYPE_SZIN :
    readSzin(gauge_xml, u, input.cfg.cfg_file);
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
  unitarityCheck(u);

  // Calculate some gauge invariant observables just for info.
  Double w_plaq, s_plaq, t_plaq, link;
  MesPlq(u, w_plaq, s_plaq, t_plaq, link);

  push(xml_out, "Observables");
  Write(xml_out, w_plaq);
  Write(xml_out, s_plaq);
  Write(xml_out, t_plaq);
  Write(xml_out, link);
  pop(xml_out);

  xml_out.flush();

  /*
   * Construct fermionic BC. Need one for LatticeFermion and multi1d<LatticeFermion>
   * Note, the handle is on an ABSTRACT type
   */
  Handle< FermBC<LatticeFermion> >  fbc(new SimpleFermBC<LatticeFermion>(input.param.boundary));
  Handle< FermBC<multi1d<LatticeFermion> > >  fbc_a(new SimpleFermBC<multi1d<LatticeFermion> >(input.param.boundary));

  //
  // Initialize fermion action
  //
#if 1
  Real Mass = 1/(2*input.param.Kappa) - Nd;
  UnprecWilsonFermAct S_f(fbc,Mass);
#else
  Real WilsonMass = 1.5;
  Real m_q = 0.1;
  int  N5  = 8;
  UnprecDWFermActArray S_f(fbc_a,WilsonMass, m_q, N5);
//  UnprecDWFermAct S_f(fbc_a, WilsonMass, m_q);
#endif

//  FermAct = UNPRECONDITIONED_WILSON;  // global
  input.param.invType = CG_INVERTER;  // enum


  //
  // Loop over the source color and spin, creating the source
  // and calling the relevant propagator routines. The QDP
  // terminology is that a propagator is a matrix in color
  // and spin space
  //
  LatticePropagator quark_propagator;
  XMLBufferWriter xml_buf;
  int ncg_had;

  {
    Handle<const ConnectState> state(S_f.createState(u));  // uses phase-multiplied u-fields

    quarkProp4(quark_propagator, xml_buf, quark_prop_source,
  	       S_f, state, input.param.invType, input.param.RsdCG, input.param.MaxCG, ncg_had);
  }

  xml_out << xml_buf;

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
    writeSzinQprop(quark_propagator, input.prop.prop_file, input.param.Kappa);
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
