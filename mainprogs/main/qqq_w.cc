/*
 *  $Id: qqq_w.cc,v 1.9 2004-02-13 15:27:14 sbasak Exp $
 *
 *  This is the test program for the routine that reads in a quark propagator,
 *  stored in SZIN format, and computes the generalised quark propagators
 *  that will be the basis of our baryon code
 */

#include <iostream>
#include <cstdio>

#define MAIN

#include "chroma.h"

using namespace QDP;

int main(int argc, char **argv)
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  multi1d<int> nrow;

  // Read in params
  XMLReader xml_in("DATA");

  CfgType cfg_type;
  PropType prop_type;

  int j_decay;
  int version; 		// The input-parameter version

  string cfg_file;
  string prop_file;
  string qqq_file;
  string xml_in_root = "/qqq";

  QDPIO::cout << "Reading in input parameters!" << endl;

  {
	  
    XMLReader inputtop(xml_in, xml_in_root);

    try
    {

      read(inputtop, "IO_version/version", version);
      
      switch(version) 	// The parameters we read in IO version
      {
      case 1:
      {
	XMLReader paramtop(inputtop, "param");

	read(paramtop, "nrow", nrow);

	read(paramtop, "j_decay", j_decay);

	read(paramtop, "cfg_type", cfg_type);

	read(paramtop, "prop_type", prop_type);

      }
      break;

      default:
	QDP_error_exit("Unknown io version", version);

      }

      // Read in the gauge configuration file name
      read(inputtop, "Cfg/cfg_file", cfg_file);
      // Read in the propagator filename
      read(inputtop, "Prop/prop_file", prop_file);
      // Read in the generalized propagator filename
      read(inputtop, "Barcomp/qqq_file", qqq_file);
    }
    catch(const string& e)
    {
      QDP_error_exit("Error reading in qqq: %s", e.c_str());
    }
  }


  QDPIO::cout << "Input parameters read in!" << endl;

  Layout::setLattSize(nrow);
  Layout::create();


  // Read a gauge field
  multi1d<LatticeColorMatrix> u(Nd);
  XMLReader gauge_xml;

  switch (cfg_type) 
  {
  case CFG_TYPE_SZIN:
    readSzin(gauge_xml, u, cfg_file);
    break;

  case CFG_TYPE_NERSC:
    readArchiv(gauge_xml, u, cfg_file);
    break;
  default :
    QDP_error_exit("Configuration type is unsupported.");
  }


  QDPIO::cout << "Gauge field configuration read in!" << endl;

  // Initialize the slow Fourier transform phases
  SftMom phases(0, true, j_decay);

  // Now the lattice quark propagator, just a single one for this
  // example, together with the corresponding header

  LatticePropagator quark_propagator;
  // LatticePropagator quark_prop,quark_propagator;
  XMLReader prop_xml;
  PropHead          header;

  switch (prop_type)
  {
    case PROP_TYPE_SZIN:
      readSzinQprop(prop_xml, quark_propagator, prop_file);
      // readSzinQprop(prop_xml, quark_prop, prop_file);
      break;
    default :
      QDP_error_exit("Lattice propagator type is unsupported.");
  }

  QDPIO::cout << "Quark propagator read in!" << endl;

  XMLFileWriter xml_out("GP.xml");
  push(xml_out,"qqq");

  xml_out << xml_in;  // save a copy of the input

  write(xml_out, "config_info", gauge_xml);

  int t0 = 0;       // set source here to 0
  int bc_spec = 1;  // periodic
  multiNd<Complex> barprop;


  // Generalized propagator calculation

  push(xml_out, qqq_file);
  //
  // quark_propagator = quark_prop * 1.0e10;
  //
  barcomp(barprop,
	  quark_propagator,
	  quark_propagator,
	  quark_propagator,
	  phases, t0, bc_spec);


  
  // Write out the file
  writeBarcomp(qqq_file, barprop, 
	       header, header, header, j_decay);



  pop(xml_out);
  xml_out.close();
  xml_in.close();

  // Time to bolt
  QDP_finalize();

  exit(0);
}

