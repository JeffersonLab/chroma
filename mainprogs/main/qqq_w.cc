/*
 *  $Id: qqq_w.cc,v 1.8 2004-02-06 04:23:11 edwards Exp $
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

  int j_decay;
  int version; 		// The input-parameter version

  string cfg_file;

  {
    XMLReader inputtop(xml_in, "qqq");

    try
    {
      read(inputtop, "IO_version/version", version);
      
      switch(version) 	// The parameters we read in IO version
      {
      case 1:
      {
	XMLReader paramtop(inputtop, "param");

	read(paramtop, "j_decay", j_decay);

	read(paramtop, "cfg_type", cfg_type);

	read(paramtop, "nrow", nrow);
      }
      break;

      default:
	QDP_error_exit("Unknown io version", version);

      }

      // Read in the gauge configuration file name
      read(inputtop, "Cfg/cfg_file", cfg_file);
    }
    catch(const string& e)
    {
      QDP_error_exit("Error reading in qqq: %s", e.c_str());
    }
  }


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

  // Initialize the slow Fourier transform phases
  SftMom phases(0, true, j_decay);

  // Now the lattice quark propagator, just a single one for this example, together
  // with the corresponding header

  LatticePropagator s_s_quark_propagator;
  PropHead          s_s_header;
  readQprop("propagator_0", s_s_quark_propagator, s_s_header);

  XMLFileWriter xml_out("GP.xml");
  push(xml_out,"qqq");

  xml_out << xml_in;  // save a copy of the input

  write(xml_out, "config_info", gauge_xml);

  int t0 = 0;       // set source here to 0
  int bc_spec = 1;  // periodic
  multiNd<Complex> barprop;

  push(xml_out, "qqqqqq_barcomp");
  barcomp(barprop,
	  s_s_quark_propagator,
	  s_s_quark_propagator,
	  s_s_quark_propagator,
	  phases, t0, bc_spec);

  // Write out the file
  writeBarcomp("qqqqqq_barcomp", barprop, 
	       s_s_header, s_s_header, s_s_header, j_decay);
  pop(xml_out);

  /*******************************************************************/
  /***************************** P-WAVE ******************************/
  /*******************************************************************/

  {
    LatticePropagator s_p_quark_propagator;
    PropHead          s_p_header;
    readQprop("dz_propagator_0", s_p_quark_propagator, s_p_header); 

    {
      LatticePropagator p_p_quark_propagator;
      PropHead          p_p_header = s_p_header;

      push(xml_out, "qqqqDzqDzq_barcomp");
      p_p_header.sink_type=1;
      p_p_header.sink_direction=2;  
      D_j(u,  s_p_quark_propagator, p_p_quark_propagator, 2);
      barcomp(barprop,
	      s_s_quark_propagator,
	      s_s_quark_propagator,
	      p_p_quark_propagator,
	      phases, t0, bc_spec);
      writeBarcomp("qqqqDzqDzq_barcomp", barprop, 
		   s_s_header, s_s_header, p_p_header, j_decay);
      pop(xml_out);
      
      push(xml_out, "qqqqDxqDzq_barcomp");
      p_p_header.sink_type=1;
      p_p_header.sink_direction=0;  
      D_j(u,  s_p_quark_propagator, p_p_quark_propagator, 0);
      barcomp(barprop,
	      s_s_quark_propagator,
	      s_s_quark_propagator,
	      p_p_quark_propagator,
	      phases, t0, bc_spec);
      writeBarcomp("qqqqDxqDzq_barcomp", barprop, 
		   s_s_header, s_s_header, p_p_header, j_decay);
      pop(xml_out);
      
      push(xml_out, "qqqqDyqDzq_barcomp");
      p_p_header.sink_type=1;
      p_p_header.sink_direction=1;  
      D_j(u,  s_p_quark_propagator, p_p_quark_propagator, 1);
      barcomp(barprop, 
	      s_s_quark_propagator,
	      s_s_quark_propagator,
	      p_p_quark_propagator,
	      phases, t0, bc_spec);
      writeBarcomp("qqqqDyqDzq_barcomp", barprop, 
		   s_s_header, s_s_header, p_p_header, j_decay);
      pop(xml_out);
    }

    {  
      LatticePropagator p_s_quark_propagator;
      PropHead          p_s_header = s_s_header;

      push(xml_out, "qqDzqqqDzq_barcomp");
      p_s_header.sink_type=1;
      p_s_header.sink_direction=2;  
      D_j(u,  s_s_quark_propagator, p_s_quark_propagator, 2);
      barcomp(barprop,
	      s_s_quark_propagator,
	      p_s_quark_propagator,
	      s_p_quark_propagator,
	      phases, t0, bc_spec);
      writeBarcomp("qqDzqqqDzq_barcomp", barprop, 
		   s_s_header, p_s_header, s_p_header, j_decay);
      pop(xml_out);
      
      push(xml_out, "qqDxqqqDzq_barcomp");
      p_s_header.sink_type=1;
      p_s_header.sink_direction=0;  
      D_j(u,  s_s_quark_propagator, p_s_quark_propagator, 0);
      barcomp(barprop,
	      s_s_quark_propagator,
	      p_s_quark_propagator,
	      s_p_quark_propagator,
	      phases, t0, bc_spec);
      writeBarcomp("qqDxqqqDzq_barcomp", barprop, 
		   s_s_header, p_s_header, s_p_header, j_decay);
      pop(xml_out);
      
      push(xml_out, "qqDyqqqDzq_barcomp");
      p_s_header.sink_type=1;
      p_s_header.sink_direction=1;  
      D_j(u,  s_s_quark_propagator, p_s_quark_propagator, 1);
      barcomp(barprop,
	      s_s_quark_propagator,
	      p_s_quark_propagator,
	      s_p_quark_propagator,
	      phases, t0, bc_spec);
      writeBarcomp("qqDyqqqDzq_barcomp", barprop, 
		   s_s_header, p_s_header, s_p_header, j_decay);
      pop(xml_out);
    }
  }

  /*******************************************************************/
  /***************************** D-WAVE ******************************/
  /*******************************************************************/

  {
    LatticePropagator s_dydz_quark_propagator;
    PropHead          s_dydz_header;
    readQprop("dydz_propagator_0", s_dydz_quark_propagator, s_dydz_header);
    // dydz_propagator_0 contains s_dydz_quark_propagators

    {
      LatticePropagator dydz_dydz_quark_propagator;
      PropHead          dydz_dydz_header = s_dydz_header;

      push(xml_out, "qqqqDyDzqDyDzq_barcomp");
      dydz_dydz_header.sink_type=2;
      dydz_dydz_header.sink_direction=12; // mean dydz
      DjDk(u, s_dydz_quark_propagator, dydz_dydz_quark_propagator, 12);
      //                                                    1->dy 2->dz
      barcomp(barprop,
	      s_s_quark_propagator,
	      s_s_quark_propagator,
	      dydz_dydz_quark_propagator,
	      phases, t0, bc_spec);
      writeBarcomp("qqqqDyDzqDyDzq_barcomp", barprop, 
		   s_s_header, s_s_header, dydz_dydz_header, j_decay);
      pop(xml_out);
    }

    {
      LatticePropagator dydz_s_quark_propagator;
      PropHead          dydz_s_header = s_s_header;

      push(xml_out, "qqDyDzqqqDyDzq_barcomp");
      dydz_s_header.sink_type=2;
      dydz_s_header.sink_direction=12;  // mean dydz
      DjDk(u, s_s_quark_propagator, dydz_s_quark_propagator, 12);
      barcomp(barprop,
	      s_s_quark_propagator,
	      dydz_s_quark_propagator, 
	      s_dydz_quark_propagator, 
	      phases, t0, bc_spec);
      writeBarcomp("qqDyDzqqqDyDzq_barcomp", barprop, 
		   s_s_header, dydz_s_header, s_dydz_header, j_decay);
      pop(xml_out);
    }
  }


  {
    LatticePropagator s_dzdz_quark_propagator;
    PropHead          s_dzdz_header;
    readQprop("dzdz_propagator_0", s_dzdz_quark_propagator, s_dzdz_header);
    // dzdz_propagator_0 contains s_dzdz_quark_propagators

    {
      LatticePropagator dzdz_dzdz_quark_propagator;
      PropHead          dzdz_dzdz_header = s_dzdz_header;

      push(xml_out, "qqqqDzDzqDzDzq_barcomp");
      dzdz_dzdz_header.sink_type=2;
      dzdz_dzdz_header.sink_direction=22; // mean dzdz
      DjDk(u, s_dzdz_quark_propagator, dzdz_dzdz_quark_propagator, 22);
      barcomp(barprop,
	      s_s_quark_propagator,
	      s_s_quark_propagator,
	      dzdz_dzdz_quark_propagator, 
	      phases, t0, bc_spec);
      writeBarcomp("qqqqDzDzqDzDzq_barcomp", barprop, 
		   s_s_header, s_s_header, dzdz_dzdz_header, j_decay);
      pop(xml_out);
    }

    {
      LatticePropagator dzdz_s_quark_propagator;
      PropHead          dzdz_s_header = s_s_header;

      push(xml_out, "qqDzDzqqqDzDzq_barcomp");
      dzdz_s_header.sink_type=2;
      dzdz_s_header.sink_direction=22;
      DjDk(u, s_s_quark_propagator, dzdz_s_quark_propagator, 22);
      barcomp(barprop,
	      s_s_quark_propagator,
	      dzdz_s_quark_propagator,
	      s_dzdz_quark_propagator, 
	      phases, t0, bc_spec);
      writeBarcomp("qqDzDzqqqDzDzq_barcomp", barprop, 
		   s_s_header, dzdz_s_header, s_dzdz_header, j_decay);
      pop(xml_out);
    }

    {
      LatticePropagator dydy_dzdz_quark_propagator;
      PropHead          dydy_dzdz_header = s_dzdz_header;

      push(xml_out, "qqqqDyDyqDzDzq_barcomp");
      dydy_dzdz_header.sink_type=2;
      dydy_dzdz_header.sink_direction=11;  // mean dydy
      DjDk(u, s_dzdz_quark_propagator, dydy_dzdz_quark_propagator, 11);
      barcomp(barprop,
	      s_s_quark_propagator,
	      s_s_quark_propagator,
	      dydy_dzdz_quark_propagator, 
	      phases, t0, bc_spec);
      writeBarcomp("qqqqDyDyqDzDzq_barcomp", barprop, 
		   s_s_header, s_s_header, dydy_dzdz_header, j_decay);
      pop(xml_out);
    }
      
    {
      LatticePropagator dydy_s_quark_propagator;
      PropHead          dydy_s_header = s_s_header;

      push(xml_out, "qqDyDyqqqDzDzq_barcomp");
      dydy_s_header.sink_type=2;
      dydy_s_header.sink_direction=11;  // mean dydy
      DjDk(u, s_s_quark_propagator, dydy_s_quark_propagator, 11);
      barcomp(barprop,
	      s_s_quark_propagator,
	      dydy_s_quark_propagator,
	      s_dzdz_quark_propagator,
	      phases, t0, bc_spec);
      writeBarcomp("qqDyDyqqqDzDzq_barcomp", barprop, 
		   s_s_header, dydy_s_header, s_dzdz_header, j_decay);
      pop(xml_out);
    }
  }


  pop(xml_out);  // make_source
  xml_out.close();
  xml_in.close();

  // Time to bolt
  QDP_finalize();

  exit(0);
}
