/*
 *  $Id: my_qqq.cc,v 1.2 2004-02-11 12:51:35 bjoo Exp $
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

  multi1d<int> nrow(Nd);

  NmlReader nml_in("DATA");
  int version;

  // Now get info about the lattice sizes etc..

  push(nml_in,"IO_version");
  read(nml_in, "version", version);
  pop(nml_in);


  int source_type;
  int source_direction;
  int wf_type;

  switch(version){      // The parameters we read in IO version

  case 101:

    push(nml_in,"param");       // Push into param group

    read(nml_in, "source_type", source_type);  // S-wave, P-wave etc
    read(nml_in, "source_direction",source_direction);
    read(nml_in, "wf_type", wf_type);      // Point, Gaussian etc

    break;

  default:

    QDP_error_exit("Unknown io version", version);
  }

  // Now get the lattice sizes etc
  read(nml_in, "nrow", nrow);

  nml_in.close();

  Layout::setLattSize(nrow);
  Layout::create();

  // Useful parameters that should be read from an input file
  int j_decay = Nd-1;
  int length = Layout::lattSize()[j_decay]; // define the temporal direction


  // Now the lattice quark propagator, just a single one for this
  // example, together with the corresponding header

  LatticePropagator s_s_quark_propagator;
  PropHead          s_s_header;

  XMLReader xml_prop;

  //  readQprop("propagator_0", s_s_quark_propagator, s_s_header);

  readSzinQprop(xml_prop, s_s_quark_propagator, "propagator_0");

  /*
   *  If we are using the SZIN propagators, we have to "create" the
   *  header
   */

  s_s_header.kappa = 0.0;
  s_s_header.source_smearingparam = wf_type;
  s_s_header.source_type = source_type;
  s_s_header.source_direction = 0;
  s_s_header.source_laplace_power = 0;
  s_s_header.sink_smearingparam = 0;
  s_s_header.sink_direction  = 0;
  s_s_header.sink_laplace_power = 0;
    
  NmlWriter nml("NMLDAT");

  push(nml, "qqq_barcomp");


  barcomp(s_s_quark_propagator, s_s_header,
	  s_s_quark_propagator, s_s_header,
	  s_s_quark_propagator, s_s_header,
	  0, j_decay, 1, "qqq", nml);

  nml.close();

  return 0;
}
