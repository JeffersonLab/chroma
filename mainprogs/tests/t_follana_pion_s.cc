// $Id: t_follana_pion_s.cc,v 3.0 2006-04-03 04:59:14 edwards Exp $

#include <iostream>
#include <cstdio>

#include "chroma.h"
#include "io/follana_io_s.h"
#include "meas/hadron/pions_follana_s.h"

using namespace Chroma;



int main(int argc, char *argv[])
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);

  XMLReader xml_in("input.xml");

  // XML File will look like this
  // 
  // <? xml version="1.0" ?>
  // <test>
  //    <nrow>X Y Z T</nrow>
  //    <propStem>stem_for_prop_filenames</propStem>
  // </test>
  //
  multi1d<int> nrow(Nd);
  string filename_stem;
  string filename;


  read(xml_in, "/test/nrow", nrow);
  read(xml_in, "/test/propStem", filename_stem);



  Layout::setLattSize(nrow);
  Layout::create();

  QDPIO::cout << "Lattice: Lx = " << nrow[0] << " Ly = " << nrow[1] << " Lz = " << nrow[2]
	      << " Lt =" << nrow[3] << endl;


  QDPIO::cout << "About to read propagators" << endl;

  int i;
  multi1d<LatticePropagator> props(NUM_STAG_PROPS);

  for(i = 0;  i < NUM_STAG_PROPS; i++) { 
    ostringstream istring;
    istring << "." << i;

    filename = filename_stem + istring.str();;
    QDPIO::cout << "Reading Propagator from file " << filename.c_str() << endl;

    // Try and read the propagator;

    props[i] = zero;
    readQpropFollana((char *)filename.c_str(), props[i], true);
  }

  QDPIO::cout << "Computing the meaning of life..." << endl;
  
  multi2d<DComplex> pions(NUM_STAG_PIONS, nrow[3]);
  staggeredPionsFollana(props, pions, Nd-1);

  XMLFileWriter xml_out("output.xml");

  push(xml_out, "follanaIO");
  for(i=0; i < NUM_STAG_PIONS; i++) { 
    ostringstream tag;
    tag << "pion" << i;
    push(xml_out, tag.str());
    write(xml_out, "pions_i", pions[i]);
    pop(xml_out);
  }

  pop(xml_out);

  QDPIO::cout << "That's all folks" << endl;
  // Time to bolt
  Chroma::finalize();
  exit(0);
}
