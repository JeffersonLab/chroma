// $Id: t_mesons_w.cc,v 1.8 2005-03-02 00:44:19 edwards Exp $
//
//! \file
//  \brief Test the Wilson mesons() routine
//
// $Log: t_mesons_w.cc,v $
// Revision 1.8  2005-03-02 00:44:19  edwards
// Changed to new Chroma initialize/finalize format. Changed
// all XMLReader("DATA") to use a command-line param arg.
// Changed all XMLFileWriter(XMLDAT) to use the singleton instance.
//
// Revision 1.7  2005/01/14 20:13:09  edwards
// Removed all using namespace QDP/Chroma from lib files. The library
// should now be 100% in the Chroma namespace. All mainprogs need a
// using namespace Chroma.
//
// Revision 1.6  2004/02/11 12:51:35  bjoo
// Stripped out Read() and Write()
//
// Revision 1.5  2003/10/09 20:36:49  edwards
// Changed all cout/cerr to QDPIO::cout/cerr. Changed QDP_info calls
// to QDPIO::cout.
//
// Revision 1.4  2003/09/11 00:46:04  edwards
// Changed all programs to exit(0) instead of return 0
//
// Revision 1.3  2003/06/24 03:24:44  edwards
// Changed from nml to xml.
//
// Revision 1.2  2003/03/14 05:14:32  flemingg
// rewrite of mesons_w.cc to use the new SftMom class.  mesons_w.cc still
// needs to be cleaned up once the best strategy is resolved.  But for now,
// the library and test program compiles and runs.
//

#include <iostream>
#include <cstdio>
#include <time.h>

#include "chroma.h"

using namespace Chroma;

int main(int argc, char *argv[])
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);

  // Setup the layout
  const int foo[] = {4,4,4,4};
  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements
  Layout::setLattSize(nrow);
  Layout::create();

  XMLFileWriter xml("t_mesons_w.xml");
  push(xml,"t_mesons_w");

  push(xml,"lattice");
  write(xml,"Nd", Nd);
  write(xml,"Nc", Nc);
  write(xml,"Ns", Ns);
  write(xml,"nrow", nrow);
  pop(xml);

  LatticePropagator quark_prop_1, quark_prop_2;
  gaussian(quark_prop_1);
  gaussian(quark_prop_2);

  // propagation direction
  int j_decay = Nd-1 ;

  // source timeslice
  int t0 = 0 ;

  clock_t start_clock = clock() ;
  // create averaged Fourier phases with (mom)^2 <= 10
  SftMom phases(10, true, j_decay) ;

  mesons(quark_prop_1, quark_prop_2, phases, t0, xml,
         "Point_Point_Wilson_Mesons") ;

  clock_t end_clock = clock() ;

  QDPIO::cout << "Test took " << end_clock - start_clock << " clocks.\n" ;

  pop(xml);

  // Time to bolt
  Chroma::finalize();

  exit(0);
}

