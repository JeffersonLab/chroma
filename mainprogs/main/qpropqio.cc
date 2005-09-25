// $Id: qpropqio.cc,v 2.0 2005-09-25 21:04:46 edwards Exp $
/*! \file
 *  \brief Reads a scidac prop and converts only the volume format
 */

#include "chroma.h"

using namespace Chroma;


//! SciDAC prop conversion program
/*! \defgroup qpropqio SciDAC prop conversion program
 *  \ingroup main
 *
 * Main program for transforming propagator formats
 */

int main(int argc, char *argv[])
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);

  START_CODE();

  // Instantiate xml reader for DATA
  XMLReader xml_in(Chroma::getXMLInputFileName());

  // Read data
  InlinePropagatorParams input(xml_in, "/qpropqio");
  InlinePropagator  meas(input);

  // Specify lattice size, shape, etc.
  Layout::setLattSize(input.param.nrow);
  Layout::create();

  // Gauge field not needed, fake one
  multi1d<LatticeColorMatrix> u(Nd);
  u = zero;
  XMLBufferWriter config_xml;
  push(config_xml, "config");
  pop(config_xml);

  // Output
  XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();

  unsigned long cur_update = 0;
  meas(u, config_xml, cur_update, xml_out);

  xml_out.flush();
        
  END_CODE();

  // Time to bolt
  Chroma::finalize();

  exit(0);
}
