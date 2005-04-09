// $Id: make_source.cc,v 1.40 2005-04-09 23:14:55 edwards Exp $
/*! \file
 *  \brief Main code for source generation
 */

#include <iostream>
#include <cstdio>

#include "chroma.h"

using namespace Chroma;

//! Propagator source generation
/*! \defgroup make_source Propagator source generation
 *  \ingroup main
 *
 * Main program for propagator source generation.
 */

int main(int argc, char **argv)
{
  Chroma::initialize(&argc, &argv);

  START_CODE();

  // Instantiate xml reader for DATA
  XMLReader xml_in(Chroma::getXMLInputFileName());

  // Read data
  InlineMakeSourceParams input(xml_in, "/make_source");
  InlineMakeSource  meas(input);

  Layout::setLattSize(input.param.nrow);
  Layout::create();

  QDPIO::cout << "MAKE_SOURCE: propagator source constructor" << endl;

  // Read gauge field info
  Cfg_t  cfg;
  try
  {
    read(xml_in, "/make_source/Cfg", cfg);
  }
  catch(const string& e)
  {
    QDP_error_exit("Error reading in make_source: %s", e.c_str());
  }

  // Start up the gauge field
  multi1d<LatticeColorMatrix> u(Nd);
  XMLBufferWriter config_xml;
  {
    XMLReader gauge_file_xml, gauge_xml;

    QDPIO::cout << "Initialize Gauge field" << endl;
    gaugeStartup(gauge_file_xml, gauge_xml, u, cfg);
    QDPIO::cout << "Gauge field initialized!" << endl;

    config_xml << gauge_xml;
  }

  // Check if the gauge field configuration is unitarized
  unitarityCheck(u);

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
