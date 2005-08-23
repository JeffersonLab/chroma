// $Id: qqq_w.cc,v 1.27 2005-08-23 19:26:51 edwards Exp $
/*! \file
 *  \brief Main code for generalized quark propagator
 *
 *  This is the main program for the routine that reads in a quark propagator,
 *  stored in SZIN format, and computes the generalised quark propagators
 *  that will be the basis of our baryon code
 */

#include "chroma.h"

using namespace Chroma;

//! QQQ measurements
/*! \defgroup qqqmain Measure baryon QQQ components
 *  \ingroup main
 *
 * Main program to measure baryon QQQ components
 */

int main(int argc, char **argv)
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);

  START_CODE();

  // Instantiate xml reader for DATA
  XMLReader xml_in(Chroma::getXMLInputFileName());

  // Read data
  InlineQQQParams input(xml_in, "/qqq_w");
  InlineQQQ  meas(input);

  // Specify lattice size, shape, etc.
  Layout::setLattSize(input.param.nrow);
  Layout::create();

  // Read gauge field info
  Cfg_t  cfg;
  try
  {
    read(xml_in, "/qqq_w/Cfg", cfg);
  }
  catch(const string& e)
  {
    QDP_error_exit("Error reading in qqq_w: %s", e.c_str());
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

