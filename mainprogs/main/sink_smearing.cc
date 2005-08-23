// $Id: sink_smearing.cc,v 1.19 2005-08-23 19:26:51 edwards Exp $
/*! \file
 * \brief Main program for sink-smearing quark propagators
 *
 * This test program reads in: quark propagator(s);
 *                   does    : sink smearing on quark propagator(s);
 *                   returns : quark propagator(s).
 *
 * Quark propagator read in by this program will not be changed.
 *
 * sink quark smearing: arbitrary powers of laplacian;
 *                      gaussian smearing;
 *                      displacement.
 *       
 * gauge link: either "bare" link or ape-smeared link.
 *       
 */

#include "chroma.h"

using namespace Chroma;


//! Sink smear propagators
/*! \defgroup sinksmearmain Sink smear propagators
 *  \ingroup main
 *
 * Main program for sink-smearing quark propagators
 */

int main(int argc, char **argv)
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);

  START_CODE();

  // Instantiate xml reader for DATA
  XMLReader xml_in(Chroma::getXMLInputFileName());

  // Read data
  InlineSinkSmearParams input(xml_in, "/sink_smearing");
  InlineSinkSmear  meas(input);

  // Specify lattice size, shape, etc.
  Layout::setLattSize(input.param.nrow);
  Layout::create();

  // Read gauge field info
  Cfg_t  cfg;
  try
  {
    read(xml_in, "/sink_smearing/Cfg", cfg);
  }
  catch(const string& e)
  {
    QDP_error_exit("Error reading in sink_smearing: %s", e.c_str());
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

