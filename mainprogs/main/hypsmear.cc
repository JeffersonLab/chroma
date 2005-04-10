/*
 *  $Id: hypsmear.cc,v 1.23 2005-04-10 16:42:10 edwards Exp $
 *
 *  This is the top-level routine for HYP smearing.
 *  It is a wrapper for Urs' and Robert's implmenetation of the HYP
 *  smearing algorithm
 */


#include <iostream>
#include <cstdio>

#include "chroma.h"

using namespace Chroma;


int main(int argc, char *argv[])
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);

  START_CODE();

  // Instantiate xml reader for DATA
  XMLReader xml_in(Chroma::getXMLInputFileName());

  // Read data
  InlineHypSmearParams input(xml_in, "/hypsmear");
  InlineHypSmear  meas(input);

  // Specify lattice size, shape, etc.
  Layout::setLattSize(input.param.nrow);
  Layout::create();

  // Read gauge field info
  Cfg_t  cfg;
  try
  {
    read(xml_in, "/hypsmear/Cfg", cfg);
  }
  catch(const string& e)
  {
    QDP_error_exit("Error reading in hypsmear: %s", e.c_str());
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
