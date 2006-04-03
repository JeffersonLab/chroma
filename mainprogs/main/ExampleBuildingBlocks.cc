// $Id: ExampleBuildingBlocks.cc,v 3.0 2006-04-03 04:59:13 edwards Exp $
/*! \file
 *  \brief Main code for generating building blocks
 */

//###################################################################################//
//###################################################################################//
//                                                                                   //
// ExampleBuildingBlocks.cc                                                          //
//                                                                                   //
//###################################################################################//
//###################################################################################//
//                                                                                   //
// description:                                                                      //
//                                                                                   //
// Note that the file name patterns for the u and d building blocks must be a string //
// of the form "...%c%i...%c%i...%c%i..." where the first %c%i corresponds to qz,    //
// the second to qy, and the third to qx.                                            //
//                                                                                   //
// history:                                                                          //
//                                                                                   //
// There were at least four versions of "MIT" code.  Andrew Pochinsky has a c        //
// version. Dmitri Dolgov has a c++ version.  Dru B. Renner has c and c++ versions.  //
// All were independent and checked against one another.  Of course, all were        //
// developed under the guidance of John W. Negele.  The code here is just the        //
// "Building Blocks" portion of the MIT code.                                        //
//                                                                                   //
// authors:                                                                          //
//                                                                                   //
// Dru B. Renner, dru@mit.edu, 2002 - port of Building Blocks (MIT) code to qdp++    //
//                                                                                   //
// There are others who have contributed since the code has been migrated to qdp++.  //
// The cvs log entries indicate these other authors.                                 //
//                                                                                   //
//###################################################################################//
//###################################################################################//

#include <iostream>
#include <string.h>
#include <assert.h>
#include "chroma.h"

using namespace Chroma;

//###################################################################################//
// Main Function                                                                     //
//###################################################################################//

//! Building blocks
/*! \defgroup buildingblocks Generate building blocks
 *  \ingroup main
 *
 * Main program for generating building blocks
 */

int main( int argc, char** argv )
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);

  START_CODE();

  // Instantiate xml reader for DATA
  XMLReader xml_in(Chroma::getXMLInputFileName());

  // Read data
  InlineBuildingBlocksParams input(xml_in, "/ExampleBuildingBlocks");
  InlineBuildingBlocks  meas(input);

  // Specify lattice size, shape, etc.
  Layout::setLattSize(input.param.nrow);
  Layout::create();

  // Read gauge field info
  Cfg_t  cfg;
  try
  {
    read(xml_in, "/ExampleBuildingBlocks/Cfg", cfg);
  }
  catch(const string& e)
  {
    QDP_error_exit("Error reading in ExampleBuildingBlocks: %s", e.c_str());
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

  xml_in.close();
  xml_out.close();

  END_CODE();

  // Time to bolt
  Chroma::finalize();

  exit( 0 );
}

//###################################################################################//
//###################################################################################//
