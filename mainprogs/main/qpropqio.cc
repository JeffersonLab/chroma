// $Id: qpropqio.cc,v 3.0 2006-04-03 04:59:14 edwards Exp $
/*! \file
 *  \brief Reads a scidac prop and converts only the volume format
 */

#include "chroma.h"

using namespace Chroma;


struct QpropQIO_t
{
  InlinePropagatorParams param;
  multi1d<int>     nrow;
};

// Reader for input parameters
void read(XMLReader& xml, const string& path, QpropQIO_t& input)
{
  XMLReader inputtop(xml, path);

  // Read the input
  try
  {
    // Propagator parameters
    InlinePropagatorParams foo(inputtop, "Param");
    input.param = foo;

    // Read the lattice size
    read(inputtop, "nrow", input.nrow);
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading data: " << e << endl;
    throw;
  }
}



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
  QpropQIO_t input;
  read(xml_in, "/qpropqio", input);
  InlinePropagator  meas(input.param);

  // Specify lattice size, shape, etc.
  Layout::setLattSize(input.nrow);
  Layout::create();

  // Output
  XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();

  unsigned long cur_update = 0;
  meas(cur_update, xml_out);

  xml_out.flush();
        
  END_CODE();

  // Time to bolt
  Chroma::finalize();

  exit(0);
}
