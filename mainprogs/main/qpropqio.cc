// $Id: qpropqio.cc,v 1.1 2005-03-22 16:53:38 edwards Exp $
/*! \file
 *  \brief Reads a scidac prop and converts only the volume format
 */

#include "chroma.h"

using namespace Chroma;

/*
 * Input 
 */

// Parameters which must be determined from the XML input
// and written to the XML output
struct Param_t
{
  multi1d<int> nrow;		// Lattice dimension
};

struct Prop_t
{
  string    prop_in_file;
  string    prop_out_file;
  QDP_volfmt_t prop_out_volfmt; // volume format (SINGLEFILE or MULTIFILE)
};

struct QpropQIO_input_t
{
  Param_t  param;
  Prop_t   prop;
};



//! Propagator parameters
void read(XMLReader& xml, const string& path, Prop_t& input)
{
  XMLReader inputtop(xml, path);

  read(inputtop, "prop_in_file", input.prop_in_file);

  read(inputtop, "prop_out_file", input.prop_out_file);
  read(inputtop, "prop_out_volfmt", input.prop_out_volfmt);  // singlefile or multifile
}


//! Parameters for running code
void read(XMLReader& xml, const string& path, Param_t& param)
{
  XMLReader paramtop(xml, path);

  int version;
  read(paramtop, "version", version);

  switch (version) 
  {
  case 1:
    /**************************************************************************/
    break;

  default :
    /**************************************************************************/
    QDPIO::cerr << "Input parameter version " << version << " unsupported." << endl;
    QDP_abort(1);
  }

  read(paramtop, "nrow", param.nrow);
}


// Reader for input parameters
void read(XMLReader& xml, const string& path, QpropQIO_input_t& input)
{
  XMLReader inputtop(xml, path);

  // Read all the input groups
  try
  {
    // Read program parameters
    read(inputtop, "Param", input.param);

    // Read in the propagator file info
    read(inputtop, "Prop", input.prop);
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading qproptransf data: " << e << endl;
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
  
  // Parameter structure for the input
  QpropQIO_input_t input;

  // Instantiate xml reader for DATA
  XMLReader xml_in(Chroma::getXMLInputFileName());

  // Read data
  read(xml_in, "/qpropqio", input);

  // Setup QDP
  Layout::setLattSize(input.param.nrow);
  Layout::create();

  QDPIO::cout << "QPROPQIO: propagator transformation utility" << endl;

  XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();
  push(xml_out, "qpropqio");

  proginfo(xml_out);    // Print out basic program info

  write(xml_out, "input", xml_in); // save a copy of the input
  xml_out.flush();
  
  /*
   * Now read them thangs...
   * NOTE: only SCIDAC format is allowed !!!
   */
  XMLReader prop_in_record_xml, prop_in_file_xml;
  LatticePropagator  prop;
  readQprop(prop_in_file_xml, prop_in_record_xml, prop, 
	    input.prop.prop_in_file, QDPIO_SERIAL);

  push(xml_out,"SciDAC_propagator");
  write(xml_out, "File_xml", prop_in_file_xml);
  write(xml_out, "Record_xml", prop_in_record_xml);
  pop(xml_out);

  // Sanity check - write out the propagator (pion) correlator in the Nd-1 direction
  {
    // Initialize the slow Fourier transform phases
    SftMom phases(0, true, Nd-1);

    multi1d<Double> prop_corr = sumMulti(localNorm2(prop), 
					 phases.getSet());

    push(xml_out, "Prop_correlator");
    write(xml_out, "prop_corr", prop_corr);
    pop(xml_out);
  }

  xml_out.flush();

  /*
   * Now write them thangs...
   */ 
  {
    XMLBufferWriter prop_out_file_xml;
    prop_out_file_xml << prop_in_file_xml;

    XMLBufferWriter prop_out_record_xml;
    prop_out_record_xml << prop_in_record_xml;
    
    // Write the source
    writeQprop(prop_out_file_xml, prop_out_record_xml, prop,
	       input.prop.prop_out_file, input.prop.prop_out_volfmt, 
	       QDPIO_SERIAL);
  }

  pop(xml_out);   // qpropqio
        
  END_CODE();

  // Time to bolt
  Chroma::finalize();

  exit(0);
}
