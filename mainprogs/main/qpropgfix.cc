// $Id: qpropgfix.cc,v 3.0 2006-04-03 04:59:13 edwards Exp $
/*! \file
 *  \brief Applies gauge transformation matrices on a propagator
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
  string    gfix_in_file;

  string    prop_out_file;
  QDP_volfmt_t prop_out_volfmt; // volume format (SINGLEFILE or MULTIFILE)
};

struct QpropGFix_input_t
{
  Param_t  param;
  Cfg_t    cfg;
  Prop_t   prop;
};



//! Propagator parameters
void read(XMLReader& xml, const string& path, Prop_t& input)
{
  XMLReader inputtop(xml, path);

  read(inputtop, "prop_in_file", input.prop_in_file);
  read(inputtop, "gfix_in_file", input.gfix_in_file);

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
void read(XMLReader& xml, const string& path, QpropGFix_input_t& input)
{
  XMLReader inputtop(xml, path);

  // Read all the input groups
  try
  {
    // Read program parameters
    read(inputtop, "Param", input.param);

    // Read in the gauge configuration info
    read(inputtop, "Cfg", input.cfg);

    // Read in the propagator file info
    read(inputtop, "Prop", input.prop);
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading qproptransf data: " << e << endl;
    throw;
  }
}


//! Applies gauge transformation matrices on a propagator
/*! \defgroup qpropgfix Gauge fix a propagator
 *  \ingroup main
 *
 * Main program for gauge fixing a propagator
 */

int main(int argc, char *argv[])
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);

  START_CODE();
  
  // Parameter structure for the input
  QpropGFix_input_t input;

  // Instantiate xml reader for DATA
  XMLReader xml_in(Chroma::getXMLInputFileName());

  // Read data
  read(xml_in, "/qpropgfix", input);

  // Setup QDP
  Layout::setLattSize(input.param.nrow);
  Layout::create();

  QDPIO::cout << "QPROPGFIX: propagator gauge fixing utility" << endl;

  XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();
  push(xml_out, "qpropgfix");

  proginfo(xml_out);    // Print out basic program info

  write(xml_out, "input", xml_in); // save a copy of the input
  xml_out.flush();
  
  /*
   * Now read them thangs...
   */
  /*
   *  Read in the gauge fixed configuration along with relevant information.
   */
  multi1d<LatticeColorMatrix> u(Nd);
  XMLReader gauge_file_xml, gauge_xml;

  // Startup gauge
  gaugeStartup(gauge_file_xml, gauge_xml, u, input.cfg);

  /*
   * Read in the gauge transformation matrices
   */
  LatticeColorMatrix  g;
  XMLReader transf_file_xml, transf_xml;
  {
    QDPFileReader from(transf_xml, input.prop.gfix_in_file, QDPIO_SERIAL);

    LatticeColorMatrixF g_f;
    read(from,transf_xml,g_f);         // Always save in single precision!
    g = g_f;

    close(from);
  }


  /*
   * Read in a Chroma prop
   */
  LatticePropagator  prop;
  XMLReader prop_in_xml, prop_in_file_xml;

  push(xml_out,"SciDAC_propagator");
  write(xml_out, "prop_in_file", input.prop.prop_in_file);

  readQprop(prop_in_file_xml, prop_in_xml, prop, 
	    input.prop.prop_in_file, QDPIO_SERIAL);

  write(xml_out, "File_xml", prop_in_file_xml);
  write(xml_out, "Record_xml", prop_in_xml);
  pop(xml_out);
    

  // Try to invert this record XML into a source struct
  // Also pull out the id of this source
  ChromaProp_t prop_header;
  PropSourceConst_t source_header;
  try
  {
    read(prop_in_xml, "/Propagator/ForwardProp", prop_header);
    read(prop_in_xml, "/Propagator/PropSource", source_header);
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error extracting forward_prop header: " << e << endl;
    throw;
  }

  // Derived from input prop
  int j_decay = source_header.j_decay;
  multi1d<int> t_srce = source_header.getTSrce();


  // Sanity check - write out the propagator (pion) correlator in the j_decay direction
  {
    // Initialize the slow Fourier transform phases
    SftMom phases(0, true, j_decay);

    multi1d<Double> prop_corr = sumMulti(localNorm2(prop), 
					 phases.getSet());

    push(xml_out, "Prop_correlator");
    write(xml_out, "prop_corr", prop_corr);
    pop(xml_out);
  }

  xml_out.flush();


  /*
   * Gauge transform the beasty
   */
  {
    LatticePropagator tmp = g * prop * adj(peekSite(g,t_srce));
    prop = tmp;
  }



  // Sanity check - write out the propagator (pion) correlator in the j_decay direction
  {
    // Initialize the slow Fourier transform phases
    SftMom phases(0, true, j_decay);

    multi1d<Double> gfix_prop_corr = sumMulti(localNorm2(prop), 
					      phases.getSet());

    push(xml_out, "GFixProp_correlator");
    write(xml_out, "gfix_prop_corr", gfix_prop_corr);
    pop(xml_out);
  }


  /*
   * Now write them thangs...
   */ 
  {
    XMLBufferWriter prop_out_file_xml;
    push(prop_out_file_xml, "propagator");
    int id = 0;    // NEED TO FIX THIS - SOMETHING NON-TRIVIAL NEEDED
    write(prop_out_file_xml, "id", id);
    pop(prop_out_file_xml);

    XMLBufferWriter prop_out_record_xml;
    push(prop_out_record_xml, "Propagator");
    write(prop_out_record_xml, "ForwardProp", prop_header);
    write(prop_out_record_xml, "PropSource", source_header);
    write(prop_out_record_xml, "Config_info", gauge_xml);
    pop(prop_out_record_xml);
    
    // Write the source
    writeQprop(prop_out_file_xml, prop_out_record_xml, prop,
	       input.prop.prop_out_file, input.prop.prop_out_volfmt, 
	       QDPIO_SERIAL);
  }

  pop(xml_out);   // qpropgfix
        
  END_CODE();

  // Time to bolt
  Chroma::finalize();

  exit(0);
}
