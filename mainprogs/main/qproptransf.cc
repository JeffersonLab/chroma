// $Id: qproptransf.cc,v 1.1 2004-03-04 03:20:50 edwards Exp $
/*! \file
 *  \brief Converts quark propagators in one format into another format.
 */

#include "chroma.h"

using namespace QDP;

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
  PropType  prop_in_type;       // propagator format
  string    prop_in_file;

  PropType  prop_out_type;      // propagator format
  string    prop_out_file;
  QDP_volfmt_t prop_out_volfmt; // volume format (SINGLEFILE or MULTIFILE)
};

struct QpropTransf_input_t
{
  Param_t  param;
  Cfg_t    cfg;
  Prop_t   prop;
};



//! Propagator parameters
void read(XMLReader& xml, const string& path, Prop_t& input)
{
  XMLReader inputtop(xml, path);

  read(inputtop, "prop_in_type", input.prop_in_type);
  read(inputtop, "prop_in_file", input.prop_in_file);

  read(inputtop, "prop_out_type", input.prop_out_type);
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
void read(XMLReader& xml, const string& path, QpropTransf_input_t& input)
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


//! Many-to-many propagator transformation routine
/*! \defgroup qproptransf Tranformation routine
 *  \ingroup main
 *
 * Main program for transforming propagator formats
 */

int main(int argc, char *argv[])
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  START_CODE("qproptransf");
  
  // Parameter structure for the input
  QpropTransf_input_t input;

  // Hmm, at this moment I've not decided between what style of input
  // to use - all XML or mostly ascii and some XML when needed.
#if 0
  // Instantiate xml reader for DATA
  XMLReader xml_in("DATA");

  // Read data
  read(xml_in, "/qproptransf", input);
#else
  input.param.nrow.resize(Nd);

  QDPIO::cout << "Enter lattice size\n";
  QDPIO::cin >> input.param.nrow;
#endif

  // Setup QDP
  Layout::setLattSize(input.param.nrow);
  Layout::create();

  XMLFileWriter xml_out("qproptransf.xml");
  push(xml_out, "qproptransf");

  proginfo(xml_out);    // Print out basic program info

//  xml_out << xml_in;  // save a copy of the input
//  write(xml_out, "config_info", gauge_xml);
  xml_out.flush();


  int input_type;
  QDPIO::cout << "Enter input propagator format\n"
	      << "  (1) SZIN\n"
	      << "  (2) SciDAC\n"
	      << "  (3) Kentucky\n";
  QDPIO::cin >> input_type;
  
  int output_type;
  QDPIO::cout << "Enter output propagator field type\n"
	      << "  (1) SZIN\n"
	      << "  (2) SciDAC\n";
  QDPIO::cin >> output_type;
  
  string prop_in_file;
  QDPIO::cout << "Enter input file name\n";
  QDPIO::cin >> prop_in_file;

  string prop_out_file;
  QDPIO::cout << "Enter output file name\n";
  QDPIO::cin >> prop_out_file;
  
  /*
   * Now read them thangs...
   */
  XMLReader prop_in_xml, prop_in_file_xml;
  LatticePropagator  prop;

  switch (input_type)
  {
  case 1:
    // SZIN
    push(xml_out,"SZIN_propagator");
    write(xml_out, "input_type", input_type);
    write(xml_out, "prop_in_file", prop_in_file);

    readSzinQprop(prop_in_xml, prop, prop_in_file);

    write(xml_out, "propagator_info", prop_in_xml);
    pop(xml_out);
    break;

  case 2:
    // SciDAC
    push(xml_out,"SciDAC_propagator");
    write(xml_out, "input_type", input_type);
    write(xml_out, "prop_in_file", prop_in_file);

    readQprop(prop_in_file_xml, prop_in_xml, prop, 
	      prop_in_file, QDPIO_SERIAL);

    write(xml_out, "File_xml", prop_in_file_xml);
    write(xml_out, "Record_xml", prop_in_xml);
    pop(xml_out);
    break;

  default:
    QDP_error_exit("unknown input type", input_type);
  }
    

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
  switch (output_type)
  {
  case 1:
  {
    // SZIN
    Real Kappa(3.14159265359);

    push(xml_out,"SZIN_propagator");
    write(xml_out, "output_type", output_type);
    write(xml_out, "prop_out_file", prop_out_file);
    pop(xml_out);

    writeSzinQprop(prop, prop_out_file, Kappa);
  }
  break;

  case 2:
    // SciDAC
    QDPIO::cerr << "SciDAC output not implemented\n";
    QDP_abort(1);
    break;

  default:
    QDPIO::cerr << "unknown output type = " << output_type << endl;
    QDP_abort(1);
  }

  pop(xml_out);
        
  xml_out.close();
//  xml_in.close();

  END_CODE("qproptransf");

  // Time to bolt
  QDP_finalize();

  exit(0);
}
