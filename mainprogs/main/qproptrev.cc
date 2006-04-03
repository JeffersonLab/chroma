// $Id: qproptrev.cc,v 3.0 2006-04-03 04:59:14 edwards Exp $
/*! \file
 *  \brief Time-reverse a propagator
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

struct QpropTRev_input_t
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
void read(XMLReader& xml, const string& path, QpropTRev_input_t& input)
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


//! Applies gauge transformation matrices on a propagator
/*! \defgroup qproptrev Tranformation routine
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
  QpropTRev_input_t input;

  // Instantiate xml reader for DATA
  XMLReader xml_in(Chroma::getXMLInputFileName());

  // Read data
  read(xml_in, "/qproptrev", input);

  // Setup QDP
  Layout::setLattSize(input.param.nrow);
  Layout::create();

  QDPIO::cout << "QPROPTREV: propagator gauge fixing utility" << endl;

  XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();
  push(xml_out, "qproptrev");

  proginfo(xml_out);    // Print out basic program info

  write(xml_out, "input", xml_in); // save a copy of the input
  xml_out.flush();
  
  /*
   * Now read them thangs...
   */
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
  int j_decay  = source_header.j_decay;
  int t_source = source_header.t_source;


  // Initialize the slow Fourier transform phases
  // This is used to get the time-slice subsets in the j_decay dir
  SftMom phases(0, true, j_decay);

  // Length of lattice in j_decay direction and 3pt correlations fcns
  int length = phases.numSubsets();

  // Sanity check - write out the propagator (pion) correlator in the j_decay direction
  {
    multi1d<Double> prop_corr = sumMulti(localNorm2(prop), 
					 phases.getSet());

    push(xml_out, "Prop_correlator");
    write(xml_out, "prop_corr", prop_corr);
    pop(xml_out);
  }

  xml_out.flush();


  /*
   * Charge-conjugation and time-reverse the beasty
   */
  {
    /* Time-charge reverse the quark propagators */
    /* S_{CT} = gamma_5 gamma_4 = gamma_1 gamma_2 gamma_3 = Gamma(7) */
    LatticePropagator prop_tmp = - (Gamma(7) * prop * Gamma(7));

    // This is a really dumb way to implement this. Shift each slice around.
    // Do nothing on the t=0 slice
    prop[phases.getSet()[0]] = prop_tmp;

    LatticePropagator tmp1, tmp2, tmp3;

    for(int t=1; t < length; ++t)
    {
      int tp = (length - t) % length;

      if (t < tp)
      {
	tmp1[phases.getSet()[tp]] = prop_tmp;
	for(int k=tp-1; k >= t; --k)
	{
	  tmp2[phases.getSet()[k]] = shift(tmp1, FORWARD, j_decay);
	  tmp1[phases.getSet()[k]] = tmp2;
	}

	prop[phases.getSet()[t]] = tmp1;
      }
      else if (t == tp)
      {
	prop[phases.getSet()[t]] = prop_tmp;
      }
      else if (t > tp)
      {
	tmp1[phases.getSet()[tp]] = prop_tmp;
	for(int k=tp+1; k <= t; ++k)
	{
	  tmp2[phases.getSet()[k]] = shift(tmp1, BACKWARD, j_decay);
	  tmp1[phases.getSet()[k]] = tmp2;
	}

	prop[phases.getSet()[t]] = tmp1;
      }
    }
  }



  // Sanity check - write out the propagator (pion) correlator in the j_decay direction
  {
    multi1d<Double> trev_prop_corr = sumMulti(localNorm2(prop), 
					      phases.getSet());

    push(xml_out, "TRevProp_correlator");
    write(xml_out, "trev_prop_corr", trev_prop_corr);
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

    source_header.t_source = (length - source_header.t_source) % length;

    XMLBufferWriter prop_out_record_xml;
    push(prop_out_record_xml, "Propagator");
    write(prop_out_record_xml, "ForwardProp", prop_header);
    write(prop_out_record_xml, "PropSource", source_header);
    {
      QDPIO::cout << "Create config info" << endl;
      XMLReader gauge_xml(prop_in_xml, "/Propagator/Config_info");
      ostringstream gauge_str;
      gauge_xml.print(gauge_str);
      write(prop_out_record_xml, "Config_info", gauge_str);
      QDPIO::cout << "Done config info" << endl;
    }
    pop(prop_out_record_xml);
    
    // Write the source
    writeQprop(prop_out_file_xml, prop_out_record_xml, prop,
	       input.prop.prop_out_file, input.prop.prop_out_volfmt, 
	       QDPIO_SERIAL);
  }

  pop(xml_out);   // qproptrev
        
  END_CODE();

  // Time to bolt
  Chroma::finalize();

  exit(0);
}
