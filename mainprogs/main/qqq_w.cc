// $Id: qqq_w.cc,v 1.10 2004-02-23 03:13:58 edwards Exp $
/*! \file
 *  \brief Main code for generalized quark propagator
 *
 *  This is the test program for the routine that reads in a quark propagator,
 *  stored in SZIN format, and computes the generalised quark propagators
 *  that will be the basis of our baryon code
 */

#include <iostream>
#include <cstdio>

#include "chroma.h"

using namespace QDP;

/*
 *  This is the reader for the input parameters
 */

/*
 * Input 
 */
// Parameters which must be determined from the XML input
// and written to the XML output
struct Param_t
{
  multi1d<int> nrow;		// Lattice dimension
};

//! Propagators
struct Prop_t
{
  string           prop_file;  // The file is expected to be in SciDAC format!
};

//! Barcomp info
struct Barcomp_t
{
  string           qqq_file;  // The file is expected to be in SciDAC format!
};

//! Mega-structure of all input
struct QQQ_input_t
{
  Param_t          param;
  Cfg_t            cfg;
  Prop_t           prop;
  Barcomp_t        barcomp;
};



//! Input propagator file
void read(XMLReader& xml, const string& path, Prop_t& input)
{
  XMLReader inputtop(xml, path);

  read(inputtop, "prop_file", input.prop_file);
}


//! Output barcomp file
void read(XMLReader& xml, const string& path, Barcomp_t& input)
{
  XMLReader inputtop(xml, path);

  read(inputtop, "qqq_file", input.qqq_file);
}


//! Parameters for running code
void read(XMLReader& xml, const string& path, Param_t& param)
{
  XMLReader paramtop(xml, path);

  int version;
  read(paramtop, "version", version);

  switch (version) 
  {
  case 2:
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
void read(XMLReader& xml, const string& path, QQQ_input_t& input)
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

    // Read in the barcomp file info
    read(inputtop, "Barcomp", input.barcomp);
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading qqq data: " << e << endl;
    throw;
  }
}



int main(int argc, char **argv)
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  // Parameter structure for the input
  QQQ_input_t input;

  // Instantiate xml reader for DATA
  XMLReader xml_in("DATA");

  // Read data
  read(xml_in, "/qqq_w", input);

  // Specify lattice size, shape, etc.
  Layout::setLattSize(input.param.nrow);
  Layout::create();

  QDPIO::cout << " QQQ: Generalized propagator generation" << endl;


  // Read a gauge field
  multi1d<LatticeColorMatrix> u(Nd);
  XMLReader gauge_xml;

  switch (input.cfg.cfg_type) 
  {
  case CFG_TYPE_SZIN:
    readSzin(gauge_xml, u, input.cfg.cfg_file);
    break;

  case CFG_TYPE_NERSC:
    readArchiv(gauge_xml, u, input.cfg.cfg_file);
    break;
  default :
    QDP_error_exit("Configuration type is unsupported.");
  }

  QDPIO::cout << "Gauge field read!" << endl;

  // Check if the gauge field configuration is unitarized
  unitarityCheck(u);


  // Output
  XMLFileWriter xml_out("XMLDAT");
  push(xml_out,"qqq");

  xml_out << xml_in;  // save a copy of the input
  write(xml_out, "config_info", gauge_xml);
  xml_out.flush();


  // Calculate some gauge invariant observables just for info.
  Double w_plaq, s_plaq, t_plaq, link;
  MesPlq(u, w_plaq, s_plaq, t_plaq, link);
  
  push(xml_out, "Observables");
  write(xml_out, "w_plaq", w_plaq);
  write(xml_out, "s_plaq", s_plaq);
  write(xml_out, "t_plaq", t_plaq);
  write(xml_out, "link", link);


  /*
   * Read the quark propagator and extract headers
   *
   * For now, only 1 propagator is supported.
   */
  XMLReader prop_file_xml, prop_record_xml;
  LatticePropagator quark_propagator;
  ChromaProp_t prop_header;
  PropSource_t source_header;
  PropSink_t   sink_header;
  {
    QDPIO::cout << "Attempt to read forward propagator" << endl;
    readQprop(prop_file_xml, 
	      prop_record_xml, quark_propagator,
	      input.prop.prop_file, QDPIO_SERIAL);
    QDPIO::cout << "Forward propagator successfully read" << endl;
   
    // Try to invert this record XML into a ChromaProp struct
    // Also pull out the id of this source
    try
    {
      read(prop_record_xml, "/SinkSmear/PropSink", sink_header);
      read(prop_record_xml, "/SinkSmear/ForwardProp", prop_header);
      read(prop_record_xml, "/SinkSmear/PropSource", source_header);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << "Error extracting forward_prop header: " << e << endl;
      throw;
    }
  }

  // Save prop input
  write(xml_out, "PropSource", source_header);
  write(xml_out, "ForwardProp", prop_header);
  write(xml_out, "PropSink", sink_header);

  // Derived from input prop
  int j_decay = source_header.j_decay;
  multi1d<int> boundary = prop_header.boundary;
  multi1d<int> t_source = source_header.t_source;
  int t0      = t_source[j_decay];
  int bc_spec = boundary[j_decay];

  // Initialize the slow Fourier transform phases
  SftMom phases(0, true, j_decay);

  // Sanity check - write out the propagator (pion) correlator in the j_decay direction
  {
    multi1d<Double> prop_corr = sumMulti(localNorm2(quark_propagator), 
					 phases.getSet());

    push(xml_out, "SinkSmearedProp_correlator");
    write(xml_out, "sink_smeared_prop_corr", prop_corr);
    pop(xml_out);
  }

  /*
   * Generalized propagator calculation
   */
  multiNd<Complex> barprop;

  //
  // quark_propagator = quark_prop * 1.0e10;
  //
  barcomp(barprop,
	  quark_propagator,
	  quark_propagator,
	  quark_propagator,
	  phases, t0, bc_spec);

  // Convert the data into a mult1d
  multi1d<Complex> barprop_1d;
  convertBarcomp(barprop_1d, barprop, j_decay);

  // Save the qqq output
  // ONLY SciDAC output format is supported!
  {
    XMLBufferWriter file_xml;
    push(file_xml, "qqq");
    int id = 0;    // NEED TO FIX THIS - SOMETHING NON-TRIVIAL NEEDED
    write(file_xml, "id", id);
    pop(file_xml);

    XMLBufferWriter record_xml;
    push(record_xml, "QQQ");

    // For the moment, the single propagator header is replicated.
    // However, 3 distinct propagators are supported
    push(record_xml, "Propagator1");
    write(record_xml, "PropSink", sink_header);
    write(record_xml, "ForwardProp", prop_header);
    write(record_xml, "PropSource", source_header);
    pop(record_xml);

    push(record_xml, "Propagator2");
    write(record_xml, "PropSink", sink_header);
    write(record_xml, "ForwardProp", prop_header);
    write(record_xml, "PropSource", source_header);
    pop(record_xml);

    push(record_xml, "Propagator3");
    write(record_xml, "PropSink", sink_header);
    write(record_xml, "ForwardProp", prop_header);
    write(record_xml, "PropSource", source_header);
    pop(record_xml);

    write(record_xml, "Config_info", gauge_xml);
    pop(record_xml);

    // Write the scalar data
    QDPFileWriter to(file_xml, input.barcomp.qqq_file, 
		     QDPIO_SINGLEFILE, QDPIO_SERIAL, QDPIO_OPEN);
//   write(to,record_xml,barprop_1d);    // OOOPS, need to make a function in QDP !!!
    QDPIO::cout << "OOOOPSS, need a new QDP function here" << endl;
    close(to);
  }

  pop(xml_out);    // qqq

  xml_out.close();
  xml_in.close();

  // Time to bolt
  QDP_finalize();

  END_CODE("seqprop");

  exit(0);
}

