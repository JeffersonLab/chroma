// $Id: qqq_w.cc,v 1.23 2005-02-27 22:47:19 edwards Exp $
/*! \file
 *  \brief Main code for generalized quark propagator
 *
 *  This is the test program for the routine that reads in a quark propagator,
 *  stored in SZIN format, and computes the generalised quark propagators
 *  that will be the basis of our baryon code
 */

#include "chroma.h"

using namespace Chroma;

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
  bool         Dirac_basis;     // Use the Dirac basis for output?
  multi1d<int> nrow;		// Lattice dimension
};

//! Propagators and Barcomp info
struct Prop_t
{
  multi1d<string>  prop_file;  // The file is expected to be in SciDAC format!
  string           qqq_file;   // The file is expected to be in SciDAC format!
};

//! Sink-smear
struct QQQParam_t
{
  Param_t          qqq_param;
  Prop_t           prop;
};

//! Mega-structure of all input
struct QQQ_input_t
{
  multi1d<QQQParam_t>  param;
  Cfg_t                cfg;
};


//! Parameters for running code
void read(XMLReader& xml, const string& path, Param_t& param)
{
  XMLReader paramtop(xml, path);

  int version;
  read(paramtop, "version", version);

  switch (version) 
  {
  case 4:
    break;

  default:
    /**************************************************************************/
    QDPIO::cerr << "Input parameter version " << version << " unsupported." << endl;
    QDP_abort(1);
  }

  read(paramtop, "Dirac_basis", param.Dirac_basis);
  read(paramtop, "nrow", param.nrow);
}


//! Input propagator files and output barcomp file
void read(XMLReader& xml, const string& path, Prop_t& input)
{
  XMLReader inputtop(xml, path);

  read(inputtop, "prop_file", input.prop_file);
  read(inputtop, "qqq_file", input.qqq_file);
}


//! Input propagator files and output barcomp file
void read(XMLReader& xml, const string& path, QQQParam_t& input)
{
  XMLReader inputtop(xml, path);

  read(inputtop, "QQQParam", input.qqq_param);
  read(inputtop, "Prop", input.prop);
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
  ChromaInitialize(&argc, &argv);

  START_CODE();

  // Parameter structure for the input
  QQQ_input_t input;

  // Instantiate xml reader for DATA
  XMLReader xml_in("./DATA");

  // Read data
  read(xml_in, "/qqq_w", input);

  if (input.param.size() == 0)
    QDP_error_exit("QQ input empty");

  // Specify lattice size, shape, etc.
  Layout::setLattSize(input.param[0].qqq_param.nrow);
  Layout::create();

  QDPIO::cout << " QQQ: Generalized propagator generation" << endl;

  // Read a gauge field
  multi1d<LatticeColorMatrix> u(Nd);
  XMLReader gauge_file_xml, gauge_xml;

  // Start up the gauge field
  gaugeStartup(gauge_file_xml, gauge_xml, u, input.cfg);

  QDPIO::cout << "Gauge field read!" << endl;

  // Check if the gauge field configuration is unitarized
  unitarityCheck(u);


  // Output
  XMLFileWriter& xml_out = TheXMLOutputWriter::Instance();
  push(xml_out,"qqq");

  proginfo(xml_out);    // Print out basic program info

  // Write out the input
  write(xml_out, "Input", xml_in);

  // Write out the config info
  write(xml_out, "Config_info", gauge_xml);

  push(xml_out, "Output_version");
  write(xml_out, "out_version", 4);
  pop(xml_out);

  // Calculate some gauge invariant observables just for info.
  Double w_plaq, s_plaq, t_plaq, link;
  MesPlq(u, w_plaq, s_plaq, t_plaq, link);
  
  push(xml_out, "Observables");
  write(xml_out, "w_plaq", w_plaq);
  write(xml_out, "s_plaq", s_plaq);
  write(xml_out, "t_plaq", t_plaq);
  write(xml_out, "link", link);

  xml_out.flush();

  // Check the barcomp size
  if (input.param.size() == 0)
  {
    QDPIO::cerr << "Expecting some propagator and barcomp files" << endl;
    QDP_abort(1);
  }

  //
  // Big loop over all the barcomp input
  //
  XMLArrayWriter xml_array(xml_out,input.param.size());
  push(xml_array, "Barcomp_measurements");

  for(int loop = 0; loop < xml_array.size(); ++loop)
  {
    push(xml_array);
    write(xml_array, "loop", loop);

    QDPIO::cout << "Barcomp loop = " << loop << endl;

    const QQQParam_t&  param = input.param[loop];  // make a short-cut reference

    // Check to make sure there are 3 files
    const int Nprops = 3;
    if (param.prop.prop_file.size() != Nprops)
    {
      QDPIO::cerr << "Error on input params - expecting 3 filenames" << endl;
      QDP_abort(1);
    }

    write(xml_array, "propFiles", param.prop.prop_file);
    write(xml_array, "barcompFile", param.prop.qqq_file);

    /*
     * Read the quark propagators and extract headers
     */
    multi1d<LatticePropagator> quark_propagator(Nprops);
    QQQBarcomp_t  qqq;
    qqq.Dirac_basis = false;
    qqq.forward_props.resize(Nprops);
    for(int i=0; i < Nprops; ++i)
    {
      // Optimize the read - if the filename is the same, no need to read
      if ((i == 0) || (param.prop.prop_file[i] != param.prop.prop_file[0]))
      {
	XMLReader prop_file_xml, prop_record_xml;

	QDPIO::cout << "Attempt to read forward propagator XX" 
		    << param.prop.prop_file[i] << "XX" << endl;
	readQprop(prop_file_xml, 
		  prop_record_xml, quark_propagator[i],
		  param.prop.prop_file[i], QDPIO_SERIAL);
	QDPIO::cout << "Forward propagator successfully read" << endl;
   
	// Try to invert this record XML into a ChromaProp struct
	// Also pull out the id of this source
	try
	{
	  read(prop_record_xml, "/SinkSmear/PropSink", qqq.forward_props[i].sink_header);
	  read(prop_record_xml, "/SinkSmear/ForwardProp", qqq.forward_props[i].prop_header);
	  read(prop_record_xml, "/SinkSmear/PropSource", qqq.forward_props[i].source_header);
	}
	catch (const string& e) 
	{
	  QDPIO::cerr << "Error extracting forward_prop header: " << e << endl;
	  QDP_abort(1);
	}
      }
      else
      {
	QDPIO::cout << "Same prop: optimize away read of propagator "  << i << endl;
	quark_propagator[i] = quark_propagator[0];
	qqq.forward_props[i].sink_header = qqq.forward_props[0].sink_header;
	qqq.forward_props[i].prop_header = qqq.forward_props[0].prop_header;
	qqq.forward_props[i].source_header = qqq.forward_props[0].source_header;
      }
    }

    // Save prop input
    write(xml_array, "Propagator_input", qqq);

    // Derived from input prop
    int j_decay = qqq.forward_props[0].source_header.j_decay;
    multi1d<int> boundary = getFermActBoundary(qqq.forward_props[0].prop_header.fermact);
    multi1d<int> t_source = qqq.forward_props[0].source_header.t_source;
    int t0      = t_source[j_decay];
    int bc_spec = boundary[j_decay];

    // Initialize the slow Fourier transform phases
    SftMom phases(0, true, j_decay);

    // Sanity check - write out the propagator (pion) correlator in the j_decay direction
    for(int i=0; i < Nprops; ++i)
    {
      multi1d<Double> prop_corr = sumMulti(localNorm2(quark_propagator[i]), 
					   phases.getSet());

      push(xml_array, "SinkSmearedProp_correlator");
      write(xml_array, "correlator_num", i);
      write(xml_array, "sink_smeared_prop_corr", prop_corr);
      pop(xml_array);
    }

    /*
     * Generalized propagator calculation
     */
    multiNd<Complex> barprop;

    // Switch to Dirac-basis if desired.
    if (param.qqq_param.Dirac_basis)
    {
      qqq.Dirac_basis = true;

      // The spin basis matrix
      SpinMatrix U = DiracToDRMat();

      LatticePropagator q_tmp;
      for(int i=0; i < Nprops; ++i)
      {
	q_tmp = adj(U) * quark_propagator[i] * U;   // DeGrand-Rossi ---> Dirac
	quark_propagator[i] = q_tmp;
      }
    }

    // Compute generation propagator
    barcomp(barprop,
	    quark_propagator[0],
	    quark_propagator[1],
	    quark_propagator[2],
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
      write(record_xml, ".", qqq);  // do not write the outer group
      write(record_xml, "Config_info", gauge_xml);
      pop(record_xml);  // QQQ

      // Write the scalar data
      QDPFileWriter to(file_xml, param.prop.qqq_file, 
		       QDPIO_SINGLEFILE, QDPIO_SERIAL, QDPIO_OPEN);
      write(to,record_xml,barprop_1d);
      close(to);
    }

    pop(xml_array);  // array element

  } // end for(loop)

  pop(xml_array);  // Barcomp_measurements
  pop(xml_out);    // qqq

  END_CODE();

  // Time to bolt
  ChromaFinalize();



  exit(0);
}

