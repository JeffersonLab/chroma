// $Id: inline_qqq_w.cc,v 1.1 2005-04-06 04:34:53 edwards Exp $
/*! \file
 * \brief Inline construction of qqq_w
 *
 * QQQ calcs
 */

#include "meas/inline/hadron/inline_qqq_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "meas/hadron/barcomp_w.h"
#include "io/barcomp_io.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "util/ferm/diractodr.h"

namespace Chroma 
{ 
  namespace InlineQQQEnv 
  { 
    AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					    const std::string& path) 
    {
      return new InlineQQQ(InlineQQQParams(xml_in, path));
    }

    const std::string name = "QQQ";
    const bool registered = TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
  };


  //! Param input
  void read(XMLReader& xml, const string& path, InlineQQQParams::Param_t& input)
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

    read(paramtop, "Dirac_basis", input.Dirac_basis);
    read(paramtop, "nrow", input.nrow);
  }

  //! Param output
  void write(XMLWriter& xml, const string& path, const InlineQQQParams::Param_t& input)
  {
    push(xml, path);

    int version = 4;
    write(xml, "version", version);
    write(xml, "Dirac_basis", input.Dirac_basis);
    write(xml, "nrow", input.nrow);

    pop(xml);
  }

  //! Propagator input
  void read(XMLReader& xml, const string& path, InlineQQQParams::Prop_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "prop_file", input.prop_file);
    read(inputtop, "qqq_file", input.qqq_file);
  }

  //! Propagator output
  void write(XMLWriter& xml, const string& path, const InlineQQQParams::Prop_t& input)
  {
    push(xml, path);

    write(xml, "prop_file", input.prop_file);
    write(xml, "qqq_file", input.qqq_file);

    pop(xml);
  }


  // Param stuff
  InlineQQQParams::InlineQQQParams() { frequency = 0; }

  InlineQQQParams::InlineQQQParams(XMLReader& xml_in, const std::string& path) 
  {
    try 
    {
      XMLReader paramtop(xml_in, path);

      if (paramtop.count("Frequency") == 1)
	read(paramtop, "Frequency", frequency);
      else
	frequency = 1;

      // Parameters for source construction
      read(paramtop, "Param", param);

      // Read in the output propagator/source configuration info
      read(paramtop, "Prop", prop);
    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << "Caught Exception reading XML: " << e << endl;
      QDP_abort(1);
    }
  }


  void
  InlineQQQParams::write(XMLWriter& xml_out, const std::string& path) 
  {
    push(xml_out, path);
    
    // Parameters for source construction
    Chroma::write(xml_out, "Param", param);

    // Write out the output propagator/source configuration info
    Chroma::write(xml_out, "Prop", prop);

    pop(xml_out);
  }


  void 
  InlineQQQ::operator()(const multi1d<LatticeColorMatrix>& u,
			XMLBufferWriter& gauge_xml,
			unsigned long update_no,
			XMLWriter& xml_out) 
  {
    START_CODE();

    push(xml_out,"qqq");
    write(xml_out, "update_no", update_no);

    QDPIO::cout << " QQQ: Generalized propagator generation" << endl;

    // Write out the input
    params.write(xml_out, "Input");

    proginfo(xml_out);    // Print out basic program info

    // Write out the input
    params.write(xml_out, "Input");

    // Write out the config info
    write(xml_out, "Config_info", gauge_xml);

    push(xml_out, "Output_version");
    write(xml_out, "out_version", 4);
    pop(xml_out);

    // Calculate some gauge invariant observables just for info.
    MesPlq(xml_out, "Observables", u);

    // Check to make sure there are 3 files
    const int Nprops = 3;
    if (params.prop.prop_file.size() != Nprops)
    {
      QDPIO::cerr << "Error on input params - expecting 3 filenames" << endl;
      QDP_abort(1);
    }

    write(xml_out, "propFiles", params.prop.prop_file);
    write(xml_out, "barcompFile", params.prop.qqq_file);

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
      if ((i == 0) || (params.prop.prop_file[i] != params.prop.prop_file[0]))
      {
	XMLReader prop_file_xml, prop_record_xml;

	QDPIO::cout << "Attempt to read forward propagator XX" 
		    << params.prop.prop_file[i] << "XX" << endl;
	readQprop(prop_file_xml, 
		  prop_record_xml, quark_propagator[i],
		  params.prop.prop_file[i], QDPIO_SERIAL);
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
    write(xml_out, "Propagator_input", qqq);

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

      push(xml_out, "SinkSmearedProp_correlator");
      write(xml_out, "correlator_num", i);
      write(xml_out, "sink_smeared_prop_corr", prop_corr);
      pop(xml_out);
    }

    /*
     * Generalized propagator calculation
     */
    multiNd<Complex> barprop;

    // Switch to Dirac-basis if desired.
    if (params.param.Dirac_basis)
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
      QDPFileWriter to(file_xml, params.prop.qqq_file, 
		       QDPIO_SINGLEFILE, QDPIO_SERIAL, QDPIO_OPEN);
      write(to,record_xml,barprop_1d);
      close(to);
    }

    pop(xml_out);    // qqq

    END_CODE();
  } 

};
