// $Id: inline_qqbar_w.cc,v 1.1 2005-08-31 05:50:00 edwards Exp $
/*! \file
 * \brief Inline construction of qqbar
 *
 * QQbar calcs
 */

#include "meas/inline/hadron/inline_qqbar_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "meas/hadron/mescomp_w.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "util/ferm/diractodr.h"

namespace Chroma 
{ 
  namespace InlineQQbarEnv 
  { 
    AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					    const std::string& path) 
    {
      return new InlineQQbar(InlineQQbarParams(xml_in, path));
    }

    const std::string name = "QQBAR";
    const bool registered = TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
  };


  //! Param input
  void read(XMLReader& xml, const string& path, InlineQQbarParams::Param_t& input)
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
  void write(XMLWriter& xml, const string& path, const InlineQQbarParams::Param_t& input)
  {
    push(xml, path);

    int version = 4;
    write(xml, "version", version);
    write(xml, "Dirac_basis", input.Dirac_basis);
    write(xml, "nrow", input.nrow);

    pop(xml);
  }

  //! Propagator input
  void read(XMLReader& xml, const string& path, InlineQQbarParams::Prop_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "prop_file", input.prop_file);
    read(inputtop, "qqbar_file", input.qqbar_file);
  }

  //! Propagator output
  void write(XMLWriter& xml, const string& path, const InlineQQbarParams::Prop_t& input)
  {
    push(xml, path);

    write(xml, "prop_file", input.prop_file);
    write(xml, "qqbar_file", input.qqbar_file);

    pop(xml);
  }


  // Param stuff
  InlineQQbarParams::InlineQQbarParams() { frequency = 0; }

  InlineQQbarParams::InlineQQbarParams(XMLReader& xml_in, const std::string& path) 
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
  InlineQQbarParams::write(XMLWriter& xml_out, const std::string& path) 
  {
    push(xml_out, path);
    
    // Parameters for source construction
    Chroma::write(xml_out, "Param", param);

    // Write out the output propagator/source configuration info
    Chroma::write(xml_out, "Prop", prop);

    pop(xml_out);
  }


  void 
  InlineQQbar::operator()(const multi1d<LatticeColorMatrix>& u,
			XMLBufferWriter& gauge_xml,
			unsigned long update_no,
			XMLWriter& xml_out) 
  {
    START_CODE();

    push(xml_out,"qqbar");
    write(xml_out, "update_no", update_no);

    QDPIO::cout << " QQbar: Generalized propagator generation" << endl;

    // Write out the input
    params.write(xml_out, "Input");

    proginfo(xml_out);    // Print out basic program info

    // Write out the input
    params.write(xml_out, "Input");

    // Write out the config info
    write(xml_out, "Config_info", gauge_xml);

    push(xml_out, "Output_version");
    write(xml_out, "out_version", 1);
    pop(xml_out);

    // Calculate some gauge invariant observables just for info.
    MesPlq(xml_out, "Observables", u);

    // Check to make sure there are 2 files
    const int Nprops = 2;
    if (params.prop.prop_file.size() != Nprops)
    {
      QDPIO::cerr << "Error on input params - expecting 2 filenames" << endl;
      QDP_abort(1);
    }

    write(xml_out, "propFiles", params.prop.prop_file);
    write(xml_out, "mescompFile", params.prop.qqbar_file);

    /*
     * Read the quark propagators and extract headers
     */
    multi1d<LatticePropagator> quark_propagator(Nprops);
    QQbarMescomp_t  qqbar;
    qqbar.Dirac_basis = false;
    qqbar.forward_props.resize(Nprops);
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
	  read(prop_record_xml, "/SinkSmear/PropSink", qqbar.forward_props[i].sink_header);
	  read(prop_record_xml, "/SinkSmear/ForwardProp", qqbar.forward_props[i].prop_header);
	  read(prop_record_xml, "/SinkSmear/PropSource", qqbar.forward_props[i].source_header);
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
	qqbar.forward_props[i].sink_header = qqbar.forward_props[0].sink_header;
	qqbar.forward_props[i].prop_header = qqbar.forward_props[0].prop_header;
	qqbar.forward_props[i].source_header = qqbar.forward_props[0].source_header;
      }
    }

    // Save prop input
    write(xml_out, "Propagator_input", qqbar);

    // Derived from input prop
    int j_decay = qqbar.forward_props[0].source_header.j_decay;
    multi1d<int> boundary = getFermActBoundary(qqbar.forward_props[0].prop_header.fermact);
    multi1d<int> t_source = qqbar.forward_props[0].source_header.t_source;
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
    multiNd<Complex> mesprop;

    // Switch to Dirac-basis if desired.
    if (params.param.Dirac_basis)
    {
      qqbar.Dirac_basis = true;

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
    mescomp(mesprop,
	    quark_propagator[0],
	    quark_propagator[1],
	    phases, t0);

    // Convert the data into a mult1d
    multi1d<Complex> mesprop_1d;
    convertMescomp(mesprop_1d, mesprop, j_decay);

    // Save the qqbar output
    // ONLY SciDAC output format is supported!
    {
      XMLBufferWriter file_xml;
      push(file_xml, "qqbar");
      int id = 0;    // NEED TO FIX THIS - SOMETHING NON-TRIVIAL NEEDED
      write(file_xml, "id", id);
      pop(file_xml);

      XMLBufferWriter record_xml;
      push(record_xml, "QQbar");
      write(record_xml, ".", qqbar);  // do not write the outer group
      write(record_xml, "Config_info", gauge_xml);
      pop(record_xml);  // QQbar

      // Write the scalar data
      QDPFileWriter to(file_xml, params.prop.qqbar_file, 
		       QDPIO_SINGLEFILE, QDPIO_SERIAL, QDPIO_OPEN);
      write(to,record_xml,mesprop_1d);
      close(to);
    }

    pop(xml_out);    // qqbar

    QDPIO::cout << "QQbar finished" << endl;

    END_CODE();
  } 

};
