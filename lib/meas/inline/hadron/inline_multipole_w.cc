// $Id: inline_multipole_w.cc,v 1.3 2005-04-19 17:11:07 edwards Exp $
/*! \file
 *  \brief Inline multipole measurement
 */

#include "meas/inline/hadron/inline_multipole_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "meas/hadron/multipole_w.h"
#include "io/qprop_io.h"
#include "meas/inline/make_xml_file.h"

namespace Chroma 
{ 
  namespace InlineMultipoleEnv 
  { 
    AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					    const std::string& path) 
    {
      return new InlineMultipole(InlineMultipoleParams(xml_in, path));
    }

    const std::string name = "MULTIPOLE";
    const bool registered = TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
  };


  // Reader for input parameters
  void read(XMLReader& xml, const string& path, InlineMultipoleParams::Param_t& param)
  {
    XMLReader inputtop(xml, path);

    int version;
    read(inputtop, "version", version);

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

    read(inputtop, "max_L", param.max_L);
  }


  // Writer for input parameters
  void write(XMLWriter& xml, const string& path, const InlineMultipoleParams::Param_t& param)
  {
    push(xml, path);

    int version = 1;
    write(xml, "version", version);
    write(xml, "max_L", param.max_L);

    pop(xml);
  }


  //! Propagator input
  void read(XMLReader& xml, const string& path, InlineMultipoleParams::Prop_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "seqprop_file", input.seqprop_file);
    read(inputtop, "GammaInsertion", input.GammaInsertion);
  }

  //! Propagator output
  void write(XMLWriter& xml, const string& path, const InlineMultipoleParams::Prop_t& input)
  {
    push(xml, path);

    write(xml, "seqprop_file", input.seqprop_file);
    write(xml, "GammaInsertion", input.GammaInsertion);

    pop(xml);
  }

  //! Multipole parameters
  void read(XMLReader& xml, const string& path, InlineMultipoleParams::Multipole_out_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "prop_file", input.prop_file);
    read(inputtop, "seqprops", input.seqprops);
  }

  //! Multipole parameters
  void write(XMLWriter& xml, const string& path, const InlineMultipoleParams::Multipole_out_t& input)
  {
    push(xml, path);

    write(xml, "prop_file", input.prop_file);
    write(xml, "seqprops", input.seqprops);

    pop(xml);
  }



  // Param stuff
  InlineMultipoleParams::InlineMultipoleParams() { frequency = 0; }

  InlineMultipoleParams::InlineMultipoleParams(XMLReader& xml_in, const std::string& path) 
  {
    try 
    {
      XMLReader inputtop(xml_in, path);

      read(inputtop, "Frequency", frequency);

      // Read program parameters
      read(inputtop, "Param", param);

      // Read in the multipole setup
      read(inputtop, "Multipole", pole);

      // Possible alternate XML file pattern
      if (inputtop.count("xml_file") != 0) 
      {
	read(inputtop, "xml_file", xml_file);
      }
    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << "Caught Exception reading XML: " << e << endl;
      QDP_abort(1);
    }
  }

  //! Write params
  void
  InlineMultipoleParams::write(XMLWriter& xml, const std::string& path) 
  {
    push(xml, path);
    
    Chroma::write(xml, "Param", param);
    Chroma::write(xml, "Multipole", pole);
    QDP::write(xml, "xml_file", xml_file);

    pop(xml);
  }


  // Function call
  void 
  InlineMultipole::operator()(const multi1d<LatticeColorMatrix>& u,
			      XMLBufferWriter& gauge_xml,
			      unsigned long update_no,
			      XMLWriter& xml_out) 
  {
    // If xml file not empty, then use alternate
    if (params.xml_file != "")
    {
      Handle<XMLFileWriter> xml(makeXMLFileWriter(params.xml_file, update_no));
      func(u, gauge_xml, update_no, *xml);
    }
    else
    {
      func(u, gauge_xml, update_no, xml_out);
    }
  }


  // Real work done here
  void 
  InlineMultipole::func(const multi1d<LatticeColorMatrix>& u,
			XMLBufferWriter& gauge_xml,
			unsigned long update_no,
			XMLWriter& xml_out) 
  {
    push(xml_out, "multipole");
    write(xml_out, "update_no", update_no);

    QDPIO::cout << " MULTIPOLE: Multipole measurements for Wilson-like fermions" << endl;

    proginfo(xml_out);    // Print out basic program info

    // Write out the input
    params.write(xml_out, "Input");

    // Write out the config info
    write(xml_out, "Config_info", gauge_xml);

    push(xml_out, "Output_version");
    write(xml_out, "out_version", 1);
    pop(xml_out);

    // First calculate some gauge invariant observables just for info.
    // This is really cheap.
    MesPlq(xml_out, "Observables", u);

    //
    // Read the quark propagator and extract headers
    //
    XMLReader prop_file_xml, prop_record_xml;
    LatticePropagator quark_propagator;
    ChromaProp_t prop_header;
    PropSource_t source_header;
    {
      QDPIO::cout << "Attempt to read forward propagator" << endl;
      readQprop(prop_file_xml, 
		prop_record_xml, quark_propagator,
		params.pole.prop_file, QDPIO_SERIAL);
   
      // Try to invert this record XML into a ChromaProp struct
      // Also pull out the id of this source
      try
      {
	read(prop_record_xml, "/Propagator/ForwardProp", prop_header);
	read(prop_record_xml, "/Propagator/PropSource", source_header);
      }
      catch (const string& e) 
      {
	QDPIO::cerr << "Error extracting forward_prop header: " << e << endl;
	QDP_abort(1);
      }

      // Save propagator input
      write(xml_out, "Propagator_file_info", prop_file_xml);
      write(xml_out, "Propagator_record_info", prop_record_xml);
    }
    QDPIO::cout << "Forward propagator successfully read" << endl;

    // Derived from input prop
    int  j_decay = source_header.j_decay;
    multi1d<int> t_source = source_header.t_source;

    // Sanity check - write out the norm2 of the forward prop in the j_decay direction
    // Use this for any possible verification
    {
      // Initialize the slow Fourier transform phases
      SftMom phases(0, true, j_decay);

      multi1d<Double> forward_prop_corr = sumMulti(localNorm2(quark_propagator), 
						   phases.getSet());

      push(xml_out, "Forward_prop_correlator");
      write(xml_out, "forward_prop_corr", forward_prop_corr);
      pop(xml_out);
    }


    //
    // Big nested loop over seqprops
    //
    push(xml_out, "Multipole_measurements");

    XMLArrayWriter  xml_seq_src(xml_out, params.pole.seqprops.size());
    push(xml_seq_src, "Sequential_source");

    for (int seq_src_ctr = 0; seq_src_ctr < params.pole.seqprops.size(); ++seq_src_ctr) 
    {
      push(xml_seq_src);
      write(xml_seq_src, "seq_src_ctr", seq_src_ctr);

      // Read the sequential propagator
      // Read the quark propagator and extract headers
      LatticePropagator seq_quark_prop;
      SeqSource_t seqsource_header;
      {
	XMLReader seqprop_file_xml, seqprop_record_xml;
	readQprop(seqprop_file_xml, 
		  seqprop_record_xml, seq_quark_prop,
		  params.pole.seqprops[seq_src_ctr].seqprop_file, QDPIO_SERIAL);

	// Try to invert this record XML into a ChromaProp struct
	// Also pull out the id of this source
	// NEED SECURITY HERE - need a way to cross check props. Use the ID.
	try
	{
	  read(seqprop_record_xml, "/SequentialProp/SeqSource", seqsource_header);
	}
	catch (const string& e) 
	{
	  QDPIO::cerr << "Error extracting seqprop header: " << e << endl;
	  QDP_abort(1);
	}

	// Save seqprop input
	write(xml_seq_src, "SequentialProp_file_info", seqprop_file_xml);
	write(xml_seq_src, "SequentialProp_record_info", seqprop_record_xml);
      }
      QDPIO::cout << "Sequential propagator successfully read" << endl;

      // Sanity check - write out the norm2 of the forward prop in the j_decay direction
      // Use this for any possible verification
      {
	// Initialize the slow Fourier transform phases
	SftMom phases(0, true, Nd-1);
      
	multi1d<Double> backward_prop_corr = sumMulti(localNorm2(seq_quark_prop), 
						      phases.getSet());
      
	push(xml_seq_src, "Backward_prop_correlator");
	write(xml_seq_src, "backward_prop_corr", backward_prop_corr);
	pop(xml_seq_src);
      }

      // Derived from input seqprop
      string seq_src = seqsource_header.seq_src;
      QDPIO::cout << "Seqsource name = " << seqsource_header.seq_src << endl;
      int           t_sink   = seqsource_header.t_sink;
      multi1d<int>  sink_mom = seqsource_header.sink_mom;

      write(xml_seq_src, "hadron_type", "HADRON");
      write(xml_seq_src, "seq_src", seq_src);
      write(xml_seq_src, "t_source", t_source);
      write(xml_seq_src, "t_sink", t_sink);
      write(xml_seq_src, "sink_mom", sink_mom);
      write(xml_seq_src, "max_L", params.param.max_L);
      write(xml_seq_src, "GammaInsertion", params.pole.seqprops[seq_src_ctr].GammaInsertion);
	
      // Now the 3pt contractions
      multipole(quark_propagator, seq_quark_prop, params.pole.seqprops[seq_src_ctr].GammaInsertion, 
		params.param.max_L, j_decay, t_source[j_decay], xml_seq_src, "Multipole");

      pop(xml_seq_src);   // elem
    } // end loop over sequential sources

    pop(xml_seq_src);  // Sequential_source

    pop(xml_out);  // Multipole_measurements

    // Close the namelist output file XMLDAT
    pop(xml_out);     // multipole

    QDPIO::cout << "Multipole ran successfully" << endl;

    END_CODE();
  } 

};
