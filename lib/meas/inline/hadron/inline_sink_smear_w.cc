// $Id: inline_sink_smear_w.cc,v 1.1 2005-04-06 04:34:54 edwards Exp $
/*! \file
 * \brief Inline construction of sink_smear
 *
 * Sink smear propagators
 */

#include "meas/inline/hadron/inline_sink_smear_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "meas/smear/ape_smear.h"
#include "meas/smear/gaus_smear.h"
#include "meas/smear/displacement.h"
#include "meas/smear/laplacian.h"
#include "meas/hadron/D_j_w.h"
#include "meas/hadron/DjDk_w.h"
#include "util/ferm/transf.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"

namespace Chroma 
{ 
  namespace InlineSinkSmearEnv 
  { 
    AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					    const std::string& path) 
    {
      return new InlineSinkSmear(InlineSinkSmearParams(xml_in, path));
    }

    const std::string name = "SINK_SMEAR";
    const bool registered = TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
  };


  //! Propagator input
  void read(XMLReader& xml, const string& path, InlineSinkSmearParams::Prop_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "prop_file", input.prop_file);
    read(inputtop, "smeared_prop_file", input.smeared_prop_file);
    read(inputtop, "smeared_prop_volfmt", input.smeared_prop_volfmt);
  }

  //! Propagator output
  void write(XMLWriter& xml, const string& path, const InlineSinkSmearParams::Prop_t& input)
  {
    push(xml, path);

    write(xml, "prop_file", input.prop_file);
    write(xml, "smeared_prop_file", input.smeared_prop_file);
    write(xml, "smeared_prop_volfmt", input.smeared_prop_volfmt);

    pop(xml);
  }


  // Param stuff
  InlineSinkSmearParams::InlineSinkSmearParams() { frequency = 0; }

  InlineSinkSmearParams::InlineSinkSmearParams(XMLReader& xml_in, const std::string& path) 
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
  InlineSinkSmearParams::write(XMLWriter& xml_out, const std::string& path) 
  {
    push(xml_out, path);
    
    // Parameters for source construction
    Chroma::write(xml_out, "Param", param);

    // Write out the output propagator/source configuration info
    Chroma::write(xml_out, "Prop", prop);

    pop(xml_out);
  }


  void 
  InlineSinkSmear::operator()(const multi1d<LatticeColorMatrix>& u,
			      XMLBufferWriter& gauge_xml,
			      unsigned long update_no,
			      XMLWriter& xml_out) 
  {
    push(xml_out, "sink_smear");
    write(xml_out, "update_no", update_no);

    QDPIO::cout << " SINK_SMEAR: Sink smearing for propagators" << endl;

    // Write out the input
    params.write(xml_out, "Input");

    // Write out the config header
    write(xml_out, "Config_info", gauge_xml);

    // Calculate some gauge invariant observables just for info.
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
		params.prop.prop_file, QDPIO_SERIAL);
      QDPIO::cout << "Forward propagator successfully read" << endl;
   
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
    }

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



    /*
     * Smear the gauge field if needed
     */
    multi1d<LatticeColorMatrix> u_smr(Nd);
    u_smr = u;

    if (params.param.sink_type == SNK_TYPE_SHELL_SINK &&
	params.param.link_smear_num > 0)
    {
      int BlkMax = 100;	// Maximum number of blocking/smearing iterations
      Real BlkAccu = 1.0e-5;	// Blocking/smearing accuracy

      for(int i=0; i < params.param.link_smear_num; ++i)
      {
	multi1d<LatticeColorMatrix> u_tmp(Nd);

	for(int mu = 0; mu < Nd; ++mu)
	  if ( mu != j_decay )
	    APE_Smear(u_smr, u_tmp[mu], mu, 0,
		      params.param.link_smear_fact, BlkAccu, BlkMax, 
		      j_decay);
	  else
	    u_tmp[mu] = u_smr[mu];
      
	u_smr = u_tmp;
      }
      QDPIO::cout << "Gauge field APE-smeared!" << endl;
    }


    /*
     * Now apply appropriate sink smearing
     */
    if(params.param.sink_type == SNK_TYPE_SHELL_SINK)
    {
      // There should be a call to maksrc2 or some-such for general source smearing
      gausSmear(u_smr, quark_propagator, 
		params.param.sinkSmearParam.wvf_param, 
		params.param.sinkSmearParam.wvfIntPar, 
		j_decay);
      laplacian(u_smr, quark_propagator, 
		j_decay, params.param.laplace_power);
      displacement(u_smr, quark_propagator,
		   params.param.disp_length, params.param.disp_dir);
    }
	
    // Apply any sink constructions
    switch(params.param.wave_state)
    {
    case WAVE_TYPE_S_WAVE:
      break;
    case WAVE_TYPE_P_WAVE:
    {
      LatticePropagator temp;
      temp = quark_propagator;
      D_j(u_smr, temp, quark_propagator, params.param.direction);
    }
    break;
    case WAVE_TYPE_D_WAVE:   
    {
      LatticePropagator temp;
      temp = quark_propagator;
      DjDk(u_smr, temp, quark_propagator, params.param.direction);
    }
    break;
    default: 
      QDPIO::cerr << "Invalid wave_state" << endl;
      QDP_abort(1);
      break;
    } 
  

    /*
     * Save sink smeared propagator
     */
    // Save some info
    write(xml_out, "PropSource", source_header);
    write(xml_out, "ForwardProp", prop_header);
    write(xml_out, "PropSink", params.param);

    // Sanity check - write out the propagator (pion) correlator in the Nd-1 direction
    {
      // Initialize the slow Fourier transform phases
      SftMom phases(0, true, j_decay);

      multi1d<Double> prop_corr = sumMulti(localNorm2(quark_propagator), 
					   phases.getSet());

      push(xml_out, "SinkSmearedProp_correlator");
      write(xml_out, "sink_smeared_prop_corr", prop_corr);
      pop(xml_out);
    }


    // Save the propagator
    // ONLY SciDAC output format is supported!
    {
      XMLBufferWriter file_xml;
      push(file_xml, "sink_smear");
      int id = 0;    // NEED TO FIX THIS - SOMETHING NON-TRIVIAL NEEDED
      write(file_xml, "id", id);
      pop(file_xml);

      XMLBufferWriter record_xml;
      push(record_xml, "SinkSmear");
      write(record_xml, "PropSink", params.param);
      write(record_xml, "ForwardProp", prop_header);
      write(record_xml, "PropSource", source_header);
      write(record_xml, "Config_info", gauge_xml);
      pop(record_xml);

      // Write the source
      writeQprop(file_xml, record_xml, quark_propagator,
		 params.prop.smeared_prop_file, 
		 params.prop.smeared_prop_volfmt, QDPIO_SERIAL);
    }

    pop(xml_out);  // sink_smear

    END_CODE();
  } 

};
