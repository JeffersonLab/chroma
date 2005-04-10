// $Id: inline_make_source_w.cc,v 1.2 2005-04-10 17:06:22 edwards Exp $
/*! \file
 * \brief Inline construction of make_source
 *
 * Construct source for propagator calculations
 */

#include "meas/inline/hadron/inline_make_source_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "meas/smear/ape_smear.h"
#include "meas/smear/gaus_smear.h"
#include "meas/smear/displacement.h"
#include "meas/smear/laplacian.h"
#include "meas/hadron/srcfil.h"
#include "meas/hadron/walfil_w.h"
#include "meas/sources/p_src_w.h"
#include "meas/sources/d_src_w.h"
#include "util/ferm/transf.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"

namespace Chroma 
{ 
  namespace InlineMakeSourceEnv 
  { 
    AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					    const std::string& path) 
    {
      return new InlineMakeSource(InlineMakeSourceParams(xml_in, path));
    }

    const std::string name = "MAKE_SOURCE";
    const bool registered = TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
  };


  //! Propagator input
  void read(XMLReader& xml, const string& path, InlineMakeSourceParams::Prop_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "source_file", input.source_file);
    read(inputtop, "source_volfmt", input.source_volfmt);
  }

  //! Propagator output
  void write(XMLWriter& xml, const string& path, const InlineMakeSourceParams::Prop_t& input)
  {
    push(xml, path);

    write(xml, "source_file", input.source_file);
    write(xml, "source_volfmt", input.source_volfmt);

    pop(xml);
  }


  // Param stuff
  InlineMakeSourceParams::InlineMakeSourceParams() { frequency = 0; }

  InlineMakeSourceParams::InlineMakeSourceParams(XMLReader& xml_in, const std::string& path) 
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
  InlineMakeSourceParams::write(XMLWriter& xml_out, const std::string& path) 
  {
    push(xml_out, path);
    
    // Parameters for source construction
    Chroma::write(xml_out, "Param", param);

    // Write out the output propagator/source configuration info
    Chroma::write(xml_out, "Prop", prop);

    pop(xml_out);
  }


  void 
  InlineMakeSource::operator()(const multi1d<LatticeColorMatrix>& u,
			       XMLBufferWriter& gauge_xml,
			       unsigned long update_no,
			       XMLWriter& xml_out) 
  {
    START_CODE();

    push(xml_out, "make_source");
    write(xml_out, "update_no", update_no);

    QDPIO::cout << "MAKE_SOURCE: propagator source constructor" << endl;

    switch(params.param.source_type)
    {
    case SRC_TYPE_POINT_SOURCE:
      QDPIO::cout << "Point source" << endl;
      break;
    case SRC_TYPE_WALL_SOURCE:
      QDPIO::cout << "Wall source" << endl;
      break;
    case SRC_TYPE_SHELL_SOURCE:
      if (params.param.sourceSmearParam.wvf_kind != WVF_KIND_GAUGE_INV_GAUSSIAN)
      {
	QDPIO::cout << "Unsupported source smearing type" << endl;
	QDP_abort(1);
      }
      QDPIO::cout << "Smeared source wvf_param= " 
		  << params.param.sourceSmearParam.wvf_param 
		  << ": wvfIntPar= "
		  << params.param.sourceSmearParam.wvfIntPar << endl
		  << "Power of Laplacian operator= " << params.param.laplace_power << endl
		  << "Displacement length= " << params.param.disp_length
		  <<": Displacement direction= " << params.param.disp_dir << endl;
      break;
    default:
      QDPIO::cout << "Unknown source_type" << endl;
      QDP_abort(1);
    }

    proginfo(xml_out);    // Print out basic program info

    // Write out the input
    params.write(xml_out, "Input");

    // Write out the config header
    write(xml_out, "Config_info", gauge_xml);

    // Calculate some gauge invariant observables just for info.
    MesPlq(xml_out, "Observables", u);

    // Smear the gauge field if needed
    multi1d<LatticeColorMatrix> u_smr(Nd);
    u_smr = u;

    if (params.param.source_type == SRC_TYPE_SHELL_SOURCE && 
	params.param.link_smear_num > 0)
    {
      QDPIO::cout << "Smear gauge field" << endl;

      int BlkMax = 100;	// Maximum number of blocking/smearing iterations
      Real BlkAccu = 1.0e-5;	// Blocking/smearing accuracy

      for(int i=0; i < params.param.link_smear_num; ++i)
      {
	multi1d<LatticeColorMatrix> u_tmp(Nd);

	for(int mu = 0; mu < Nd; ++mu)
	  if ( mu != params.param.j_decay )
	    APE_Smear(u_smr, u_tmp[mu], mu, 0, 
		      params.param.link_smear_fact, BlkAccu, BlkMax, 
		      params.param.j_decay);
	  else
	    u_tmp[mu] = u_smr[mu];

	u_smr = u_tmp;
      }

      QDPIO::cout << "Gauge field APE-smeared!" << endl;
    }


    //
    // Loop over the source color and spin, creating the source
    // and calling the relevant propagator routines. The QDP
    // terminology is that a propagator is a matrix in color
    // and spin space
    //
    // For this calculation, a smeared source is used. A point
    // source is first constructed and then smeared. If a user
    // only wanted a point source, then remove the smearing stuff
    //
    LatticePropagator quark_source;

    switch (params.param.source_type)
    {
      //
      // Gauge inv. point or shell sources within some S_WAVE, P_WAVE, state etc.
      //
    case SRC_TYPE_POINT_SOURCE:
    case SRC_TYPE_SHELL_SOURCE:
    {
      for(int color_source = 0; color_source < Nc; ++color_source)
      {
	QDPIO::cout << "color = " << color_source << endl; 

	LatticeColorVector src_color_vec = zero;

	// Make a point source at coordinates t_source
	srcfil(src_color_vec, params.param.t_source, color_source);

	// Smear the colour source if specified
	if(params.param.source_type == SRC_TYPE_SHELL_SOURCE)
	{
	  // There should be a call to maksrc2 or some-such for general source smearing

	  // displace the point source first, then smear
	  // displacement has to be taken along negative direction.
	  displacement(u_smr,src_color_vec,
		       (-1)*params.param.disp_length, params.param.disp_dir);

	  if(params.param.wave_state == WAVE_TYPE_P_WAVE)
	    p_src(u_smr, src_color_vec, params.param.direction);

	  if(params.param.wave_state == WAVE_TYPE_D_WAVE)   /* added */
	    d_src(u_smr, src_color_vec, params.param.direction);

	  gausSmear(u_smr, src_color_vec, 
		    params.param.sourceSmearParam.wvf_param, 
		    params.param.sourceSmearParam.wvfIntPar, 
		    params.param.j_decay);

	  laplacian(u_smr, src_color_vec, 
		    params.param.j_decay, params.param.laplace_power);

	  //power = 1 for one laplacian operator
	}

	for(int spin_source = 0; spin_source < Ns; ++spin_source)
	{
	  QDPIO::cout << "spin = " << spin_source << endl; 

	  // Insert a ColorVector into spin index spin_source
	  // This only overwrites sections, so need to initialize first
	  LatticeFermion chi = zero;

	  CvToFerm(src_color_vec, chi, spin_source);
      

	  /*
	   *  Move the source to the appropriate components
	   *  of quark source.
	   */
	  FermToProp(chi, quark_source, color_source, spin_source);
	}
      }
    }
    break;

    case SRC_TYPE_WALL_SOURCE:
    {
      for(int color_source = 0; color_source < Nc; ++color_source)
      {
	for(int spin_source = 0; spin_source < Ns; ++spin_source)
	{
	  // Wall fill a fermion source. Insert it into the propagator source
	  LatticeFermion chi;
	  walfil(chi, 
		 params.param.t_source[params.param.j_decay], 
		 params.param.j_decay, 
		 color_source, spin_source);
	  FermToProp(chi, quark_source, color_source, spin_source);
	}
      }
    }
    break;
  
    default:
      QDPIO::cout << "Unsupported source type" << endl;
      QDP_abort(1);
    }
  

    // Sanity check - write out the norm2 of the source in the j_decay direction.
    // Use this for any possible verification.
    {
      // Initialize the slow Fourier transform phases
      SftMom phases(0, true, params.param.j_decay);

      multi1d<Double> source_corr = sumMulti(localNorm2(quark_source),
					     phases.getSet());

      push(xml_out, "Source_correlator");
      write(xml_out, "source_corr", source_corr);
      pop(xml_out);
    }
 

    // Now write the source
    // ONLY SciDAC output format is supported!!!
    {
      XMLBufferWriter file_xml;
      push(file_xml, "make_source");
      int id = 0;    // NEED TO FIX THIS - SOMETHING NON-TRIVIAL NEEDED
      write(file_xml, "id", id);
      pop(file_xml);

      XMLBufferWriter record_xml;
      push(record_xml, "MakeSource");
      write(record_xml, "PropSource", params.param);
      write(record_xml, "Config_info", gauge_xml);
      pop(record_xml);
    
//    // Save a copy in the output
//    push(xml_out, "Source_info");
//    xml_out << file_xml;
//    xml_out << record_xml;
//    pop(xml_out);

      // Write the source
      writeQprop(file_xml, record_xml, quark_source,
		 params.prop.source_file, params.prop.source_volfmt, QDPIO_SERIAL);

      QDPIO::cout << "Source successfully written" << endl;
    }
    
    QDPIO::cout << "Make_source ran successfully" << endl;

    pop(xml_out);  // make_source

    END_CODE();
  } 

};
