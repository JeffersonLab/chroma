// $Id: make_source.cc,v 1.37 2005-02-28 03:34:46 edwards Exp $
/*! \file
 *  \brief Main code for source generation
 */

#include <iostream>
#include <cstdio>

#include "chroma.h"

using namespace Chroma;


struct Prop_t
{
  string          source_file;
  QDP_volfmt_t    source_volfmt;
};

struct Propagator_input_t
{
  PropSource_t    param;
  Cfg_t           cfg;
  Prop_t          prop;
};


//! Propagator input
void read(XMLReader& xml, const string& path, Prop_t& input)
{
  XMLReader inputtop(xml, path);

  read(inputtop, "source_file", input.source_file);
  read(inputtop, "source_volfmt", input.source_volfmt);
}


// Reader for input parameters
void read(XMLReader& xml, const string& path, Propagator_input_t& input)
{
  XMLReader inputtop(xml, path);

  // Read the input
  try
  {
    // Parameters for source construction
    read(inputtop, "Param", input.param);

    // Read in the gauge configuration info
    read(inputtop, "Cfg", input.cfg);

    // Read in the output propagator/source configuration info
    read(inputtop, "Prop", input.prop);
  }
  catch(const string& e)
  {
    QDP_error_exit("Error reading in make_source: %s", e.c_str());
  }
}




//! Propagator source generation
/*! \defgroup make_source Propagator source generation
 *  \ingroup main
 *
 * Main program for propagator source generation.
 */

int main(int argc, char **argv)
{
  // Put the machine into a known state
  //  QDP_initialize(&argc, &argv);

  ChromaInitialize(&argc, &argv);

  START_CODE();

  // Input parameter structure
  Propagator_input_t  input;

  // Instantiate xml reader for DATA
  XMLReader xml_in("./DATA");

  // Read data
  read(xml_in, "/make_source", input);

  Layout::setLattSize(input.param.nrow);
  Layout::create();

  QDPIO::cout << "MAKE_SOURCE: propagator source constructor" << endl;

  switch(input.param.source_type)
  {
  case SRC_TYPE_POINT_SOURCE:
    QDPIO::cout << "Point source" << endl;
    break;
  case SRC_TYPE_WALL_SOURCE:
    QDPIO::cout << "Wall source" << endl;
    break;
  case SRC_TYPE_SHELL_SOURCE:
    if (input.param.sourceSmearParam.wvf_kind != WVF_KIND_GAUGE_INV_GAUSSIAN)
    {
      QDPIO::cout << "Unsupported source smearing type" << endl;
      QDP_abort(1);
    }
    QDPIO::cout << "Smeared source wvf_param= " 
		<< input.param.sourceSmearParam.wvf_param 
		<< ": wvfIntPar= "
		<< input.param.sourceSmearParam.wvfIntPar << endl
		<< "Power of Laplacian operator= " << input.param.laplace_power << endl
		<< "Displacement length= " << input.param.disp_length
		<<": Displacement direction= " << input.param.disp_dir << endl;
    break;
  default:
    QDPIO::cout << "Unknown source_type" << endl;
    QDP_abort(1);
  }

  // Start up the gauge field
  multi1d<LatticeColorMatrix> u(Nd);
  XMLReader gauge_file_xml, gauge_xml;

  QDPIO::cout << "Initialize Gauge field" << endl;
  gaugeStartup(gauge_file_xml, gauge_xml, u, input.cfg);
  QDPIO::cout << "Gauge field initialized!" << endl;

  // Output
  XMLFileWriter&  xml_out = TheXMLOutputWriter::Instance();
  push(xml_out,"make_source");

  // Write out the input
  write(xml_out, "Input", xml_in);

  // Write out the config header
  write(xml_out, "Config_info", gauge_xml);

  xml_out.flush();


  // Check if the gauge field configuration is unitarized
  unitarityCheck(u);

  // Calculate some gauge invariant observables just for info.
  MesPlq(xml_out, "Observables", u);
  xml_out.flush();

  // Smear the gauge field if needed
  multi1d<LatticeColorMatrix> u_smr(Nd);
  u_smr = u;

  if (input.param.source_type == SRC_TYPE_SHELL_SOURCE && 
      input.param.link_smear_num > 0)
  {
    QDPIO::cout << "Smear gauge field" << endl;

    int BlkMax = 100;	// Maximum number of blocking/smearing iterations
    Real BlkAccu = 1.0e-5;	// Blocking/smearing accuracy

    for(int i=0; i < input.param.link_smear_num; ++i)
    {
      multi1d<LatticeColorMatrix> u_tmp(Nd);

      for(int mu = 0; mu < Nd; ++mu)
	if ( mu != input.param.j_decay )
	  APE_Smear(u_smr, u_tmp[mu], mu, 0, 
		    input.param.link_smear_fact, BlkAccu, BlkMax, 
		    input.param.j_decay);
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

  switch (input.param.source_type)
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
      srcfil(src_color_vec, input.param.t_source, color_source);

      // Smear the colour source if specified
      if(input.param.source_type == SRC_TYPE_SHELL_SOURCE)
      {
	// There should be a call to maksrc2 or some-such for general source smearing

	// displace the point source first, then smear
	// displacement has to be taken along negative direction.
	displacement(u_smr,src_color_vec,
		     (-1)*input.param.disp_length, input.param.disp_dir);

        if(input.param.wave_state == WAVE_TYPE_P_WAVE)
	  p_src(u_smr, src_color_vec, input.param.direction);

        if(input.param.wave_state == WAVE_TYPE_D_WAVE)   /* added */
	  d_src(u_smr, src_color_vec, input.param.direction);

        gausSmear(u_smr, src_color_vec, 
		  input.param.sourceSmearParam.wvf_param, 
		  input.param.sourceSmearParam.wvfIntPar, 
		  input.param.j_decay);

        laplacian(u_smr, src_color_vec, 
		  input.param.j_decay, input.param.laplace_power);

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
	       input.param.t_source[input.param.j_decay], 
	       input.param.j_decay, 
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
  

  // Sanity check - write out the norm2 of the source in the Nd-1 direction.
  // NOTE: this is NOT necessarily the time direction.
  // Use this for any possible verification.
  {
    // Initialize the slow Fourier transform phases
    SftMom phases(0, true, Nd-1);

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
    write(record_xml, "PropSource", input.param);
    write(record_xml, "Config_info", gauge_xml);
    pop(record_xml);
    
//    // Save a copy in the output
//    push(xml_out, "Source_info");
//    xml_out << file_xml;
//    xml_out << record_xml;
//    pop(xml_out);

    // Write the source
    writeQprop(file_xml, record_xml, quark_source,
	       input.prop.source_file, input.prop.source_volfmt, QDPIO_SERIAL);
  }


  pop(xml_out);  // make_source

  END_CODE();

  // Time to bolt
  ChromaFinalize();

  exit(0);
}
