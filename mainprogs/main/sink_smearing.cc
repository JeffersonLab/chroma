// $Id: sink_smearing.cc,v 1.6 2004-02-26 05:26:51 edwards Exp $
/*! \file
 * \brief Main program for sink-smearing quark propagators
 *
 * This test program reads in: quark propagator(s);
 *                   does    : sink smearing on quark propagator(s);
 *                   returns : quark propagator(s).
 *
 * Quark propagator read in by this program will not be changed.
 *
 * sink quark smearing: arbitrary powers of laplacian;
 *                      gaussian smearing;
 *                      displacement.
 *       
 * gauge link: either "bare" link or ape-smeared link.
 *       
 */

#include <iostream>
#include <cstdio>

#include "chroma.h"

using namespace QDP;

//! Propagators
struct Prop_t
{
  string        prop_file;

  QDP_volfmt_t  smeared_prop_volfmt; // Volume format is  SINGLEFILE or MULTIFILE
  string        smeared_prop_file;
};

//! Mega-structure of parameters
struct SinkSmear_input_t
{
  PropSink_t       param;
  Cfg_t            cfg;
  Prop_t           prop;
};


//! Propagator parameters
void read(XMLReader& xml, const string& path, Prop_t& input)
{
  XMLReader inputtop(xml, path);

  read(inputtop, "prop_file", input.prop_file);
  read(inputtop, "smeared_prop_file", input.smeared_prop_file);
  read(inputtop, "smeared_prop_volfmt", input.smeared_prop_volfmt);
}



// Reader for input parameters
void read(XMLReader& xml, const string& path, SinkSmear_input_t& input)
{
  XMLReader inputtop(xml, path);

  // Read all the input groups
  try
  {
    // Read program parameters
    read(inputtop, "Param", input.param);

    // Read in the gauge configuration info
    read(inputtop, "Cfg", input.cfg);

    // Read in the propagator(s) info
    read(inputtop, "Prop", input.prop);
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading prop data: " << e << endl;
    throw;
  }
}



int main(int argc, char **argv)
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  // Input parameter structure
  SinkSmear_input_t  input;

  // Instantiate xml reader for DATA
  XMLReader xml_in("DATA");

  // Read data
  read(xml_in, "/sink_smearing", input);

  // Specify lattice size, shape, etc.
  Layout::setLattSize(input.param.nrow);
  Layout::create();

  QDPIO::cout << " SINK_SMEAR: Sink smearing for propagators" << endl;

  // Read in the configuration along with relevant information.
  multi1d<LatticeColorMatrix> u(Nd);
  XMLReader gauge_xml;

  switch (input.cfg.cfg_type) 
  {
  case CFG_TYPE_SZIN :
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
  push(xml_out,"sink_smear");

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
  pop(xml_out);

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
	      input.prop.prop_file, QDPIO_SERIAL);
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
      throw;
    }
  }

  // Derived from input prop
  int  j_decay = source_header.j_decay;
  multi1d<int> boundary = prop_header.boundary;
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

  if (input.param.link_smear_num > 0)
  {
    int BlkMax = 100;	// Maximum number of blocking/smearing iterations
    Real BlkAccu = 1.0-5;	// Blocking/smearing accuracy

    for(int i=0; i < input.param.link_smear_num; ++i)
    {
      multi1d<LatticeColorMatrix> u_tmp(Nd);

      for(int mu = 0; mu < Nd; ++mu)
	if ( mu != j_decay )
	  APE_Smear(u_smr, u_tmp[mu], mu, 0,
		    input.param.link_smear_fact, BlkAccu, BlkMax, 
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
  if(input.param.sink_type == SNK_TYPE_SHELL_SINK)
  {
    // There should be a call to maksrc2 or some-such for general source smearing
    gausSmear(u_smr, quark_propagator, 
	      input.param.sinkSmearParam.wvf_param, 
	      input.param.sinkSmearParam.wvfIntPar, 
	      j_decay);
    laplacian(u_smr, quark_propagator, 
	      j_decay, input.param.laplace_power);
    displacement(u_smr, quark_propagator,
		 input.param.disp_length, input.param.disp_dir);
  }
	
  // Apply any sink constructions
  switch(input.param.wave_state)
  {
  case WAVE_TYPE_S_WAVE:
    break;
  case WAVE_TYPE_P_WAVE:
  {
    LatticePropagator temp;
    temp = quark_propagator;
    D_j(u_smr, temp, quark_propagator, input.param.direction);
  }
  break;
  case WAVE_TYPE_D_WAVE:   
  {
    LatticePropagator temp;
    temp = quark_propagator;
    DjDk(u_smr, temp, quark_propagator, input.param.direction);
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
  write(xml_out, "PropSink", input.param);

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

  xml_out.flush();

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
    write(record_xml, "PropSink", input.param);
    write(record_xml, "ForwardProp", prop_header);
    write(record_xml, "PropSource", source_header);
    write(record_xml, "Config_info", gauge_xml);
    pop(record_xml);

    // Write the source
    writeQprop(file_xml, record_xml, quark_propagator,
	       input.prop.smeared_prop_file, input.prop.smeared_prop_volfmt, QDPIO_SERIAL);
  }

  pop(xml_out);  // sink_smear

  xml_out.close();
  xml_in.close();

  // Time to bolt
  QDP_finalize();

  exit(0);
}

