// $Id: spectrumQll_w.cc,v 1.1 2004-09-08 14:00:14 kostas Exp $
//
//! \file
//  \brief Computing the spectum of heavy light baryons
//

#include <iostream>
#include <cstdio>

#include "chroma.h"


using namespace QDP;


/*
 * Input 
 */
// Parameters which must be determined from the XML input
// and written to the XML output
struct Param_t
{
  bool Pt_snk;             // point sink
  bool Sl_snk;             // shell sink
  bool Wl_snk;             // wall sink

  WvfKind       wvf_kind;  // Wave function kind: gauge invariant
  multi1d<Real> wvf_param; // Array of width's or other parameters
  //   for "shell" source/sink wave function
  multi1d<int> wvfIntPar;  // Array of iter numbers to approx. Gaussian or
  //   terminate CG inversion for Wuppertal smearing

  multi1d<int> nrow;

  multi1d<int> Qsrc_coord ;
};


//! Propagators
struct Prop_t
{
  multi1d<string> prop_files;  // The files are expected to be in SciDAC format!
};


//! Mega-structure of parameters
struct Spectrum_input_t
{
  Param_t          param;
  Cfg_t            cfg;
  Prop_t           prop;
};


//! Propagator parameters
void read(XMLReader& xml, const string& path, Prop_t& input)
{
  XMLReader inputtop(xml, path);

  read(inputtop, "prop_files", input.prop_files);
}


//! Reader for parameters
void read(XMLReader& xml, const string& path, Param_t& param)
{
  XMLReader paramtop(xml, path);

  int version;
  read(paramtop, "version", version);
  
  read(paramtop, "Pt_snk", param.Pt_snk);
  read(paramtop, "Sl_snk", param.Sl_snk);

  read(paramtop, "wvf_kind", param.wvf_kind);
  read(paramtop, "wvf_param", param.wvf_param);
  read(paramtop, "wvfIntPar", param.wvfIntPar);

  if (param.wvf_param.size() != param.wvfIntPar.size())
  {
    QDPIO::cerr << "wvf_param size inconsistent with wvfintpar size" << endl;
    QDP_abort(1);
  }

  read(paramtop, "nrow", param.nrow);
  read(paramtop, "Qsrc_coord", param.Qsrc_coord);
}



// Reader for input parameters
void read(XMLReader& xml, const string& path, Spectrum_input_t& input)
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


//
// Main program
//
int main(int argc, char **argv)
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  START_CODE();

  // Input parameter structure
  Spectrum_input_t  input;

  // Instantiate xml reader for DATA
  XMLReader xml_in("DATA");

  // Read data
  read(xml_in, "/spectrumQll_w", input);

  // Specify lattice size, shape, etc.
  Layout::setLattSize(input.param.nrow);
  Layout::create();

  QDPIO::cout << " SPECTRUMQLL_W: Spectroscopy of heavy light baryons with Wilson fermions" << endl;

  /*
   * Sanity checks
   */
  if (input.param.wvf_param.size() != input.prop.prop_files.size())
  {
    QDPIO::cerr << "wvf_param size inconsistent with prop_files size" << endl;
    QDP_abort(1);
  }

  QDPIO::cout << endl << "     Gauge group: SU(" << Nc << ")" << endl;

  QDPIO::cout << "     volume: " << input.param.nrow[0];
  for (int i=1; i<Nd; ++i) {
    QDPIO::cout << " x " << input.param.nrow[i];
  }
  QDPIO::cout << endl;

  // Read in the configuration along with relevant information.
  multi1d<LatticeColorMatrix> u(Nd);
  XMLReader gauge_file_xml, gauge_xml;

  // Startup gauge
  gaugeStartup(gauge_file_xml, gauge_xml, u, input.cfg);

  // Check if the gauge field configuration is unitarized
  unitarityCheck(u);

  // Instantiate XML writer for XMLDAT
  XMLFileWriter xml_out("XMLDAT");
  push(xml_out, "spectrumQll_w");

  proginfo(xml_out);    // Print out basic program info

  // Write out the input
  write(xml_out, "Input", xml_in);

  // Write out the config info
  write(xml_out, "Config_info", gauge_xml);

  push(xml_out, "Output_version");
  write(xml_out, "out_version", 1);
  pop(xml_out);

  xml_out.flush();


  // First calculate some gauge invariant observables just for info.
  Double w_plaq, s_plaq, t_plaq, link;
  multi1d<DComplex> pollp(Nd);

  MesPlq(u, w_plaq, s_plaq, t_plaq, link);
  for(int mu = 0; mu < Nd; ++mu)
    polylp(u, pollp[mu], mu);

  push(xml_out, "Observables");
  write(xml_out, "w_plaq", w_plaq);
  write(xml_out, "s_plaq", s_plaq);
  write(xml_out, "t_plaq", t_plaq);
  write(xml_out, "link", link);
  write(xml_out, "pollp", pollp);
  pop(xml_out);

  xml_out.flush();

  // Keep an array of all the xml output buffers
  XMLArrayWriter xml_array(xml_out,input.prop.prop_files.size());
  push(xml_array, "Wilson_hadron_measurements");


  // Now loop over the various fermion masses
  for (int loop=0; loop < input.prop.prop_files.size(); ++loop)
  {
    // Read the quark propagator and extract headers
    LatticePropagator quark_propagator;
    ChromaProp_t prop_header;
    PropSource_t source_header;
    {
      XMLReader prop_file_xml, prop_record_xml;
      readQprop(prop_file_xml, 
		prop_record_xml, quark_propagator,
		input.prop.prop_files[loop], QDPIO_SERIAL);

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
    int j_decay = source_header.j_decay;
    Real Mass   = prop_header.FermActHandle->getMass();
    multi1d<int> boundary = prop_header.boundary;
    multi1d<int> t_source = source_header.t_source;

   
    //  phases with NO momenta
    SftMom phases(0, true, j_decay);

    // Next array element - name auto-written
    push(xml_array);
    write(xml_array, "loop", loop);
    write(xml_array, "Mass_mes", Mass);
    write(xml_array, "t_source", t_source);

    // Save prop input
    write(xml_array, "ForwardProp", prop_header);
    write(xml_array, "PropSource", source_header);

    // Sanity check - write out the norm2 of the forward prop in the j_decay direction
    // Use this for any possible verification
    {
      multi1d<Double> forward_prop_corr = sumMulti(localNorm2(quark_propagator), 
						   phases.getSet());

      push(xml_array, "Forward_prop_correlator");
      write(xml_array, "forward_prop_corr", forward_prop_corr);
      pop(xml_array);
    }

    // Determine what kind of source to use
    bool Pt_src = false;
    bool Sl_src = false;
    bool Wl_src = false;

    switch (source_header.source_type)
    {
    case SRC_TYPE_POINT_SOURCE:
      Pt_src = true;
      break;

    case SRC_TYPE_SHELL_SOURCE:
      Sl_src = true;
      break;

    case SRC_TYPE_WALL_SOURCE:
      Wl_src = true;
      break;

    default:
      QDPIO::cerr << "Unsupported source type" << endl;
      QDP_abort(1);
    }

    
    if (input.param.Pt_snk) 
      {
	if (Pt_src)
	  Qll(u,quark_propagator,Qsrc_coord, phases,xml_array, "Point_Point_Wilson_QllBaryons");
	if (Sl_src)
	  Qll(u,quark_propagator,Qsrc_coord, phases,xml_array, "Shell_Point_Wilson_QllBaryons");
	if (Wl_src)
	  Qll(u,quark_propagator,Qsrc_coord, phases,xml_array, "Wall_Point_Wilson_QllBaryons");
      } // end if (Pt_snk)
    
    // Convolute the quark propagator with the sink smearing function.
    // Make a copy of the quark propagator and then overwrite it with
    // the convolution. 
    if (input.param.Sl_snk) 
      {
	LatticePropagator quark_prop_smr;
	quark_prop_smr = quark_propagator;
	sink_smear2(u, quark_prop_smr, 
		    input.param.wvf_kind, 
		    input.param.wvf_param[loop],
		    input.param.wvfIntPar[loop], 
		    j_decay);
	if (Pt_src)
	  Qll(u,quark_prop_smr,Qsrc_coord, phases,xml_array, "Point_Shell_Wilson_QllBaryons");
	if (Sl_src)
	  Qll(u,quark_prop_smr,Qsrc_coord, phases,xml_array, "Shell_Shell_Wilson_QllBaryons");
	if (Wl_src)
	  Qll(u,quark_prop_smr,Qsrc_coord, phases,xml_array, "Wall_Shell_Wilson_QllBaryons");
      } // end if (Sl_snk)
    
    pop(xml_array);  // array element
    
  } // end for(loop)
  
  pop(xml_array);  // Wilson_spectroscopy
  pop(xml_out);  // spectrum_w
  
  xml_out.close();
  xml_in.close();
  
  END_CODE();
  
  // Time to bolt
  QDP_finalize();

  exit(0);
}
