// $Id: inline_polar_source_w.cc,v 1.4 2007-02-25 22:39:28 edwards Exp $
// Code for sequential source construction

#error "THIS CODE IS NOT READY YET"

#include <iostream>
#include <cstdio>

#include "chroma.h"

using namespace Chroma;


/*
 * Input 
 */

//! Propagators
struct Prop_t
{
  string           prop_file;  // The file is expected to be in SciDAC format!
  string           seqsource_file;  // The file is expected to be in SciDAC format!
  QDP_volfmt_t     seqsource_volfmt;
};


//! Mega-structure of all input
struct SeqSource_input_t
{
  PropSource_t      param;
  Cfg_t             cfg;
  Prop_t            prop;
};


//! Propagator parameters
void read(XMLReader& xml, const string& path, Prop_t& input)
{
  XMLReader inputtop(xml, path);

  read(inputtop, "prop_file", input.prop_file);
  read(inputtop, "seqsource_file", input.seqsource_file);
  read(inputtop, "seqsource_volfmt", input.seqsource_volfmt);
}


// Reader for input parameters
void read(XMLReader& xml, const string& path, SeqSource_input_t& input)
{
  XMLReader inputtop(xml, path);

  // Read the input
  try
  {
    // The parameters holds the version number
    read(inputtop, "Param", input.param);

    // Read in the gauge configuration info
    read(inputtop, "Cfg", input.cfg);

    // Read in the forward_prop/seqsource info
    read(inputtop, "Prop", input.prop);
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading data: " << e << endl;
    throw;
  }
}


//! Sequential source generation
/*
 *  \defgroup seqsource Sequential source generation
 *  \ingroup main
 *
 *  Read quark propagators, convert into sequential sources
 *
 */

int main(int argc, char **argv)
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  START_CODE();

  // Input parameter structure
  SeqSource_input_t  input;

  // Instantiate xml reader for DATA
  XMLReader xml_in("DATA");

  // Read data
  read(xml_in, "/seqsource", input);

  // Specify lattice size, shape, etc.
  Layout::setLattSize(input.param.nrow);
  Layout::create();

  // Sanity checks
  QDPIO::cout << endl << "     Gauge group: SU(" << Nc << ")" << endl;

  QDPIO::cout << "     Computing sequential source of type "
	      << input.param.source_type << endl;
  
  QDPIO::cout << "     Volume: " << input.param.nrow[0];
  for (int i=1; i<Nd; ++i) {
    QDPIO::cout << " x " << input.param.nrow[i];
  }
  QDPIO::cout << endl;


  // Read in the configuration along with relevant information.
  multi1d<LatticeColorMatrix> u(Nd);
  XMLReader gauge_file_xml, gauge_xml;

  // Startup gauge
  gaugeStartup(gauge_file_xml, gauge_xml, u, input.cfg);

  // Instantiate XML writer for XMLDAT
  XMLFileWriter xml_out("XMLDAT");
  push(xml_out, "seqsource");

  proginfo(xml_out);    // Print out basic program info

  // Write out the input
  write(xml_out, "Input", xml_in);

  // Write out the config header
  write(xml_out, "Config_info", gauge_xml);

  push(xml_out, "Output_version");
  write(xml_out, "out_version", 1);
  pop(xml_out);

  xml_out.flush();


  // Check if the gauge field configuration is unitarized
  unitarityCheck(u);

  // Calculate some gauge invariant observables just for info.
  Double w_plaq, s_plaq, t_plaq, link;
  MesPlq(u, w_plaq, s_plaq, t_plaq, link);

  push(xml_out, "Observables");
  write(xml_out, "w_plaq", w_plaq);
  write(xml_out, "s_plaq", s_plaq);
  write(xml_out, "t_plaq", t_plaq);
  write(xml_out, "link", link);
  pop(xml_out);

  xml_out.flush();


  //
  // Read the quark propagator and extract headers
  //
  LatticePropagator quark_propagator;
  ChromaProp_t prop_header;
  PropSource_t source_header;
  {
    XMLReader prop_file_xml, prop_record_xml;

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

    // Save prop input
    write(xml_out, "Propagator_info", prop_record_xml);
  }

  // Derived from input prop
  int  j_decay = source_header.j_decay;
// multi1d<int> boundary = prop_header.boundary;   // not currently needed
  multi1d<int> t_source = source_header.t_source;

  // Initialize the slow Fourier transform phases
  SftMom phases(0, true, j_decay);

  // Sanity check - write out the norm2 of the forward prop in the j_decay direction
  // Use this for any possible verification
  {
    multi1d<Double> forward_prop_corr = sumMulti(localNorm2(quark_propagator), 
						 phases.getSet());

    push(xml_out, "Forward_prop_correlator");
    write(xml_out, "forward_prop_corr", forward_prop_corr);
    pop(xml_out);
  }

  //------------------ Start main body of calculations -----------------------------

  //
  // Construct the sequential source
  //
  LatticePropagator quark_prop_src = zero;
  LatticePropagator quark_prop_tmp;
  LatticeColorMatrix uu;
  multi1d<LatticeColorMatrix> ua(Nd);
  LatticeComplex a;
  Real alocr, aloci;
  Complex aloc;
  multi1d<int> size = input.param.nrow;
  multi1d<int> coords(Nd);

// Create a LatticeComplex em field

  for (int i=0; i<size[0]; ++i) {
  for (int j=0; j<size[1]; ++j) {
  for (int k=0; k<size[2]; ++k) {
  for (int l=0; l<size[3]; ++l) {

// If we want to do a magnetic insertion (note the pretty arbitrary,
// but not completely stupid, strength of the external gauge field -
// in final results, this strength is an overall factor and will be
// divided out in data analysis):

//    alocr=0.0;
//    aloci = (float) j;
//
//    if ( j > (size[1]/2) ) {
//      aloci = aloci - (float) size[1];
//    };
//
//    aloci=aloci*10.0*6.2831852/(size[1]*size[0]);


// If we want to do an electric insertion:

//    alocr=0.0;
//    aloci = (float) l;
//    aloci=aloci*10.0*6.2831852/(size[2]*size[3]);


//  If we want to do an A^2 insertion (electric):

    alocr = (float) (l*l);
    aloci = 0.0;
    alocr=alocr*10.0*6.2831852/(size[2]*size[3]);
    alocr=alocr*10.0*6.2831852/(size[2]*size[3]);
    alocr=alocr*(-0.5000000);


//  So let's write that into the field:

    aloc = cmplx(alocr,aloci);
    coords[0]=i;
    coords[1]=j;
    coords[2]=k;
    coords[3]=l;
    pokeSite(a, aloc, coords);
  };
  };
  };
  };

// ok, so now we have a LatticeComplex containing our em field
// next, make a LatticeColorMatrix out of it

  ua[0] = zero;
  ua[1] = zero;
  ua[2] = zero;
  ua[3] = zero;

// Pick which component of the gauge field we're making nontrivial
// (for electric insertion, ua[2], for magnetic insertion, ua[0]);

  for (int elem=0; elem<Nc; ++elem) {
//    pokeColor(ua[0], a, elem, elem);
//    pokeColor(ua[1], a, elem, elem);
    pokeColor(ua[2], a, elem, elem);
//    pokeColor(ua[3], a, elem, elem);
  };

// and now put it all together ...

  int i=1;

  for (int j=0; j<Nd; ++j) {

  uu = ua[j]*u[j];
  quark_prop_tmp = quark_prop_src;
  quark_prop_src = quark_prop_tmp - uu*shift(quark_propagator, FORWARD, j);
  quark_prop_tmp = quark_prop_src;
  quark_prop_src = quark_prop_tmp + uu*(Gamma(i)*shift(quark_propagator, FORWARD, j));
  quark_prop_tmp = quark_prop_src;
  quark_prop_src = quark_prop_tmp - adj(shift(uu, BACKWARD, j))*shift(quark_propagator, BACKWARD, j);
  quark_prop_tmp = quark_prop_src;
  quark_prop_src = quark_prop_tmp - adj(shift(uu, BACKWARD, j))*(Gamma(i)*shift(quark_propagator, BACKWARD, j));
  i=i*2;

  };

  Double factor=0.5000000;
  quark_prop_tmp = quark_prop_src;

//  If we want to apply the Wilson Dirac operator (in the absence of
//  external fields - this is just for testing):
//  Double kappa=0.11000000;
//  quark_prop_src = (factor/kappa)*quark_propagator+factor*quark_prop_tmp;

//  If we want to insert a coupling to the external em field:
  quark_prop_src = factor*quark_prop_tmp;

//  If we just want to apply the identity operation (just for test
//  purposes):
//  quark_prop_src = quark_propagator;


  // Sanity check - write out the norm2 of the propagator source in the j_decay direction
  // Use this for any possible verification
  {
    multi1d<Double> seqsource_corr = sumMulti(localNorm2(quark_prop_src), 
					      phases.getSet());
	
    push(xml_out, "SeqSource_correlator");
    write(xml_out, "seqsource_corr", seqsource_corr);
    pop(xml_out);
  }


  /*
   *  Write the sequential source out to disk
   */
  {
    XMLBufferWriter file_xml;
    push(file_xml, "make_source");
    write(file_xml, "id", uniqueId());  // NOTE: new ID form
    pop(file_xml);

    XMLBufferWriter record_xml;
    push(record_xml, "MakeSource");
    write(record_xml, "PropSource", input.param);
    write(record_xml, "Config_info", gauge_xml);
    pop(record_xml);  // SequentialSource

    // Write the seqsource
    writeQprop(file_xml, record_xml, quark_prop_src,
	       input.prop.seqsource_file, 
	       input.prop.seqsource_volfmt, QDPIO_SERIAL);

    QDPIO::cout << "Sequential source successfully written" << endl;
  }

  pop(xml_out);    // seqsource

  xml_out.close();
  xml_in.close();

  // Time to bolt
  QDP_finalize();

  END_CODE();

  exit(0);
}

