// $Id: collect_propcomp.cc,v 1.1 2004-04-20 13:08:12 bjoo Exp $
// $Log: collect_propcomp.cc,v $
// Revision 1.1  2004-04-20 13:08:12  bjoo
// Added multi mass component based solves and propagator collection program
//
// Revision 1.3  2004/04/16 22:03:59  bjoo
// Added Zolo 4D test files and tidyed up
//
// Revision 1.2  2004/04/16 20:18:03  bjoo
// Zolo seems to work now
//
// Revision 1.1  2004/04/16 17:04:49  bjoo
// Added multi_propagator for Zolo4D multi massery. Seems to work even
//
// Revision 1.49  2004/04/15 14:43:25  bjoo
// Added generalised qprop_io FermAct reading
//
// Revision 1.48  2004/04/06 04:20:33  edwards
// Added SZINQIO support.
//
// Revision 1.47  2004/04/01 18:10:22  edwards
// Added support for non-relativistic quark props.
//
// Revision 1.46  2004/02/23 03:13:58  edwards
// Major overhaul of input/output model! Now using EXCLUSIVELY
// SciDAC propagator format for propagators. Now, Param part of input
// files directly matches source/sink/propagator/seqprop headers
// of propagators. All ``known'' input of a propagator is derived
// from its header(s) and used for subsequent calculations.
//
// Revision 1.45  2004/02/11 12:51:35  bjoo
// Stripped out Read() and Write()
//
// Revision 1.44  2004/02/07 04:51:58  edwards
// Removed SSE hack - pushed it down into the SSE code where it belongs.
//
// Revision 1.43  2004/02/06 22:31:00  edwards
// Put in sse hack for the short term.
//
// Revision 1.42  2004/02/06 17:39:05  edwards
// Added a flush to xml_out.
//
// Revision 1.41  2004/02/05 04:18:56  edwards
// Changed call of quarkProp4 to write to xml_out instead of xml buffer.
//
// Revision 1.40  2004/02/04 17:01:55  edwards
// Changed getSubset() to the now correct getSet().
//
// Revision 1.39  2004/01/31 23:22:01  edwards
// Added proginfo call.
//
// Revision 1.38  2004/01/30 21:35:22  kostas
// added calls to calculate mres for dwf
// 
/*! \file
 *  \brief Main code for propagator generation
 */

#include <iostream>
#include <cstdio>
#include <iomanip>
using namespace std;
#include "chroma.h"

using namespace QDP;



/*
 * Input 
 */
struct Prop_t
{
  string          source_file;
  string          prop_file;
  QDP_volfmt_t    prop_volfmt;
};


struct Component_t { 
  int color;
  int spin;
};

struct PropagatorComponent_input_t
{
  ChromaMultiProp_t     param;
  Cfg_t            cfg;
  Prop_t           prop;
  multi1d<Component_t> components;
};


void read(XMLReader& xml, const string& path, Component_t &comp)
{
  XMLReader top(xml,path);
  try {
    read(top, "color", comp.color);
    read(top, "spin",  comp.spin);
  }
  catch( const string& e ) {
    QDPIO::cerr << "Caught Exception : " << e << endl;
    QDP_abort(1);
  }  
  if( comp.color < 0 || comp.color >= Nc ) { 
    QDPIO::cerr << "Component color >= Nc. color = " << comp.color << endl;
    QDP_abort(1);
  }

  if( comp.spin < 0 || comp.spin >= Ns ) { 
    QDPIO::cerr << "Component spin >= Ns.  spin = " << comp.spin << endl;
    QDP_abort(1);
  }
}

void write(XMLWriter& xml, const string& path, const Component_t &comp)
{
  
  push( xml, path );

  write(xml, "color", comp.color);
  write(xml, "spin",  comp.spin);
  
  pop( xml );
}


// Propagator parameters
void read(XMLReader& xml, const string& path, Prop_t& input)
{
  XMLReader inputtop(xml, path);

  read(inputtop, "source_file", input.source_file);
  read(inputtop, "prop_file", input.prop_file);
  read(inputtop, "prop_volfmt", input.prop_volfmt);  // singlefile or multifile
}


// Reader for input parameters
void read(XMLReader& xml, const string& path, PropagatorComponent_input_t& input)
{
  XMLReader inputtop(xml, path);

  // Read the input
  try
  {
    // The parameters holds the version number
    read(inputtop, "Param", input.param);

    // Read in the gauge configuration info
    read(inputtop, "Cfg", input.cfg);

    // Read in the source/propagator info
    read(inputtop, "Prop", input.prop);

    read(inputtop, "Components", input.components);
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading data: " << e << endl;
    throw;
  }
}

/*!
 * Main program for propagator generation. 
 */

int main(int argc, char **argv)
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  // Input parameter structure
  PropagatorComponent_input_t  input;

  // Instantiate xml reader for DATA
  XMLReader xml_in("DATA");

  // Read data
  read(xml_in, "/multiPropagatorComp", input);

  // Specify lattice size, shape, etc.
  Layout::setLattSize(input.param.nrow);
  Layout::create();

  QDPIO::cout << "multiPropagatorComp" << endl;

  // Read in the configuration along with relevant information.
  multi1d<LatticeColorMatrix> u(Nd);
  XMLReader gauge_file_xml, gauge_xml;

  switch (input.cfg.cfg_type) 
  {
  case CFG_TYPE_SZIN :
    readSzin(gauge_xml, u, input.cfg.cfg_file);
    break;

  case CFG_TYPE_SZINQIO:
    readGauge(gauge_file_xml, gauge_xml, u, input.cfg.cfg_file, QDPIO_SERIAL);
    break;

  case CFG_TYPE_NERSC:
    readArchiv(gauge_xml, u, input.cfg.cfg_file);
    break;
  default :
    QDP_error_exit("Configuration type is unsupported.");
  }


  // Read in the source along with relevant information.
  LatticePropagator quark_prop_source;
  XMLReader source_file_xml, source_record_xml;

  // ONLY SciDAC mode is supported for propagators!!
  readQprop(source_file_xml, 
	    source_record_xml, quark_prop_source,
	    input.prop.source_file, QDPIO_SERIAL);

  // Try to invert this record XML into a source struct
  // Also pull out the id of this source
  PropSource_t source_header;

  try
  {
    read(source_record_xml, "/MakeSource/PropSource", source_header);
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error extracting source_header: " << e << endl;
    throw;
  }


  // Instantiate XML writer for XMLDAT
  XMLFileWriter xml_out("XMLDAT");
  push(xml_out, "collectPropcomp");

  proginfo(xml_out);    // Print out basic program info

  // Write out the input
  write(xml_out, "Input", xml_in);

  // Write out the config header
  write(xml_out, "Config_info", gauge_xml);

  // Write out the source header
  write(xml_out, "Source_file_info", source_file_xml);
  write(xml_out, "Source_record_info", source_record_xml);

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

  // Sanity check - write out the norm2 of the source in the Nd-1 direction
  // Use this for any possible verification
  {
    // Initialize the slow Fourier transform phases
    SftMom phases(0, true, Nd-1);

    multi1d<Double> source_corr = sumMulti(localNorm2(quark_prop_source), 
					   phases.getSet());

    push(xml_out, "Source_correlator");
    write(xml_out, "source_corr", source_corr);
    pop(xml_out);
  }

  xml_out.flush();

  // Loop over the source color and spin, creating the source
  // and calling the relevant propagator routines. The QDP
  // terminology is that a propagator is a matrix in color
  // and spin space
  //
  int num_mass = input.param.MultiMasses.size();

  LatticeFermion psi;
  LatticePropagator quark_prop;

  XMLReader file_xml;
  XMLReader record_xml;

  for(int m =0; m < num_mass; m++) { 
    for(int spin=0; spin < Ns; spin++) { 
      for(int color=0; color < Nc; color++) {
	ostringstream filename ;
	filename << input.prop.prop_file << "_component_s" << spin
		 << "_c" << color << "_" 
		 << setw(3) << setfill('0') << m;
	  
	QDPIO::cout << "Attempting to read " << filename.str() << endl;
	  
	// Write the source
	readFermion(file_xml, record_xml, psi,
		    filename.str(), QDPIO_SERIAL);
	
	FermToProp(psi, quark_prop, color, spin);
      }
    }

    SftMom phases(0, true, Nd-1);

    multi1d<Double> prop_corr = sumMulti(localNorm2(quark_prop), 
					 phases.getSet());
    
    push(xml_out, "Prop_correlator");
    write(xml_out, "Number", m);
    write(xml_out, "Mass", input.param.MultiMasses[m]);
    write(xml_out, "prop_corr", prop_corr);
    pop(xml_out);
    
    
    xml_out.flush();

    XMLBufferWriter file_xml;
    push(file_xml, "propagator");
    int id = 0;    // NEED TO FIX THIS - SOMETHING NON-TRIVIAL NEEDED
    write(file_xml, "id", id);
    pop(file_xml);

    

    XMLBufferWriter record_xml;
    push(record_xml, "Propagator");

    // Jiggery pokery. Substitute the ChromaMultiProp_t with a 
    // ChromaProp. This is a pisser because of the FermActParams
    // THIS IS NOT TOTALLY KOSHER AS IT CHANGES THE MASS IN INPUT
    // PARAM as well. However, at this stage we have no further need
    // for input param.
    // I Will eventually write Copy Constructors.

    ChromaProp_t out_param(input.param, m);
 
    write(record_xml, "ForwardProp", out_param);
    write(record_xml, "PropSource", source_header);
    write(record_xml, "Config_info", gauge_xml);
    pop(record_xml);

    ostringstream outfile;
    outfile << input.prop.prop_file << "_" << setw(3) << setfill('0') << m;

    QDPIO::cout << "Attempting to write " << outfile.str() << endl;
   
    // Write the source
    writeQprop(file_xml, record_xml, quark_prop,
	       outfile.str(), input.prop.prop_volfmt, QDPIO_SERIAL);
  }


    
  pop(xml_out);  // propagator
  
  xml_out.close();
  xml_in.close();
  
  // Time to bolt
  QDP_finalize();
  
  exit(0);
}


