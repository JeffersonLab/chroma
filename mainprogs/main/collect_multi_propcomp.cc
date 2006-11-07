// $Id: collect_multi_propcomp.cc,v 3.1 2006-11-07 21:51:25 edwards Exp $
/*! \file
 *  \brief Main code for propagator generation
 */

#include <iostream>
#include <cstdio>
#include <iomanip>
#include "chroma.h"

using namespace Chroma;



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
  Chroma::initialize(&argc, &argv);

  START_CODE();

  // Input parameter structure
  PropagatorComponent_input_t  input;

  // Instantiate xml reader for DATA
  XMLReader xml_in(Chroma::getXMLInputFileName());

  // Read data
  read(xml_in, "/multiPropagatorComp", input);

  // Specify lattice size, shape, etc.
  Layout::setLattSize(input.param.nrow);
  Layout::create();

  QDPIO::cout << "multiPropagatorComp" << endl;

  // Read in the configuration along with relevant information.
  multi1d<LatticeColorMatrix> u(Nd);
  XMLReader gauge_file_xml, gauge_xml;

  gaugeStartup(gauge_file_xml, gauge_xml, u, input.cfg);


  // Read in the source along with relevant information.
  LatticePropagator quark_prop_source;
  XMLReader source_file_xml, source_record_xml;

 // ONLY SciDAC mode is supported for propagators!!
  readQprop(source_file_xml, 
	    source_record_xml, quark_prop_source,
	    input.prop.source_file, QDPIO_SERIAL);

  // Try to invert this record XML into a source struct
  // Also pull out the id of this source
  bool make_sourceP = false;;
  bool seqsourceP = false;
  {
    // ONLY SciDAC mode is supported for propagators!!
    QDPIO::cout << "Attempt to read source" << endl;
    readQprop(source_file_xml, 
	      source_record_xml, quark_prop_source,
	      input.prop.source_file, QDPIO_SERIAL);
    QDPIO::cout << "Forward propagator successfully read" << endl;

    // Try to invert this record XML into a source struct
    try
    {
      // First identify what kind of source might be here
      if (source_record_xml.count("/MakeSource") != 0)
      {
	PropSourceConst_t source_header;

	read(source_record_xml, "/MakeSource/PropSource", source_header);
	make_sourceP = true;
      }
      else if (source_record_xml.count("/SequentialSource") != 0)
      {
	SeqSource_t seqsource_header;

	read(source_record_xml, "/SequentialSource/SeqSource", seqsource_header);
	seqsourceP = true;
      }
      else
	throw std::string("No appropriate header found");
    }
    catch (const string& e) 
    {
      QDPIO::cerr << "Error extracting source_header: " << e << endl;
      QDP_abort(1);
    }
  }    

  // Sanity check
  if (seqsourceP)
  {
    QDPIO::cerr << "Sequential propagator not supportd under multi-mass " << endl;
    QDPIO::cerr << "since source is not mass independent " << endl;
    QDP_abort(1);

  }

  // Instantiate XML writer for XMLDAT
//  XMLFileWriter xml_out(Chroma::getXMLOutputFileName());
  XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();
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
  MesPlq(xml_out, "Observables", u);
  xml_out.flush();

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
    XMLReader xml_tmp(source_record_xml, "/MakeSource");
    record_xml << xml_tmp;
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
  
  END_CODE();

  // Time to bolt
  Chroma::finalize();
  
  exit(0);
}


