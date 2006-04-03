// $Id: collect_propcomp.cc,v 3.0 2006-04-03 04:59:13 edwards Exp $
// $Log: collect_propcomp.cc,v $
// Revision 3.0  2006-04-03 04:59:13  edwards
// Major overhaul of fermion and gauge action interface. Basically,
// all fermacts and gaugeacts now carry out  <T,P,Q>  template parameters. These are
// the fermion type, the "P" - conjugate momenta, and "Q" - generalized coordinates
// in the sense of Hamilton's equations. The fermbc's have been rationalized to never
// be over multi1d<T>. The "createState" within the FermionAction is now fixed meaning
// the "u" fields are now from the coordinate type. There are now "ConnectState" that
// derive into FermState<T,P,Q> and GaugeState<P,Q>.
//
// Revision 2.4  2005/12/15 04:04:22  edwards
// New QuarkSpinType to replace the "bool nonRelProp". Now choices
// of how to compute quark prop - FULL,UPPER,LOWER.
//
// Revision 2.3  2005/11/30 04:47:27  edwards
// Changed PropSource_t to PropSourceConst_t and added a new PropSourceSmear_t.
// Renamed PropSink_t to PropSinkSmear_t .
//
// Revision 2.2  2005/11/08 06:30:59  edwards
// Moved nrow to outer structure.
//
// Revision 2.1  2005/11/08 05:40:49  edwards
// Update to use new source mechanism.
//
// Revision 2.0  2005/09/25 21:04:45  edwards
// Moved to version 2.0
//
// Revision 1.13  2005/03/02 00:44:18  edwards
// Changed to new Chroma initialize/finalize format. Changed
// all XMLReader("DATA") to use a command-line param arg.
// Changed all XMLFileWriter(XMLDAT) to use the singleton instance.
//
// Revision 1.12  2005/02/28 03:34:46  edwards
// Collapsed code surrounding MesPlq call to a single sub call.
//
// Revision 1.11  2005/01/14 20:13:08  edwards
// Removed all using namespace QDP/Chroma from lib files. The library
// should now be 100% in the Chroma namespace. All mainprogs need a
// using namespace Chroma.
//
// Revision 1.10  2005/01/12 15:23:26  bjoo
// Moved the mainprogs to use ChromaInitialize and ChromaFinalize. Howver this doesnt buy us much since the linkage hack cannot be properly hidden at the moment (causes segfaults in propagator) and I need closure about how to deal with default input streams. You do get a TheXMLOutputWriter tho
//
// Revision 1.9  2004/12/24 04:19:22  edwards
// Removed explict FermBC args to FermAct factory functions.
//
// Revision 1.8  2004/11/17 15:23:00  bjoo
// t_su3 removed from make check. Throws stringified
//
// Revision 1.7  2004/09/08 14:41:36  edwards
// Ifdef'ed out old FermActHandle code to get to compile.
//
// Revision 1.6  2004/07/28 03:08:04  edwards
// Added START/END_CODE to all routines. Changed some to not pass an
// argument.
//
// Revision 1.5  2004/04/28 14:56:11  bjoo
// Tested propagator_comp and collect_propcomp
//
// Revision 1.3  2004/04/23 11:23:38  bjoo
// Added component based propagator (non multishift) and added Zolotarev5D operator to propagator and propagator_comp. Reworked propagator collection scripts
//
// Revision 1.2  2004/04/22 16:25:25  bjoo
// Added overlap_state_info Param struct and gauge startup convenience function. Tidyed up propagator zolo4d case
//
// Revision 1.1  2004/04/20 13:08:12  bjoo
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
  ChromaProp_t     param;
  Cfg_t            cfg;
  Prop_t           prop;
  multi1d<Component_t> components;
  multi1d<int>     nrow;
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
  Chroma::initialize(&argc, &argv);

  START_CODE();

  // Input parameter structure
  PropagatorComponent_input_t  input;

  // Instantiate xml reader for DATA
  XMLReader xml_in(Chroma::getXMLInputFileName());

  // Read data
  try { 
    read(xml_in, "/propagatorComp", input);
  }
  catch ( const string& e ) { 
    QDPIO::cerr << "Caught exception " << e << endl;
    QDP_abort(1);
  }

  // Specify lattice size, shape, etc.
  Layout::setLattSize(input.nrow);
  Layout::create();

  QDPIO::cout << "propagatorComp" << endl;

  // Read in the configuration along with relevant information.
  multi1d<LatticeColorMatrix> u(Nd);
  XMLReader gauge_file_xml, gauge_xml;


 // Start up the gauge field
  gaugeStartup(gauge_file_xml, gauge_xml, u, input.cfg);

  // Read in the source along with relevant information.
  LatticePropagator quark_prop_source;
  XMLReader source_file_xml, source_record_xml;

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
  LatticeFermion psi;
  LatticePropagator quark_prop;

  XMLReader file_xml_in;
  XMLReader record_xml_in;

  int start_spin;
  int end_spin;
  if( make_sourceP ) { 
    start_spin = 0;
    end_spin = Ns;
  }
  else if (seqsourceP ) { 
    switch (input.param.quarkSpinType)
    {
    case QUARK_SPIN_TYPE_FULL:
      start_spin = 0;
      end_spin = Ns;
      break;

    case QUARK_SPIN_TYPE_UPPER:
      start_spin = 0;
      end_spin = Ns/2;
      break;

    case QUARK_SPIN_TYPE_LOWER:
      start_spin = Ns/2;
      end_spin = Ns;
      break;
    }
  }

  for(int spin=start_spin; spin < end_spin; spin++) { 
    for(int color=0; color < Nc; color++) {
      ostringstream filename ;
      filename << input.prop.prop_file << "_component_s" << spin
	       << "_c" << color ;
      
      QDPIO::cout << "Attempting to read " << filename.str() << endl;
      
      // Write the source
      readFermion(file_xml_in, record_xml_in, psi,
		  filename.str(), QDPIO_SERIAL);
      
      FermToProp(psi, quark_prop, color, spin);
    }
  }

  if( seqsourceP ) 
  {
    switch (input.param.quarkSpinType)
    {
    case QUARK_SPIN_TYPE_FULL:
      // Do nothing here
      break;

    case QUARK_SPIN_TYPE_UPPER:
    {
      for(int color_source = 0; color_source < Nc ; ++color_source) 
      {
	for(int spin_source = Ns/2; spin_source < Ns; ++spin_source) 
	{ 
         PropToFerm(quark_prop,psi, color_source, spin_source-Ns/2);
         FermToProp((-psi), quark_prop, color_source, spin_source);
	}
      }
    }
    break;

    case QUARK_SPIN_TYPE_LOWER:
    {
      // NOT SURE THIS IS CORRECT
      for(int color_source = 0; color_source < Nc ; ++color_source) 
      {
	for(int spin_source = 0; spin_source < Ns/2; ++spin_source) 
	{ 
         PropToFerm(quark_prop,psi, color_source, spin_source+Ns/2);
         FermToProp((-psi), quark_prop, color_source, spin_source);
	}
      }
    }
    break;
    }  // end switch(quarkSpinType)
  }

  SftMom phases(0, true, Nd-1);

  multi1d<Double> prop_corr = sumMulti(localNorm2(quark_prop), 
					 phases.getSet());
    
#if 0
  // Need to fix FermActHandle here
  push(xml_out, "Prop_correlator");
  write(xml_out, "Mass", input.param.FermActHandle->getMass());
  write(xml_out, "prop_corr", prop_corr);
  pop(xml_out);
#endif
  
  
  xml_out.flush();
  
  XMLBufferWriter file_xml;
  push(file_xml, "propagator");
  int id = 0;    // NEED TO FIX THIS - SOMETHING NON-TRIVIAL NEEDED
  write(file_xml, "id", id);
  pop(file_xml);
  
    
  
  XMLBufferWriter record_xml;
  if( make_sourceP) {
  
  // Jiggery pokery. Substitute the ChromaMultiProp_t with a 
  // ChromaProp. This is a pisser because of the FermActParams
  // THIS IS NOT TOTALLY KOSHER AS IT CHANGES THE MASS IN INPUT
  // PARAM as well. However, at this stage we have no further need
  // for input param.
  // I Will eventually write Copy Constructors.
    XMLReader xml_tmp(source_record_xml, "/MakeSource");
  
    push(record_xml, "Propagator");
    write(record_xml, "ForwardProp", input.param);
    record_xml << xml_tmp;
    pop(record_xml);
  }
  else if (seqsourceP) {
      XMLReader xml_tmp(source_record_xml, "/SequentialSource");

      push(record_xml, "SequentialProp");
      write(record_xml, "SeqProp", input.param);
      record_xml << xml_tmp;  // write out all the stuff under SequentialSource
      pop(record_xml);
    }


  QDPIO::cout << "Attempting to write " << input.prop.prop_file << endl;
  
  // Write the source
  writeQprop(file_xml, record_xml, quark_prop,
	     input.prop.prop_file, input.prop.prop_volfmt, QDPIO_SERIAL);


    
  pop(xml_out);  // propagator
    
  END_CODE();

  // Time to bolt
  Chroma::finalize();
  
  exit(0);
}


