// $Id: propagator.cc,v 1.94 2005-03-03 04:05:18 edwards Exp $
/*! \file
 *  \brief Main code for propagator generation
 */

#include <iostream>
#include <cstdio>

#include "chroma.h"

using namespace Chroma;


//! To insure linking of code, place the registered code flags here
/*! This is the bit of code that dictates what fermacts are in use */

bool linkage_hack()
{
  bool foo = true;

  // All actions
  foo &= WilsonTypeFermActsEnv::registered;

  return foo;
}


/*
 * Input 
 */
struct Prop_t
{
  string          source_file;
  string          prop_file;
  QDP_volfmt_t    prop_volfmt;
};

struct Propagator_input_t
{
  ChromaProp_t     param;
  Cfg_t            cfg;
  std::string      stateInfo;
  Prop_t           prop;
};


// Propagator parameters
void read(XMLReader& xml, const string& path, Prop_t& input)
{
  XMLReader inputtop(xml, path);

  read(inputtop, "source_file", input.source_file);
  read(inputtop, "prop_file", input.prop_file);
  read(inputtop, "prop_volfmt", input.prop_volfmt);  // singlefile or multifile
}


// Reader for input parameters
void read(XMLReader& xml, const string& path, Propagator_input_t& input)
{
  XMLReader inputtop(xml, path);

  // Read the input
  try
  {
    // The parameters holds the version number
    read(inputtop, "Param", input.param);

    // Read in the gauge configuration info
    read(inputtop, "Cfg", input.cfg);

    // Read any auxiliary state information
    if( inputtop.count("StateInfo") == 1 ) {
      XMLReader xml_state_info(inputtop, "StateInfo");
      std::ostringstream os;
      xml_state_info.print(os);
      input.stateInfo = os.str();
    }
    else { 
      XMLBufferWriter s_i_xml;
      push(s_i_xml, "StateInfo");
      pop(s_i_xml);
      input.stateInfo = s_i_xml.printCurrentContext();
    }

    // Read in the source/propagator info
    read(inputtop, "Prop", input.prop);
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading data: " << e << endl;
    throw;
  }
}

//! Propagator generation
/*! \defgroup propagator Propagator generation
 *  \ingroup main
 *
 * Main program for propagator generation. 
 */

int main(int argc, char **argv)
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);

  QDPIO::cout << "linkage=" << linkage_hack() << endl;

  START_CODE();

  // Input parameter structure
  Propagator_input_t  input;

  // Instantiate xml reader for DATA
  XMLReader xml_in(Chroma::getXMLInputFileName());

  // Read data
  read(xml_in, "/propagator", input);

  // Specify lattice size, shape, etc.
  Layout::setLattSize(input.param.nrow);
  Layout::create();

  QDPIO::cout << "Propagator" << endl;

  // Read in the configuration along with relevant information.
  multi1d<LatticeColorMatrix> u(Nd);
  XMLReader gauge_file_xml, gauge_xml;


  // Start up the gauge field
  gaugeStartup(gauge_file_xml, gauge_xml, u, input.cfg);


  //
  // Read in the source along with relevant information.
  // 
  LatticePropagator quark_prop_source;
  XMLReader source_file_xml, source_record_xml;
  int t0;
  int j_decay;
  bool make_sourceP = false;
  bool seqsourceP = false;
  {
    // ONLY SciDAC mode is supported for propagators!!
    QDPIO::cout << "Attempt to read source" << endl;
    readQprop(source_file_xml, 
	      source_record_xml, quark_prop_source,
	      input.prop.source_file, QDPIO_SERIAL);
    QDPIO::cout << "Source successfully read" << flush << endl;

    // Try to invert this record XML into a source struct
    try
    {
      // First identify what kind of source might be here
      if (source_record_xml.count("/MakeSource") != 0)
      {
	PropSource_t source_header;

	read(source_record_xml, "/MakeSource/PropSource", source_header);
	j_decay = source_header.j_decay;
	t0 = source_header.t_source[j_decay];
	make_sourceP = true;
      }
      else if (source_record_xml.count("/SequentialSource") != 0)
      {
	ChromaProp_t prop_header;
	PropSource_t source_header;
	SeqSource_t seqsource_header;

	read(source_record_xml, "/SequentialSource/SeqSource", seqsource_header);
	// Any source header will do for j_decay
	read(source_record_xml, "/SequentialSource/ForwardProps/elem[1]/ForwardProp", 
	     prop_header);
	read(source_record_xml, "/SequentialSource/ForwardProps/elem[1]/PropSource", 
	     source_header);
	j_decay = source_header.j_decay;
	t0 = seqsource_header.t_sink;
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
  // XMLFileWriter xml_out("XMLDAT");
  // XMLFileWriter xml_out(Chroma::getXMLOutputFileName());
  XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();
  push(xml_out, "propagator");

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

  //
  // Loop over the source color and spin, creating the source
  // and calling the relevant propagator routines. The QDP
  // terminology is that a propagator is a matrix in color
  // and spin space
  //
  LatticePropagator quark_propagator;
  int ncg_had = 0;

  //
  // Initialize fermion action
  //
  std::istringstream  xml_s(input.param.fermact);
  XMLReader  fermacttop(xml_s);
  const string fermact_path = "/FermionAction";
  string fermact;

  try
  {
    read(fermacttop, fermact_path + "/FermAct", fermact);
  }
  catch (const std::string& e) 
  {
    QDPIO::cerr << "Error reading fermact: " << e << endl;
    throw;
  }

  QDPIO::cout << "FermAct = " << fermact << endl;


  // Deal with auxiliary (and polymorphic) state information
  // eigenvectors, eigenvalues etc. The XML for this should be
  // stored as a string called "stateInfo" in the param struct.

  // Make a reader for the stateInfo
  std::istringstream state_info_is(input.stateInfo);
  XMLReader state_info_xml(state_info_is);
  string state_info_path="/StateInfo";
  //
  // Try the factories
  //
  bool success = false;
  bool mresP = true;

  if (! success)
  {
    try
    {
      QDPIO::cout << "Try the various factories" << endl;

      // Generic Wilson-Type stuff
      Handle< FermionAction<LatticeFermion> >
	S_f(TheFermionActionFactory::Instance().createObject(fermact,
							     fermacttop,
							     fermact_path));


      Handle<const ConnectState> state(S_f->createState(u,
							state_info_xml,
							state_info_path));  // uses phase-multiplied u-fields

      QDPIO::cout << "Suitable factory found: compute the quark prop" << endl;

      S_f->quarkProp(quark_propagator, xml_out, quark_prop_source,
		     t0, j_decay,
		     state, 
		     input.param.invParam, 
		     input.param.nonRelProp,
		     mresP,
		     ncg_had);
      
      success = true;
    }
    catch (const std::string& e) 
    {
      QDPIO::cout << "4D: " << e << endl;
    }
  }


  if (! success)
  {
    QDPIO::cerr << "Error: no fermact found" << endl;
    QDP_abort(1);
  }


  push(xml_out,"Relaxation_Iterations");
  write(xml_out, "ncg_had", ncg_had);
  pop(xml_out);

  // Sanity check - write out the propagator (pion) correlator in the Nd-1 direction
  {
    // Initialize the slow Fourier transform phases
    SftMom phases(0, true, Nd-1);

    multi1d<Double> prop_corr = sumMulti(localNorm2(quark_propagator), 
					 phases.getSet());

    push(xml_out, "Prop_correlator");
    write(xml_out, "prop_corr", prop_corr);
    pop(xml_out);
  }

  xml_out.flush();

  // Save the propagator
  // ONLY SciDAC output format is supported!
  {
    XMLBufferWriter file_xml;
    push(file_xml, "propagator");
    int id = 0;    // NEED TO FIX THIS - SOMETHING NON-TRIVIAL NEEDED
    write(file_xml, "id", id);
    pop(file_xml);

    XMLBufferWriter record_xml;
    if (make_sourceP)
    {
      XMLReader xml_tmp(source_record_xml, "/MakeSource");

      push(record_xml, "Propagator");
      write(record_xml, "ForwardProp", input.param);
      record_xml << input.stateInfo;  // write out the stateinfo - might be empty
      record_xml << xml_tmp;  // write out all the stuff under MakeSource
      pop(record_xml);
    } 
    else if (seqsourceP)
    {
      XMLReader xml_tmp(source_record_xml, "/SequentialSource");

      push(record_xml, "SequentialProp");
      write(record_xml, "SeqProp", input.param);
      record_xml << xml_tmp;  // write out all the stuff under SequentialSource
      pop(record_xml);
    }

    // Write the source
    writeQprop(file_xml, record_xml, quark_propagator,
	       input.prop.prop_file, input.prop.prop_volfmt, QDPIO_SERIAL);

    QDPIO::cout << "Propagator successfully written" << endl;
  }

  pop(xml_out);  // propagator

  // xml_out.close();
  xml_in.close();

  END_CODE();

  // Time to bolt
  Chroma::finalize();

  exit(0);
}
