// $Id: propagator.cc,v 1.61 2004-09-09 04:03:10 edwards Exp $
/*! \file
 *  \brief Main code for propagator generation
 */

#include <iostream>
#include <cstdio>

#include "chroma.h"

using namespace QDP;
using namespace Chroma;


// define MRES_CALCULATION in order to run the code computing the residual mass
// and the pseudoscalar to conserved axial current correlator
#define MRES_CALCULATION


//! To insure linking of code, place the registered code flags here
/*! This is the bit of code that dictates what fermacts are in use */
bool linkage_hack()
{
  bool foo = true;

  // 4D actions
  foo &= EvenOddPrecWilsonFermActEnv::registered;
  foo &= UnprecWilsonFermActEnv::registered;

  // 5D actions
  foo &= EvenOddPrecDWFermActArrayEnv::registered;
  foo &= UnprecDWFermActArrayEnv::registered;
  foo &= EvenOddPrecNEFFermActArrayEnv::registered;
  foo &= UnprecNEFFermActArrayEnv::registered;
  
  foo &= UnprecOvDWFermActArrayEnv::registered;
  foo &= UnprecOvExtFermActArrayEnv::registered;

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
  QDP_initialize(&argc, &argv);
 

  START_CODE();

  // Input parameter structure
  Propagator_input_t  input;

  // Instantiate xml reader for DATA
  XMLReader xml_in("DATA");

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
  multi1d<int> boundary;
  bool make_sourceP = false;
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
	PropSource_t source_header;

	read(source_record_xml, "/MakeSource/PropSource", source_header);
	j_decay = source_header.j_decay;
	t0 = source_header.t_source[j_decay];
	boundary = input.param.boundary;
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
	boundary = prop_header.boundary;
	seqsourceP = true;
      }
      else
	throw "No appropriate header found";
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
    for(int i=0; i < boundary.size(); ++i)
      if (boundary[i] != input.param.boundary[i])
      {
	QDPIO::cerr << "Incompatible boundary between input and seqsource" << endl;
	QDP_abort(1);
      }
  }

  // Instantiate XML writer for XMLDAT
  XMLFileWriter xml_out("XMLDAT");
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

  /*
   * Construct fermionic BC. Need one for LatticeFermion and multi1d<LatticeFermion>
   * Note, the handle is on an ABSTRACT type
   */
  Handle< FermBC<LatticeFermion> >  fbc(new SimpleFermBC<LatticeFermion>(input.param.boundary));
  Handle< FermBC<multi1d<LatticeFermion> > >  fbc_a(new SimpleFermBC<multi1d<LatticeFermion> >(input.param.boundary));

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
  catch (const std::exception& e) 
  {
    QDPIO::cerr << "Error reading fermact: " << e.what() << endl;
    throw;
  }

  QDPIO::cout << "FermAct = " << fermact << endl;

  //
  // Try each factory one-by-one
  //
  bool success = false;

  QDPIO::cout << "Try the special cases first" << endl;


#if 1
  if (fermact == EvenOddPrecDWFermActArrayEnv::name)  // FERM_ACT_DWF
  {
    QDPIO::cout << EvenOddPrecDWFermActArrayEnv::name << endl;

    EvenOddPrecDWFermActArray S_f(fbc_a, EvenOddPrecDWFermActArrayParams(fermacttop, fermact_path));
    Handle<const ConnectState> state(S_f.createState(u));  // uses phase-multiplied u-fields

#ifdef MRES_CALCULATION
    if (! input.param.nonRelProp)
      S_f.dwf_quarkProp4(quark_propagator, xml_out, quark_prop_source,
			 t0, j_decay, 
			 state, 
			 input.param.invParam, 
			 ncg_had);
    else
#endif
      S_f.quarkProp4(quark_propagator, xml_out, quark_prop_source,
		     state, 
		     input.param.invParam, 
		     input.param.nonRelProp,
		     ncg_had);
      
    success = true;
  }


  if (fermact == UnprecDWFermActArrayEnv::name)  // FERM_ACT_UNPRECONDITIONED_DWF:
  {
    QDPIO::cout << UnprecDWFermActArrayEnv::name << endl;

    UnprecDWFermActArray S_f(fbc_a, UnprecDWFermActArrayParams(fermacttop, fermact_path));
    Handle<const ConnectState> state(S_f.createState(u));  // uses phase-multiplied u-fields

#ifdef MRES_CALCULATION
    if (! input.param.nonRelProp)
      S_f.dwf_quarkProp4(quark_propagator, xml_out, quark_prop_source,
			 t0, j_decay, 
			 state, 
			 input.param.invParam, 
			 ncg_had);
    else
#endif
      S_f.quarkProp4(quark_propagator, xml_out, quark_prop_source,
		     state, 
		     input.param.invParam, 
		     input.param.nonRelProp,
		     ncg_had);
    
    success = true;
  }
#endif
		      

  QDPIO::cout << "Try the various factories" << endl;

  if (! success)
  {
    try
    {
      // Generic Wilson-Type stuff

      Handle< WilsonTypeFermAct<LatticeFermion> >
	S_f(TheWilsonTypeFermActFactory::Instance().createObject(fermact,
								 fbc,
								 fermacttop,
								 fermact_path));

      Handle<const ConnectState> state(S_f->createState(u));  // uses phase-multiplied u-fields

      S_f->quarkProp4(quark_propagator, xml_out, quark_prop_source,
		      state, 
		      input.param.invParam, 
		      input.param.nonRelProp,
		      ncg_had);
      
      success = true;
    }
    catch (const std::exception& e) 
    {
      QDPIO::cerr << "Error reading data: " << e.what() << endl;
      throw;
    }
  }


  if (! success)
  {
    try
    {
      // Generic 5D Wilson-Type stuff

      Handle< WilsonTypeFermAct<LatticeFermion> >
	S_f(TheWilsonTypeFermActFactory::Instance().createObject(fermact,
								 fbc,
								 fermacttop,
								 fermact_path));

      Handle<const ConnectState> state(S_f->createState(u));  // uses phase-multiplied u-fields

      S_f->quarkProp4(quark_propagator, xml_out, quark_prop_source,
		      state, 
		      input.param.invParam, 
		      input.param.nonRelProp,
		      ncg_had);
      
      success = true;
    }
    catch (const std::exception& e) 
    {
      QDPIO::cerr << "Error reading data: " << e.what() << endl;
      throw;
    }
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

  xml_out.close();
  xml_in.close();

  END_CODE();

  // Time to bolt
  QDP_finalize();

  exit(0);
}
