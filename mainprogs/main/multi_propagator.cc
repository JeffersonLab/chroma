// $Id: multi_propagator.cc,v 1.14 2005-02-28 03:34:46 edwards Exp $
// $Log: multi_propagator.cc,v $
// Revision 1.14  2005-02-28 03:34:46  edwards
// Collapsed code surrounding MesPlq call to a single sub call.
//
// Revision 1.13  2005/01/14 20:13:08  edwards
// Removed all using namespace QDP/Chroma from lib files. The library
// should now be 100% in the Chroma namespace. All mainprogs need a
// using namespace Chroma.
//
// Revision 1.12  2004/12/24 04:19:22  edwards
// Removed explict FermBC args to FermAct factory functions.
//
// Revision 1.11  2004/11/17 15:23:00  bjoo
// t_su3 removed from make check. Throws stringified
//
// Revision 1.10  2004/10/15 16:37:20  bjoo
// Added t_mres4d and lDeltaLs files for 4D m_res calculations
//
// Revision 1.9  2004/07/28 03:08:04  edwards
// Added START/END_CODE to all routines. Changed some to not pass an
// argument.
//
// Revision 1.8  2004/04/28 16:37:53  bjoo
// Adheres to new propagator structure
//
// Revision 1.7  2004/04/26 11:19:13  bjoo
// Added support for Reading EV-s in SZIN format. Must provide e-values in XML tho.
//
// Revision 1.6  2004/04/22 16:49:23  bjoo
// Added overlap_state_info Param struct and gauge startup convenience function. Tidyed up propagator zolo4d case
//
// Revision 1.5  2004/04/22 16:25:25  bjoo
// Added overlap_state_info Param struct and gauge startup convenience function. Tidyed up propagator zolo4d case
//
// Revision 1.4  2004/04/20 13:08:12  bjoo
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

bool linkage_hack()
{
  bool foo = true;

  // 4D actions
  foo &= UnprecWilsonFermActEnv::registered;
  foo &= OvlapPartFrac4DFermActEnv::registered;
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
  ChromaMultiProp_t     param;
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
      input.stateInfo = s_i_xml.str();
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

  // Startup the gauge fields
  gaugeStartup(gauge_file_xml, gauge_xml, u, input.cfg);

  // Read in the source along with relevant information.
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

  // Sanity check
  if (seqsourceP)
  {
    QDPIO::cerr << "Sequential propagator not supportd under multi-mass " << endl;
    QDPIO::cerr << "since source is not mass independent " << endl;
    QDP_abort(1);

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
  int num_mass = input.param.MultiMasses.size();

  multi1d<LatticePropagator> quark_propagator(num_mass);
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

  bool success;
  try {
    QDPIO::cout << "Try the various 4D factories" << endl;
    
    // Generic Wilson-Type stuff
    Handle< WilsonTypeFermAct<LatticeFermion> >
      S_f(TheWilsonTypeFermActFactory::Instance().createObject(fermact,
							       fermacttop,
							       fermact_path));
    
    
    Handle<const ConnectState> state(S_f->createState(u,
						      state_info_xml,
						      state_info_path));  // uses phase-multiplied u-fields
        
    // If this cast fails a bad cast exception is thrown.
    OverlapFermActBase& S_ov = dynamic_cast<OverlapFermActBase&>(*S_f);
    
    
    S_ov.multiQuarkProp4(quark_propagator, 
			 xml_out, 
			 quark_prop_source,
			 state, 
			 input.param.MultiMasses, 
			 input.param.invParam, 
			 1,
			 ncg_had);
		      
  }
  catch(const std::string &e) {
    QDPIO::cout << "4D Wilson like: " << e << endl;

  }
  catch(bad_cast) {
    QDPIO::cerr << "Fermion action created cannot be used for multi mass qprop" << endl;
    QDP_abort(1);
  }


  push(xml_out,"Relaxation_Iterations");
  write(xml_out, "ncg_had", ncg_had);
  pop(xml_out);

  // Sanity check - write out the propagator (pion) correlator in the Nd-1 direction
  {
    // Initialize the slow Fourier transform phases
    SftMom phases(0, true, Nd-1);
    for(int m=0; m < num_mass; m++) { 
      multi1d<Double> prop_corr = sumMulti(localNorm2(quark_propagator[m]), 
					   phases.getSet());

      push(xml_out, "Prop_correlator");
      write(xml_out, "Number", m);
      write(xml_out, "Mass", input.param.MultiMasses[m]);
      write(xml_out, "prop_corr", prop_corr);
      pop(xml_out);
    }
  }

  xml_out.flush();

  // Save the propagators
  // ONLY SciDAC output format is supported!
  for(int m=0; m < num_mass; m++) {
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

    // ChromaProp_t out_param(input.param, m);
 
    ChromaProp_t out_param;
    out_param.nonRelProp = input.param.nonRelProp;
    out_param.fermact = input.param.fermact;
    out_param.invParam.invType = input.param.invParam.invType;
    out_param.invParam.MROver = input.param.invParam.MROver;
    out_param.invParam.MaxCG = input.param.invParam.MaxCG;
    out_param.invParam.RsdCG = input.param.invParam.RsdCG[m];
    out_param.nrow =input.param.nrow ;


    write(record_xml, "ForwardProp", out_param);
    XMLReader xml_tmp(source_record_xml, "/MakeSource");
    record_xml << xml_tmp;
    pop(record_xml);

    ostringstream outfile;
    outfile << input.prop.prop_file << "_" << setw(3) << setfill('0') << m;

    QDPIO::cout << "Attempting to write " << outfile.str() << endl;
   
    // Write the source
    writeQprop(file_xml, record_xml, quark_propagator[m],
	       outfile.str(), input.prop.prop_volfmt, QDPIO_SERIAL);
  }

  pop(xml_out);  // propagator
  
  xml_out.close();
  xml_in.close();
  
  END_CODE();

  // Time to bolt
  QDP_finalize();
  
  exit(0);
}

	  
