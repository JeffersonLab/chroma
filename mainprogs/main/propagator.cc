// $Id: propagator.cc,v 1.56 2004-04-28 02:54:12 edwards Exp $
// $Log: propagator.cc,v $
// Revision 1.56  2004-04-28 02:54:12  edwards
// Added sanity check on boundary.
//
// Revision 1.55  2004/04/27 21:30:01  edwards
// Now supports input from seqsource as well as make_source.
//
// Revision 1.54  2004/04/26 11:19:13  bjoo
// Added support for Reading EV-s in SZIN format. Must provide e-values in XML tho.
//
// Revision 1.53  2004/04/23 11:23:38  bjoo
// Added component based propagator (non multishift) and added Zolotarev5D operator to propagator and propagator_comp. Reworked propagator collection scripts
//
// Revision 1.52  2004/04/22 16:25:25  bjoo
// Added overlap_state_info Param struct and gauge startup convenience function. Tidyed up propagator zolo4d case
//
// Revision 1.51  2004/04/16 22:03:59  bjoo
// Added Zolo 4D test files and tidyed up
//
// Revision 1.50  2004/04/16 14:58:30  bjoo
// Propagator now does Zolo4D
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

#include "chroma.h"

using namespace QDP;


// define MRES_CALCULATION in order to run the code computing the residual mass
// and the pseudoscalar to conserved axial current correlator
#define MRES_CALCULATION

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
  switch (input.param.FermActHandle->getFermActType())
  {
  case FERM_ACT_WILSON:
  {
    QDPIO::cout << "FERM_ACT_WILSON" << endl;

    const WilsonFermActParams& wilson_params=dynamic_cast<const WilsonFermActParams &>(*(input.param.FermActHandle));

    EvenOddPrecWilsonFermAct S_f(fbc,wilson_params.Mass,
				 wilson_params.anisoParam);
    Handle<const ConnectState> state(S_f.createState(u));  // uses phase-multiplied u-fields

    quarkProp4(quark_propagator, xml_out, quark_prop_source,
  	       S_f, state, 
	       input.param.invParam.invType, 
	       input.param.invParam.RsdCG, 
	       input.param.invParam.MaxCG, 
	       input.param.nonRelProp,
	       ncg_had);
  }
  break;

  case FERM_ACT_UNPRECONDITIONED_WILSON:
  {
    QDPIO::cout << "FERM_ACT_UNPRECONDITIONED_WILSON" << endl;

    const WilsonFermActParams& wilson_params=dynamic_cast<const WilsonFermActParams &>(*(input.param.FermActHandle));


    UnprecWilsonFermAct S_f(fbc,wilson_params.Mass);
    Handle<const ConnectState> state(S_f.createState(u));  // uses phase-multiplied u-fields

    quarkProp4(quark_propagator, xml_out, quark_prop_source,
  	       S_f, state, 
	       input.param.invParam.invType, 
	       input.param.invParam.RsdCG, 
	       input.param.invParam.MaxCG, 
	       input.param.nonRelProp,
	       ncg_had);
  }
  break;

  case FERM_ACT_DWF:
  {
    QDPIO::cout << "FERM_ACT_DWF" << endl;


    const DWFFermActParams& dwf_params=dynamic_cast<const DWFFermActParams &>(*(input.param.FermActHandle));

    EvenOddPrecDWFermActArray S_f(fbc_a,
				  dwf_params.chiralParam.OverMass, 
				  dwf_params.Mass, 
				  dwf_params.chiralParam.N5);
    Handle<const ConnectState> state(S_f.createState(u));  // uses phase-multiplied u-fields

#ifdef MRES_CALCULATION
    if (! input.param.nonRelProp)
      dwf_quarkProp4(quark_propagator, xml_out, quark_prop_source,
		     t0, j_decay, 
		     S_f, state, 
		     input.param.invParam.invType, 
		     input.param.invParam.RsdCG, 
		     input.param.invParam.MaxCG, 
		     ncg_had);
    else
#endif
      quarkProp4(quark_propagator, xml_out, quark_prop_source,
		 S_f, state, 
		 input.param.invParam.invType, 
		 input.param.invParam.RsdCG, 
		 input.param.invParam.MaxCG, 
		 input.param.nonRelProp,
		 ncg_had);
  }
  break;

  case FERM_ACT_UNPRECONDITIONED_DWF:
  {
    QDPIO::cout << "FERM_ACT_UNPRECONDITONED_DWF" << endl;

    const DWFFermActParams& dwf_params=dynamic_cast<const DWFFermActParams &>(*(input.param.FermActHandle));

    UnprecDWFermActArray S_f(fbc_a,
			     dwf_params.chiralParam.OverMass, 
			     dwf_params.Mass, 
			     dwf_params.chiralParam.N5);
    Handle<const ConnectState> state(S_f.createState(u));  // uses phase-multiplied u-fields

#ifdef MRES_CALCULATION
    if (! input.param.nonRelProp)
      dwf_quarkProp4(quark_propagator, xml_out, quark_prop_source,
		     t0, j_decay, 
		     S_f, state, 
		     input.param.invParam.invType, 
		     input.param.invParam.RsdCG, 
		     input.param.invParam.MaxCG, 
		     ncg_had);
    else
#endif
      quarkProp4(quark_propagator, xml_out, quark_prop_source,
		 S_f, state, 
		 input.param.invParam.invType, 
		 input.param.invParam.RsdCG, 
		 input.param.invParam.MaxCG, 
		 input.param.nonRelProp,
		 ncg_had);

  }
  break;


  case FERM_ACT_OVERLAP_DWF:
  {
    QDPIO::cout << "FERM_ACT_OVERLAP_DWF" << endl;
    const DWFFermActParams& dwf_params=dynamic_cast<const DWFFermActParams &>(*(input.param.FermActHandle));

    UnprecOvDWFermActArray S_f(fbc_a,
			       dwf_params.chiralParam.OverMass, 
			       dwf_params.Mass, 
			       dwf_params.chiralParam.N5);
    Handle<const ConnectState> state(S_f.createState(u));  // uses phase-multiplied u-fields

    quarkProp4(quark_propagator, xml_out, quark_prop_source,
  	       S_f, state, 
	       input.param.invParam.invType, 
	       input.param.invParam.RsdCG, 
	       input.param.invParam.MaxCG, 
	       input.param.nonRelProp,
	       ncg_had);
  }
  break;
  
  case FERM_ACT_ZOLOTAREV_4D:
  {
    QDPIO::cout << "FERM_ACT_ZOLOTAREV_4D" << endl;
    const Zolotarev4DFermActParams& zolo4d = dynamic_cast<const Zolotarev4DFermActParams& > (*(input.param.FermActHandle));
      
    // Construct Fermact -- now uses constructor from the zolo4d params
    // struct
    Zolotarev4DFermAct S(fbc, zolo4d, xml_out);

      
    // Make a state. Now calls create state with the state info 
    // params struct. Loads Evalues etc if needed, will in the 
    // future recompute them as needed
    Handle<const ConnectState> state(S.createState(u, zolo4d.StateInfo, xml_out, zolo4d.AuxFermActHandle->getMass() )  );
  
    // Call the propagator... Hooray.
    quarkProp4(quark_propagator, xml_out, quark_prop_source,
	       S, state, 
	       input.param.invParam.invType, 
	       input.param.invParam.RsdCG, 
	       input.param.invParam.MaxCG, 
	       input.param.nonRelProp,
	       ncg_had);	  
		      
  }
  break;
  case FERM_ACT_ZOLOTAREV_5D:
  {
    QDPIO::cout << "FERM_ACT_ZOLOTAREV_5D" << endl;
    const Zolotarev5DFermActParams& zolo5d = dynamic_cast<const Zolotarev5DFermActParams& > (*(input.param.FermActHandle));
      
    // Construct Fermact -- now uses constructor from the zolo4d params
    // struct
    Zolotarev5DFermActArray S(fbc_a, fbc, zolo5d, xml_out);

      
    // Make a state. Now calls create state with the state info 
    // params struct. Loads Evalues etc if needed, will in the 
    // future recompute them as needed
    Handle<const ConnectState> state(S.createState(u, zolo5d.StateInfo, xml_out, zolo5d.AuxFermActHandle->getMass()));
  
    // Call the propagator... Hooray.
    quarkProp4(quark_propagator, xml_out, quark_prop_source,
	       S, state, 
	       input.param.invParam.invType, 
	       input.param.invParam.RsdCG, 
	       input.param.invParam.MaxCG, 
	       input.param.nonRelProp,
	       ncg_had);	  
		      
  }
  break;

  default:
    QDPIO::cerr << "Unsupported fermion action" << endl;
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

  // Time to bolt
  QDP_finalize();

  exit(0);
}
