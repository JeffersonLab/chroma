// $Id: propagator.cc,v 1.51 2004-04-16 22:03:59 bjoo Exp $
// $Log: propagator.cc,v $
// Revision 1.51  2004-04-16 22:03:59  bjoo
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

// Forward declaration
void propagatorZolotarev4D(multi1d<LatticeColorMatrix>& u,
			   const Handle< FermBC<LatticeFermion> >&  fbc,
			   const LatticePropagator& quark_prop_source,
			   LatticePropagator& quark_propagator,
			   const Zolotarev4DFermActParams& zolo4d,
			   const InvertParam_t& invParam,
			   int& ncg_had,
			   XMLWriter& xml_out);
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
	       false,
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
	       false,
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

#ifndef MRES_CALCULATION
    quarkProp4(quark_propagator, xml_out, quark_prop_source,
	       S_f, state, 
	       input.param.invParam.invType, 
	       input.param.invParam.RsdCG, 
	       input.param.invParam.MaxCG, 
	       false,
	       ncg_had);
#else
    dwf_quarkProp4(quark_propagator, xml_out, quark_prop_source,
		   source_header.t_source[source_header.j_decay],
		   source_header.j_decay,
		   S_f, state, 
		   input.param.invParam.invType, 
		   input.param.invParam.RsdCG, 
		   input.param.invParam.MaxCG, 
		   ncg_had);
#endif
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
#ifndef MRES_CALCULATION
    quarkProp4(quark_propagator, xml_out, quark_prop_source,
  	       S_f, state, 
	       input.param.invParam.invType, 
	       input.param.invParam.RsdCG, 
	       input.param.invParam.MaxCG, 
	       false,
	       ncg_had);
#else
    dwf_quarkProp4(quark_propagator, xml_out, quark_prop_source,
		   source_header.t_source[source_header.j_decay],
		   source_header.j_decay,
		   S_f, state, 
		   input.param.invParam.invType, 
		   input.param.invParam.RsdCG, 
		   input.param.invParam.MaxCG, 
		   ncg_had);
#endif
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
	       false,
	       ncg_had);
  }
  break;
  
  case FERM_ACT_ZOLOTAREV_4D:
    {
      QDPIO::cout << "FERM_ACT_ZOLOTAREV_4D" << endl;
      const Zolotarev4DFermActParams& zolo4d = dynamic_cast<const Zolotarev4DFermActParams& > (*(input.param.FermActHandle));

      // Short hand...
      propagatorZolotarev4D(u,fbc,quark_prop_source, quark_propagator, zolo4d, input.param.invParam, ncg_had, xml_out);
		      
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
    push(record_xml, "Propagator");
    write(record_xml, "ForwardProp", input.param);
    write(record_xml, "PropSource", source_header);
//    record_xml << source_file_xml;
//    record_xml << source_record_xml;
    write(record_xml, "Config_info", gauge_xml);
    pop(record_xml);

    // Write the source
    writeQprop(file_xml, record_xml, quark_propagator,
	       input.prop.prop_file, input.prop.prop_volfmt, QDPIO_SERIAL);
  }

  pop(xml_out);  // propagator

  xml_out.close();
  xml_in.close();

  // Time to bolt
  QDP_finalize();

  exit(0);
}


// This is too much and too ugly to be inlined
void propagatorZolotarev4D(multi1d<LatticeColorMatrix>& u,
			   const Handle< FermBC<LatticeFermion> >&  fbc,
			   const LatticePropagator& quark_prop_source,
			   LatticePropagator& quark_propagator,
			   const Zolotarev4DFermActParams& zolo4d,
			   const InvertParam_t& invParam,
			   int& ncg_had,
			   XMLWriter& xml_out)
{

  UnprecWilsonTypeFermAct<LatticeFermion>* S_aux;
  XMLBufferWriter zolo_xml;

  // Get the auxiliary fermAct
  switch( zolo4d.AuxFermActHandle->getFermActType() ) {
  case FERM_ACT_WILSON:
    {
      // Upcast
      const WilsonFermActParams& wils = dynamic_cast<const WilsonFermActParams &>( *(zolo4d.AuxFermActHandle));

      //Get the FermAct
      S_aux = new UnprecWilsonFermAct(fbc, wils.Mass);
      if( S_aux == 0x0 ) { 
	QDPIO::cerr << "Unable to instantiate S_aux " << endl;
	QDP_abort(1);
      }

    }
    break;
  default:
    QDPIO::cerr << "Auxiliary Fermion Action Unsupported" << endl;
    QDP_abort(1);
  }
  
  // Drop AuxFermAct into a Handle immediately.
  // This should free things up at the end
  Handle<UnprecWilsonTypeFermAct<LatticeFermion> >  S_w(S_aux);

  
  if( zolo4d.StateInfo.NWilsVec == 0 ) { 
    QDPIO::cout << "Zolo4D Approx Min = " << zolo4d.StateInfo.ApproxMin << endl;
    QDPIO::cout << "Zolo4D Approx Max = " << zolo4d.StateInfo.ApproxMax << endl;
  }
  else { 
    QDPIO::cout << "Zolo4D Approx range from Eigenvalues" << endl;
  }
  
	  
  // Construct Fermact
  Zolotarev4DFermAct S(fbc, S_w, 
		       zolo4d.Mass,
		       zolo4d.RatPolyDeg, 
		       zolo4d.RsdCGInner,
		       zolo4d.MaxCGInner,
		       zolo_xml,
		       zolo4d.ReorthFreqInner);


  // If NWilsVec == 0 these dont take up much room
  multi1d<Real> lambda_lo(zolo4d.StateInfo.NWilsVec);
  multi1d<LatticeFermion> eigv_lo(zolo4d.StateInfo.NWilsVec);
  Real lambda_hi;
  
  // Get the connect state
  const ConnectState* connect_state_ptr;
  
  if( zolo4d.StateInfo.NWilsVec == 0 ) { 

    // NO EV-s so use ApproxMin and ApproxMax from StateInfo
    connect_state_ptr = S.createState(u, 
				      zolo4d.StateInfo.ApproxMin,
				      zolo4d.StateInfo.ApproxMax);
  }
  else {
    
    // Read Eigenvectors
    ChromaWilsonRitz_t ritz_header;
    readEigen(ritz_header, lambda_lo, eigv_lo, lambda_hi, 
	      zolo4d.StateInfo.eigen_io.eigen_file,
	      zolo4d.StateInfo.NWilsVec,
	      QDPIO_SERIAL);
    
    push(xml_out, "EigenSystem");
    write(xml_out, "OriginalRitzHeader", ritz_header);
    write(xml_out, "lambda_lo", lambda_lo);
    write(xml_out, "lambda_high", lambda_hi);
    
    // Test the low eigenvectors
    Handle< const ConnectState > wils_connect_state = S_w->createState(u);
    Handle< const LinearOperator<LatticeFermion> > H = S_w->gamma5HermLinOp(wils_connect_state);
	      
    multi1d<Double> check_norm(zolo4d.StateInfo.NWilsVec);
    multi1d<Double> check_norm_rel(zolo4d.StateInfo.NWilsVec);
    for(int i=0; i < zolo4d.StateInfo.NWilsVec ; i++) { 
      LatticeFermion Me;
      (*H)(Me, eigv_lo[i], PLUS);
      
      LatticeFermion lambda_e;
      
      lambda_e = lambda_lo[i]*eigv_lo[i];
      LatticeFermion r_norm = Me - lambda_e;
      check_norm[i] = sqrt(norm2(r_norm));
      check_norm_rel[i] = check_norm[i]/fabs(Double(lambda_lo[i]));

      QDPIO::cout << "Eigenpair " << i << " Resid Norm = " 
		  << check_norm[i] << " Resid Rel Norm = " << check_norm_rel[i] << endl;
    }
    write(xml_out, "eigen_norm", check_norm);
    write(xml_out, "eigen_rel_norm", check_norm_rel);

    pop(xml_out); // Eigensystem
 
    connect_state_ptr = S.createState(u,
				      lambda_lo,
				      eigv_lo,
				      lambda_hi);
  }

  // Stuff the pointer into a handle. Now, the handle owns the data.
  Handle<const ConnectState> state(connect_state_ptr);
  
  // That's the end of it
  quarkProp4(quark_propagator, xml_out, quark_prop_source,
	     S, state, 
	     invParam.invType, 
	     invParam.RsdCG, 
	     invParam.MaxCG, 
	     false,
	     ncg_had);	  
}
