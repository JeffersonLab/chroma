// $Id: t_propagator_s.cc,v 1.21 2004-11-12 11:00:58 mcneile Exp $
/*! \file
 *  \brief Main code for propagator generation
 */

#include <iostream>
#include <cstdio>

#define MAIN

// Include everything...
#include "chroma.h"

/*
 *  Here we have various temporary definitions
 */
/*
enum CfgType {
  CFG_TYPE_MILC = 0,
  CFG_TYPE_NERSC,
  CFG_TYPE_SCIDAC,
  CFG_TYPE_SZIN,
  CFG_TYPE_UNKNOWN
} ;

enum PropType {
  PROP_TYPE_SCIDAC = 2,
  PROP_TYPE_SZIN,
  PROP_TYPE_UNKNOWN
} ;

enum FermType {
  FERM_TYPE_STAGGERED,
  FERM_TYPE_UNKNOWN
};
*/


using namespace QDP;


/*
 * Input 
 */


// Parameters which must be determined from the XML input
// and written to the XML output
struct Param_t
{
  FermType     FermTypeP;
  Real         Mass;      // Staggered mass
  Real         u0;        // Tadpole Factor
 
  CfgType  cfg_type;       // storage order for stored gauge configuration
  PropType prop_type;      // storage order for stored propagator

  InvertParam_t  invParam;

  Real GFAccu, OrPara;    // Gauge fixing tolerance and over-relaxation param
  int GFMax;              // Maximum gauge fixing iterations

  multi1d<int> nrow;
  multi1d<int> boundary;
  multi1d<int> t_srce; 
};


struct Prop_t
{
  string       source_file;
  string       prop_file;
};

struct Propagator_input_t
{
  Param_t          param;
  Cfg_t            cfg;
  Prop_t           prop;
};


//
void read(XMLReader& xml, const string& path, Prop_t& input)
{
  XMLReader inputtop(xml, path);

//  read(inputtop, "source_file", input.source_file);
  read(inputtop, "prop_file", input.prop_file);
}



// Reader for input parameters
void read(XMLReader& xml, const string& path, Propagator_input_t& input)
{
  XMLReader inputtop(xml, path);


  // First, read the input parameter version.  Then, if this version
  // includes 'Nc' and 'Nd', verify they agree with values compiled
  // into QDP++

  // Read in the IO_version
  int version;
  try
  {
    read(inputtop, "IO_version/version", version);
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading data: " << e << endl;
    throw;
  }


  // Currently, in the supported IO versions, there is only a small difference
  // in the inputs. So, to make code simpler, extract the common bits 

  // Read the uncommon bits first
  try
  {
    XMLReader paramtop(inputtop, "param"); // push into 'param' group

    switch (version) 
    {
      /**************************************************************************/
    case 2 :
      /**************************************************************************/
      break;

    default :
      /**************************************************************************/

      QDPIO::cerr << "Input parameter version " << version << " unsupported." << endl;
      QDP_abort(1);
    }
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading data: " << e << endl;
    throw;
  }


  // Read the common bits
  try 
  {
    XMLReader paramtop(inputtop, "param"); // push into 'param' group

    {
      string ferm_type_str;
      read(paramtop, "FermTypeP", ferm_type_str);
      if (ferm_type_str == "STAGGERED") {
	input.param.FermTypeP = FERM_TYPE_STAGGERED;
      } 
    }

    // GTF NOTE: I'm going to switch on FermTypeP here because I want
    // to leave open the option of treating masses differently.
    switch (input.param.FermTypeP) {
    case FERM_TYPE_STAGGERED :

      QDPIO::cout << " PROPAGATOR: Propagator for Staggered fermions" << endl;

      read(paramtop, "Mass", input.param.Mass);
      read(paramtop, "u0" , input.param.u0);

#if 0
      for (int i=0; i < input.param.numKappa; ++i) {
	if (toBool(input.param.Kappa[i] < 0.0)) {
	  QDPIO::cerr << "Unreasonable value for Kappa." << endl;
	  QDPIO::cerr << "  Kappa[" << i << "] = " << input.param.Kappa[i] << endl;
	  QDP_abort(1);
	} else {
	  QDPIO::cout << " Spectroscopy Kappa: " << input.param.Kappa[i] << endl;
	}
      }
#endif

      break;

    default :
      QDP_error_exit("Fermion type not supported\n.");
    }

    // This may need changing in the future but since the plan is
    // not to write out propagators we'll leave alone for now.
    {
      string prop_type_str;
      read(paramtop, "prop_type", prop_type_str);
      if (prop_type_str == "SZIN") {
	input.param.prop_type = PROP_TYPE_SZIN;
      } else {
	QDP_error_exit("Dont know non SZIN files yet");
      }
    }

//    read(paramtop, "invType", input.param.invType);
    input.param.invParam.invType = CG_INVERTER;   //need to fix this
    read(paramtop, "RsdCG", input.param.invParam.RsdCG);
    read(paramtop, "MaxCG", input.param.invParam.MaxCG);
    read(paramtop, "GFAccu", input.param.GFAccu);
    read(paramtop, "OrPara", input.param.OrPara);
    read(paramtop, "GFMax", input.param.GFMax);

    read(paramtop, "nrow", input.param.nrow);
    read(paramtop, "boundary", input.param.boundary);
    read(paramtop, "t_srce", input.param.t_srce);
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading data: " << e << endl;
    throw;
  }


  // Read in the gauge configuration file name
  try
  {
    read(inputtop, "Cfg", input.cfg);
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

  // Force the machine into a ring for speed up purposes
//  QMP_u32_t foo[] = {1,1,1,8};
//  QMP_declare_logical_topology(foo, Nd);

  // Input parameter structure
  Propagator_input_t  input;

  // Instantiate xml reader for DATA
  XMLReader xml_in("DATA");

  // Read data
  read(xml_in, "/propagator", input);

  // Specify lattice size, shape, etc.
  Layout::setLattSize(input.param.nrow);
  Layout::create();

  // Read in the configuration along with relevant information.
  multi1d<LatticeColorMatrix> u(Nd);
  
  XMLReader gauge_xml;

  switch (input.cfg.cfg_type) 
  {
  case CFG_TYPE_NERSC :
    readArchiv(gauge_xml, u, input.cfg.cfg_file);
    break;
  default :
    QDP_error_exit("Configuration type is unsupported.");
  }

  // Read in the source along with relevant information.
  LatticeStaggeredPropagator quark_prop_source;
  XMLReader source_xml;

  switch (input.param.prop_type) 
  {
  case PROP_TYPE_SZIN :
//    readSzinQprop(source_xml, quark_prop_source, input.prop.source_file);
    quark_prop_source = zero;
    break;
  default :
    QDP_error_exit("Propagator type is unsupported.");
  }


  // Instantiate XML writer for XMLDAT
  XMLFileWriter xml_out("XMLDAT");
  push(xml_out, "hadron_corr");

  // Write out the input
  write(xml_out, "Input", xml_in);

  // Write out the config header
  write(xml_out, "Config_info", gauge_xml);

  // Write out the source header
  write(xml_out, "Source_info", source_xml);

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

  // Fix to the coulomb gauge
  int n_gf;
  int j_decay = Nd-1;

  coulGauge(u, n_gf, j_decay, input.param.GFAccu, input.param.GFMax, true, input.param.OrPara);
  QDPIO::cout << "No. of gauge fixing iterations =" << n_gf << endl;

  // Calcluate plaq on the gauge fixed field
  MesPlq(u, w_plaq, s_plaq, t_plaq, link);
  push(xml_out, "Is_this_gauge_invariant");
  write(xml_out, "w_plaq", w_plaq);
  write(xml_out, "s_plaq", s_plaq);
  write(xml_out, "t_plaq", t_plaq);
  write(xml_out, "link", link);
  pop(xml_out);

  xml_out.flush();

  // Create a fermion BC. Note, the handle is on an ABSTRACT type.
  Handle< FermBC<LatticeStaggeredFermion> >  fbc(new SimpleFermBC<LatticeStaggeredFermion>(input.param.boundary));

  //
  // Initialize fermion action
  //
  AsqtadFermAct S_f(fbc, input.param.Mass, input.param.u0);

  // Set up a state for the current u,
  // (compute fat & triple links)
  // Use S_f.createState so that S_f can pass in u0

  Handle<const ConnectState > state(S_f.createState(u));
  Handle<const EvenOddLinearOperator<LatticeStaggeredFermion> > D_asqtad(S_f.linOp(state));
  Handle<const LinearOperator<LatticeStaggeredFermion> > MdagM_asqtad(S_f.lMdagM(state));

  //
  // Loop over the source color, creating the source
  // and calling the relevant propagator routines. The QDP
  // terminology is that a staggered propagator is a matrix in color space
  // 
  //
  LatticeStaggeredPropagator quark_propagator;
  XMLBufferWriter xml_buf;
  int ncg_had = 0;
  int n_count;
  int t_length = input.param.nrow[3];

  LatticeStaggeredFermion q_source, psi;
  multi1d<LatticeStaggeredPropagator> stag_prop(8);
  multi2d<Real> pi_corr(16, input.param.nrow[3]);
  multi2d<Real> sc_corr(16, input.param.nrow[3]);

  pi_corr = zero;
  sc_corr = zero;

  for(int t_source = 0; t_source < 3; t_source += 2) {

    for(int src_ind = 0; src_ind < 8; ++src_ind){
      psi = zero;   // note this is ``zero'' and not 0
//      int t_source = 0;

      QDPIO::cout << "Inversion for source " << src_ind << endl;

      for(int color_source = 0; color_source < Nc; ++color_source) {
        QDPIO::cout << "Inversion for Color =  " << color_source << endl;
        QDPIO::cout << "Source time slice = " << t_source << endl;
        int spin_source = 0;

        q_source = zero ;
        
//  Use a wall source
	walfil(q_source, t_source, j_decay, color_source, src_ind);

        // Use the last initial guess as the current guess

        // Compute the propagator for given source color/spin 
        // int n_count;

        S_f.qprop(psi, state, q_source, input.param.invParam, n_count);
    
        ncg_had += n_count;
      
        push(xml_out,"Qprop");
        write(xml_out, "Mass" , input.param.Mass);
        write(xml_out, "RsdCG", input.param.invParam.RsdCG);
        write(xml_out, "n_count", n_count);
        pop(xml_out);

        /*
         * Move the solution to the appropriate components
         * of quark propagator.
        */
        FermToProp(psi, quark_propagator, color_source);
      }  //color_source
    
      stag_prop[src_ind] = quark_propagator;
    } // end src_ind
  
    int t_eff;

    staggered_pions pseudoscalar(t_length) ; 
    staggered_scalars  scalar_meson(t_length) ; 
    vector_meson rho(t_length) ; 
     
    pseudoscalar.compute(stag_prop, j_decay);
    scalar_meson.compute(stag_prop,  j_decay);
    rho.compute(stag_prop,  j_decay);
    
    // write the correlators to disk
    pseudoscalar.dump(t_source,xml_out);
    scalar_meson.dump(t_source,xml_out);
    rho.dump(t_source,xml_out);

    // Instantiate XML buffer to make the propagator header
    XMLBufferWriter prop_xml;
    push(prop_xml, "propagator");

    // Write out the input
    write(prop_xml, "Input", xml_in);

    // Write out the config header
    write(prop_xml, "Config_info", gauge_xml);

    // Write out the source header
    write(prop_xml, "Source_info", source_xml);

    pop(prop_xml);


  } //t_source;

  multi1d<Real> Pi(t_length), Sc(t_length);
 
  push(xml_out, "Here_are_all_16_pions_averaged_over_the_2_sources");
  for(int i=0; i < NUM_STAG_PIONS; i++) {
    Pi = pi_corr[i];
    ostringstream tag;
    tag << "pi_corr" << i;
    push(xml_out, tag.str());
    write(xml_out, "Pi", Pi);
    pop(xml_out);
  }
  pop(xml_out);

  push(xml_out, "And_the_same_for_the_scalars");
  for(int i=0; i < NUM_STAG_PIONS; i++) {
    Sc = sc_corr[i];
    ostringstream tag;
    tag << "sc_corr" << i;
    push(xml_out, tag.str());
    write(xml_out, "Sc", Sc);
    pop(xml_out);
  }
  pop(xml_out);

  pop(xml_out);
  xml_out.close();
  xml_in.close();


  // Time to bolt
  QDP_finalize();

  


  exit(0);
}
