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

enum VolSrc {
  Z2NOISE ,
  GAUSSIAN
};

typedef   VolSrc  VolSrc_type ; 

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

  bool use_gauge_invar_oper ;

  int Nsamples; 		  // Number of stochastic sources
  int CFGNO;		  // Configuration Number used for seeding rng
                          // This WILL be changed soon 
  Real GFAccu, OrPara;    // Gauge fixing tolerance and over-relaxation param
  int GFMax;              // Maximum gauge fixing iterations

  multi1d<int> nrow;
  multi1d<int> boundary;
  multi1d<int> t_srce; 

  VolSrc_type volume_source ; 

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

      QDPIO::cout << " PROPAGATOR: Disconnected loops for Staggered fermions" << endl;

      read(paramtop, "Mass", input.param.Mass);
      read(paramtop, "u0" , input.param.u0);

      break;

    default :
      QDP_error_exit("Fermion type not supported\n.");
    }

    {
      string prop_type_str;
      read(paramtop, "prop_type", prop_type_str);
      if (prop_type_str == "SZIN") {
	input.param.prop_type = PROP_TYPE_SZIN;
      } else {
	QDP_error_exit("Dont know non SZIN files yet");
      }
    }

    {
      string vol_type_str;
      read(paramtop, "volume_source", vol_type_str);
      if (vol_type_str == "Z2NOISE") 
	{
	  input.param.volume_source = Z2NOISE ;
	} 
      else if  (vol_type_str == "GAUSSIAN") 
	{
	  input.param.volume_source = GAUSSIAN ;
	}
      else
	{
	QDP_error_exit("Wrong type of volume source");
	}


    }



//    read(paramtop, "invType", input.param.invType);
    input.param.invParam.invType = CG_INVERTER;   //need to fix this
    read(paramtop, "RsdCG", input.param.invParam.RsdCG);
    read(paramtop, "MaxCG", input.param.invParam.MaxCG);
    read(paramtop, "Nsamples", input.param.Nsamples);
    read(paramtop, "CFGNO", input.param.CFGNO);
    read(paramtop, "GFAccu", input.param.GFAccu);
    read(paramtop, "OrPara", input.param.OrPara);
    read(paramtop, "GFMax", input.param.GFMax);

    read(paramtop, "nrow", input.param.nrow);
    read(paramtop, "boundary", input.param.boundary);
    read(paramtop, "t_srce", input.param.t_srce);

    read(paramtop, "use_gauge_invar_oper", input.param.use_gauge_invar_oper);

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

  // Input parameter structure
  Propagator_input_t  input;

  // Instantiate xml reader for DATA
  XMLReader xml_in("../../tests/t_asqtad_prop/DISC_DATA_v2");

  // Read data
  read(xml_in, "/propagator", input);

  // Specify lattice size, shape, etc.
  Layout::setLattSize(input.param.nrow);
  Layout::create();

  // Read in the configuration along with relevant information.
  multi1d<LatticeColorMatrix> u(Nd);

  XMLReader  gauge_file_xml,  gauge_xml;
  gaugeStartup(gauge_file_xml, gauge_xml, u, input.cfg);

  // Instantiate XML writer for output
  XMLFileWriter xml_out("t_disc_loop_s.xml");
  push(xml_out, "DISCONNECTED");

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

  // Fix to the coulomb gauge
  int n_gf;
  int j_decay = Nd-1;

  //  coulGauge(u, n_gf, j_decay, input.param.GFAccu, input.param.GFMax, true, input.param.OrPara); 
  //  QDPIO::cout << "No. of gauge fixing iterations =" << n_gf << endl;

  // Calcluate plaq on the gauge fixed field
  MesPlq(u, w_plaq, s_plaq, t_plaq, link);
  push(xml_out, "Is_this_gauge_invariant");
  write(xml_out, "w_plaq", w_plaq);
  write(xml_out, "s_plaq", s_plaq);
  write(xml_out, "t_plaq", t_plaq);
  write(xml_out, "link", link);
  pop(xml_out);

  xml_out.flush();

  // Create the fermion boundary conditions
  Handle< FermBC<LatticeStaggeredFermion> >  fbc(new SimpleFermBC<LatticeStaggeredFermion>(input.param.boundary));

  //
  // Initialize fermion action
  //
  AsqtadFermAct S_f(fbc, input.param.Mass, input.param.u0);
  Handle<const ConnectState > state(S_f.createState(u));
  Handle<const SystemSolver<LatticeStaggeredFermion> > qprop(S_f.qprop(state,input.param.invParam));


  LatticeStaggeredPropagator quark_propagator;
  XMLBufferWriter xml_buf;
  int ncg_had = 0;
  int n_count;

  int Nsamp = input.param.Nsamples;
  int t_length = input.param.nrow[3];
  bool use_gauge_invar_oper = input.param.use_gauge_invar_oper ;

  LatticeStaggeredFermion q_source, psi ;


  // the wrapped disconnected loops
  local_scalar_loop scalar_one_loop(t_length,Nsamp,u) ; 
  non_local_scalar_loop scalar_two_loop(t_length,Nsamp,u) ; 
  threelink_pseudoscalar_loop eta3_loop(t_length,Nsamp,u) ; 
  fourlink_pseudoscalar_loop eta4_loop(t_length,Nsamp,u) ; 


  // This is inefficient for memory 
  multi1d<LatticeStaggeredPropagator> stag_prop(8);
  for(int src_ind = 0; src_ind < 8; ++src_ind)
    stag_prop[src_ind] = zero ;


  int t_source =  input.param.t_srce[Nd-1] ; 

  // Connected Correlator, use a point source
  psi = zero;
  for(int color_source = 0; color_source < Nc; ++color_source) {
    int spin_source = 0;
    q_source = zero;
    srcfil(q_source, input.param.t_srce, color_source);

    // Compute the propagator
    int n_count = (*qprop)(psi, q_source);
    ncg_had += n_count;

    push(xml_out,"Qprop");
    write(xml_out, "Mass" , input.param.Mass);
    write(xml_out, "RsdCG", input.param.invParam.RsdCG);
    write(xml_out, "n_count", n_count);
    pop(xml_out);

    FermToProp(psi, quark_propagator, color_source);
  }

  // compute the connected correlators 

  staggered_pions pseudoscalar(t_length,u) ; 
  staggered_scalars  scalar_meson(t_length,u) ; 
  staggered_pion_singlet pion_singlet(t_length,u); 

  write(xml_out, "use_gauge_invar_oper", use_gauge_invar_oper);
  if( use_gauge_invar_oper )
    {
      cout << "Using gauge invariant operators "  << endl ; 
      pseudoscalar.use_gauge_invar() ;
      scalar_meson.use_gauge_invar() ;
      pion_singlet.use_gauge_invar() ;
    }
  else
    {
      cout << "Using NON-gauge invariant operators "  << endl ; 
      pseudoscalar.use_NON_gauge_invar()  ;
      scalar_meson.use_NON_gauge_invar()  ;
      pion_singlet.use_NON_gauge_invar()  ;
    }
  

  pseudoscalar.compute(stag_prop, j_decay);
  scalar_meson.compute(stag_prop, j_decay);



  // 
  //  Calculate the connected part of the taste singlet pseudoscalar
  //  meson.
  // 
  //  The code should be gauge fixed at this point
  //  or gauge invariant operators should be used.
  //


  LatticeStaggeredPropagator quark_propagator_4link ;

  // Connected Correlator, use a point source at:  1111
  psi = zero;
  for(int color_source = 0; color_source < Nc; ++color_source) {
    int spin_source = 0;
    q_source = zero;
    multi1d<int> coord(Nd);
    coord[0]=1; coord[1] = 1; coord[2] = 1; coord[3] = 1;
    srcfil(q_source, coord,color_source ) ;

    int n_count = (*qprop)(psi, q_source);

    ncg_had += n_count;
    push(xml_out,"Qprop");
    write(xml_out, "Mass" , input.param.Mass);
    write(xml_out, "RsdCG", input.param.invParam.RsdCG);
    write(xml_out, "n_count", n_count);
    pop(xml_out);

    FermToProp(psi, quark_propagator_4link, color_source);
  }


  pion_singlet.compute(quark_propagator,quark_propagator_4link,j_decay); 

  //
  //  write the connected correlators to disk
  // 


  push(xml_out, "CONNECTED");
  pseudoscalar.dump(t_source,xml_out);
  scalar_meson.dump(t_source,xml_out);
  pion_singlet.dump(t_source,xml_out ) ; 
  pop(xml_out);


  //
  // ----- compute disconnected diagrams -----
  //

  // Seed the RNG with the cfg number for now
  Seed seed;
  seed = input.param.CFGNO;
  RNG::setrn(seed);


  for(int i = 0; i < Nsamp; ++i){
    psi = zero;   // note this is ``zero'' and not 0

    // Fill the volume with random noise 

    RNG::savern(seed);
    QDPIO::cout << "SEED = " << seed << endl;
    if( input.param.volume_source == Z2NOISE  )
      gaussian(q_source);
    else if( input.param.volume_source == GAUSSIAN )
      { z2_src(q_source); }

    // Compute the solution vector for the particular source
    int n_count = (*qprop)(psi, q_source);
    ncg_had += n_count;
      
    push(xml_out,"Qprop");
    write(xml_out, "Mass" , input.param.Mass);
    write(xml_out, "RsdCG" , input.param.invParam.RsdCG);
    write(xml_out, "n_count", n_count);
    pop(xml_out);


    scalar_one_loop.compute(q_source,psi,i) ;
    scalar_two_loop.compute(q_source,psi,i) ;
    eta3_loop.compute(q_source,psi,i) ;
    eta4_loop.compute(q_source,psi,i) ;

  } // Nsamples


  // write output from the 
  scalar_one_loop.dump(xml_out) ;
  scalar_two_loop.dump(xml_out) ;
  eta3_loop.dump(xml_out) ;
  eta4_loop.dump(xml_out) ;


  pop(xml_out);

  xml_out.close();
  xml_in.close();

  // Time to bolt
  QDP_finalize();

  exit(0);
}
