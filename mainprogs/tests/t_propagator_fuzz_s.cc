// $Id: t_propagator_fuzz_s.cc,v 1.8 2004-11-20 15:50:25 mcneile Exp $
/*! \file
 *  \brief Main code for propagator generation
 *
 *  Start to add fuzzed source to the staggered 
 *  project.
 */

#include <iostream>
#include <cstdio>

#define MAIN

// Include everything...
#include "chroma.h"

#include "meas/smear/fuzz_smear.h"

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
  IO_version_t     io_version;
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
  try
  {
    read(inputtop, "IO_version/version", input.io_version.version);
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

    switch (input.io_version.version) 
    {
      /**************************************************************************/
    case 1 :
      /**************************************************************************/
      break;

    default :
      /**************************************************************************/

      QDPIO::cerr << "Input parameter version " << input.io_version.version << " unsupported." << endl;
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

    {
      string cfg_type_str;
      read(paramtop, "cfg_type", cfg_type_str);
      if (cfg_type_str == "NERSC") {
	input.param.cfg_type = CFG_TYPE_NERSC;
      } else {
	QDP_error_exit("Dont know non NERSC files yet");
      }

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
    //    read(inputtop, "Cfg", input.cfg);
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

  // Read in the configuration along with relevant information.
  multi1d<LatticeColorMatrix> u(Nd);
  
  XMLReader gauge_xml;

  string nersc = "../../tests/t_asqtad_prop/t_nersc.cfg" ; 
  switch (input.param.cfg_type) 
  {
  case CFG_TYPE_NERSC :
    /*    readArchiv(gauge_xml, u, input.cfg.cfg_file);  */
    readArchiv(gauge_xml, u, nersc );
    break;
  default :
    QDP_error_exit("Configuration type is unsupported.");
  }


  // 
  //  gauge invariance test
  //  

  // this parameter will be read from the input file
  bool do_gauge_transform ;
  do_gauge_transform = false ;
  // do_gauge_transform = true ;


  if( do_gauge_transform )
    {
      // gauge transform the gauge fields
      multi1d<LatticeColorMatrix> u_trans(Nd);

      // create a random gauge transform
       LatticeColorMatrix v ;
  
       gaussian(v);
       reunit(v) ; 

       for(int dir = 0 ; dir < Nd ; ++dir)
	 {
	   u_trans[dir] = v*u[dir]*adj(shift(v,FORWARD,dir)) ;
	   u[dir] = u_trans[dir] ;
	 }

       QDPIO::cout << "Random gauge transform done" << endl;

    } // end of gauge transform



  // Instantiate XML writer for XMLDAT
  XMLFileWriter xml_out("t_propagator_fuzz_s.xml");
  push(xml_out, "fuzzed_hadron_corr");

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

#ifdef NNNNNNNNNNNNNN
  coulGauge(u, n_gf, j_decay, input.param.GFAccu, input.param.GFMax, true, input.param.OrPara);
  QDPIO::cout << "No. of gauge fixing iterations =" << n_gf << endl;
#endif

  // 
  // Ape fuzz the gauge fields
  //

  Real sm_fact = 2.5;   // typical parameter
  int sm_numb = 10;     // number of smearing hits

  int BlkMax = 100;    // max iterations in max-ing trace
  Real BlkAccu = 1.0e-5;  // accuracy of max-ing

  multi1d<LatticeColorMatrix> u_smr(Nd);
  u_smr = u;
  for(int i=0; i < sm_numb; ++i)
  {
    multi1d<LatticeColorMatrix> u_tmp(Nd);

    for(int mu = 0; mu < Nd; ++mu)
      if ( mu != j_decay )
	APE_Smear(u_smr, u_tmp[mu], mu, 0, sm_fact, BlkAccu, BlkMax, j_decay);
      else
	u_tmp[mu] = u_smr[mu];
    
    u_smr = u_tmp;
  }

  //
  //  --- end of APE smearing -----
  //



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

  LatticeStaggeredFermion q_source, psi;
  LatticeStaggeredFermion q_source_fuzz ; 
  multi1d<LatticeStaggeredPropagator> stag_prop(8);

  // the staggered spectroscopy code is hardwired
  // for many pions
  for(int src_ind = 0; src_ind < 0; ++src_ind)
    stag_prop[src_ind] = zero ;


  int t_source = 0;

  // just look at the local pion
  for(int src_ind = 0; src_ind < 1; ++src_ind){
    psi = zero;   // note this is ``zero'' and not 0
    

      for(int color_source = 0; color_source < Nc; ++color_source) {
        QDPIO::cout << "Inversion for Color =  " << color_source << endl;
        QDPIO::cout << "Source time slice = " << t_source << endl;
        int spin_source = 0;

        q_source = zero ;
        

	//  Start to develop fuzzed source code
	enum stag_src_type { LOCAL_SRC , FUZZED_SRC } ;
	//		enum stag_src_type type_of_src = FUZZED_SRC  ;
	enum stag_src_type type_of_src = LOCAL_SRC  ;

	if( type_of_src == LOCAL_SRC )
	  {
        QDPIO::cout << "****> LOCAL SOURCE <****" << endl;

	    q_source = zero ;
	    multi1d<int> coord(Nd);
	    coord[0]=0; coord[1] = 0; coord[2] = 0; coord[3] = 0;
	    srcfil(q_source, coord,color_source ) ;
	  }
	else if( type_of_src == FUZZED_SRC )
	  {
	    int fuzz_width = 2 ; 
	    QDPIO::cout << "***> FUZZED SOURCE ****" << endl;
	    QDPIO::cout << "fuzz width = " << fuzz_width  << endl;

	    q_source = zero ;
	    multi1d<int> coord(Nd);
	    coord[0]=0; coord[1] = 0; coord[2] = 0; coord[3] = 0;
	    srcfil(q_source, coord,color_source ) ;


	    fuzz_smear(u_smr, q_source,q_source_fuzz, 
		       fuzz_width, j_decay) ; 

	    q_source = q_source_fuzz  ;
	  }



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
  
  

  int t_length = input.param.nrow[3];
      push(xml_out, "Hadrons_from_time_source");
      write(xml_out, "source_time", t_source);

      staggered_pions pseudoscalar(t_length) ; 
      pseudoscalar.compute(stag_prop, j_decay);
      pseudoscalar.dump(t_source,xml_out);
      pop(xml_out);

      pop(xml_out);



      //
      // compute some simple baryon operators
      //
      int bc_spec = 0 ;
      multi1d<int> coord(Nd);
      coord[0]=0; coord[1] = 0; coord[2] = 0; coord[3] = 0;

      multi1d<Complex>  barprop(input.param.nrow[3]) ;

      baryon_s(quark_propagator,barprop,
	       coord,j_decay, bc_spec) ;
	       
      push(xml_out, "baryon");
      write(xml_out, "nucleon", barprop);
      pop(xml_out);


      pop(xml_out); // end of document

  xml_out.close();
  xml_in.close();

  // Time to bolt
  QDP_finalize();

  exit(0);
}
