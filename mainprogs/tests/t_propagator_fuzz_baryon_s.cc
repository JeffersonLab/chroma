// $Id: t_propagator_fuzz_baryon_s.cc,v 3.2 2006-07-03 15:26:11 edwards Exp $
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

// more work -- this should be in chroma.h
#include "meas/smear/fuzz_smear.h"

/*
 *  Here we have various temporary definitions
 */
/*


enum FermType {
FERM_TYPE_STAGGERED,
FERM_TYPE_UNKNOWN
};


using namespace Chroma;


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
 
  PropType prop_type;      // storage order for stored propagator

  SysSolverCGParams  invParam;

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

      break;

    default :
      QDP_error_exit("Fermion type not supported\n.");
    }

//    read(paramtop, "invType", input.param.invType);
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


enum stag_src_enum { LOCAL_SRC  , FUZZED_SRC , GAUGE_INVAR_LOCAL_SOURCE } ;
typedef   stag_src_enum stag_src_type ;

/*
  Function to return the type of staggered
  source.
*/

stag_src_type get_stag_src(XMLReader& xml, const string& path)
{
  stag_src_type ans ; 

  try
  {
    string src_name;
    read(xml, path, src_name);
    if (src_name == "LOCAL_SRC") 
    {
      QDPIO::cout << "****> LOCAL SOURCE <****" << endl;
      ans = LOCAL_SRC ; 
    } 
    else if (src_name == "GAUGE_INVAR_LOCAL_SOURCE") 
    {
      QDPIO::cout << "****> GAUGE INVARIANT LOCAL SOURCE <****" << endl;
      ans = GAUGE_INVAR_LOCAL_SOURCE ; 
    } 
    else if (src_name == "FUZZED_SRC") 
    {
      QDPIO::cout << "***> FUZZED SOURCE ****" << endl;
      ans = FUZZED_SRC ; 
    } 

    else
    {
      QDPIO::cerr << "src_name " << src_name << " out of range " << endl;
      QDP_abort(1);
    }
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading data: " << e << endl;
    throw;
  }

  return ans ; 
}

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

void ks_compute_baryon(string name,
		       LatticeStaggeredPropagator & quark_propagator, 
		       XMLFileWriter & xml_out, 
		       int j_decay, int tlength) ;


void ks_compute_baryon(string name,
		       LatticeStaggeredPropagator & quark_propagator_a, 
		       LatticeStaggeredPropagator & quark_propagator_b, 
		       LatticeStaggeredPropagator & quark_propagator_c, 
		       XMLFileWriter & xml_out, 
		       int j_decay, int tlength) ;


int ks_compute_quark_propagator(LatticeStaggeredFermion & psi,
				stag_src_type type_of_src, 
				int fuzz_width,
				multi1d<LatticeColorMatrix> & u , 
				multi1d<LatticeColorMatrix> & u_smr,
				Handle<const SystemSolver<LatticeStaggeredFermion> > & qprop,
				XMLFileWriter & xml_out,
				Real RsdCG, Real Mass, 
				int j_decay, 
				int src_ind, int color_source) ;


void write_smearing_info(string name, stag_src_type type_of_src,
			 XMLFileWriter &xml_out, int fuzz_width ) ;

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

//! Propagator generation
/*! \defgroup t_propagator_fuzz Propagator generation
 *  \ingroup testsmain
 *
 * Main program for propagator generation. 
 */

int main(int argc, char **argv)
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);

  // Input parameter structure
  Propagator_input_t  input;

  // Get the name of the input file and read its contents

  XMLReader xml_in ; 
  string in_name = Chroma::getXMLInputFileName() ; 
  try
  {
    xml_in.open(in_name);
  }
  catch (...) 
  {
    QDPIO::cerr << "Error reading input file " << in_name << endl;
    QDPIO::cerr << "The input file name can be passed via the -i flag " << endl;
    QDPIO::cerr << "The default name is ./DATA" << endl;
    throw;
  }


  // Read data
  read(xml_in, "/propagator", input);

  // Specify lattice size, shape, etc.
  Layout::setLattSize(input.param.nrow);
  Layout::create();

  // Read in the configuration along with relevant information.
  multi1d<LatticeColorMatrix> u(Nd);
  
  XMLReader  gauge_file_xml,  gauge_xml;
  gaugeStartup(gauge_file_xml, gauge_xml, u, input.cfg);


  int fuzz_width = 2 ; 
  try
  {
    read(xml_in, "/propagator/param/fuzz_width",fuzz_width );
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading fuzzing width " << e << endl;
    throw;
  }
  QDPIO::cout << "fuzz width = " << fuzz_width  << endl;


  bool use_gauge_invar_oper ;
  use_gauge_invar_oper = false ;
  read(xml_in, "/propagator/param/use_gauge_invar_oper",use_gauge_invar_oper );


  // 
  //  gauge invariance test
  //  

  // this parameter will be read from the input file
  bool do_gauge_transform ;
  do_gauge_transform = false ;
  read(xml_in, "/propagator/param/do_gauge_transform",do_gauge_transform );

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



  // Instantiate XML writer for the output
  XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();
  push(xml_out, "fuzzed_hadron_corr");

  // Write out the input parameter file
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
  MesPlq(xml_out, "Observables", u);
  xml_out.flush();

  // Fix to the coulomb gauge
  int n_gf;
  int j_decay = Nd-1;

  if( ! use_gauge_invar_oper )
  {
    QDPIO::cout << "Starting Coulomb gauge fixing" << endl;
    coulGauge(u, n_gf, j_decay, input.param.GFAccu, input.param.GFMax, true, input.param.OrPara);
    QDPIO::cout << "No. of gauge fixing iterations =" << n_gf << endl;
  }

  // 
  // Ape fuzz the gauge fields
  //

  multi1d<LatticeColorMatrix> u_smr(Nd);


  QDPIO::cout << "Starting to APE smear the gauge configuration" << endl;
  
  Real sm_fact = 2.5;   // typical parameter
  int sm_numb = 10;     // number of smearing hits
  
  int BlkMax = 100;    // max iterations in max-ing trace
  Real BlkAccu = 1.0e-5;  // accuracy of max-ing
      
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
  MesPlq(xml_out, "Observables", u);
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

  GroupXML_t inv_param;
  {
    XMLBufferWriter xml_buf;
    write(xml_buf, "InvertParam", input.param.invParam);
    XMLReader xml_in(xml_buf);
    inv_param = readXMLGroup(xml_in, "/InvertParam", "invType");
  }
  Handle<const ConnectState > state(S_f.createState(u));
  Handle<const SystemSolver<LatticeStaggeredFermion> > qprop(S_f.qprop(state,inv_param));


  //
  // Loop over the source color, creating the source
  // and calling the relevant propagator routines.
  // 

  LatticeStaggeredPropagator quark_propagator_Lsink_Lsrc;
  LatticeStaggeredPropagator quark_propagator_Fsink_Lsrc;
  LatticeStaggeredPropagator quark_propagator_Lsink_Fsrc;
  LatticeStaggeredPropagator quark_propagator_Fsink_Fsrc;

  int ncg_had = 0;

  LatticeStaggeredFermion psi ;
  LatticeStaggeredFermion psi_fuzz;

  // pass the origin from input file
  int t_source = 0;
  QDPIO::cout << "Source time slice = " << t_source << endl;



  //
  //  Local source inversions
  //

  stag_src_type type_of_src = LOCAL_SRC ;
  QDPIO::cout << "LOCAL INVERSIONS"  << endl;

  for(int color_source = 0; color_source < Nc; ++color_source) 
  {
    psi = zero;   // note this is ``zero'' and not 0

    const int src_ind = 0 ; 
    ncg_had += ks_compute_quark_propagator(psi,type_of_src, fuzz_width,
					   u, u_smr, qprop, xml_out,
					   input.param.invParam.RsdCG,
					   input.param.Mass,
					   j_decay, 
					   src_ind, color_source) ;

    /*
     * Move the solution to the appropriate components
     * of quark propagator.
     */
    FermToProp(psi, quark_propagator_Lsink_Lsrc, color_source);

    //
    //  fuzz at the sink 
    //

    fuzz_smear(u_smr, psi,psi_fuzz, fuzz_width, j_decay) ;
    FermToProp(psi_fuzz, quark_propagator_Fsink_Lsrc, color_source);


  }  //color_source
    



  //
  //  Fuzzed source inversions
  //

  type_of_src = FUZZED_SRC ;
  QDPIO::cout << "FUZZED SOURCE INVERSIONS"  << endl;

  for(int color_source = 0; color_source < Nc; ++color_source) 
  {
    psi = zero;   // note this is ``zero'' and not 0

    const int src_ind = 0 ; 
    ncg_had += ks_compute_quark_propagator(psi,type_of_src, fuzz_width,
					   u, u_smr, qprop, xml_out,
					   input.param.invParam.RsdCG,
					   input.param.Mass,
					   j_decay, 
					   src_ind, color_source) ;

    /*
     * Move the solution to the appropriate components
     * of quark propagator.
     */
    FermToProp(psi, quark_propagator_Lsink_Fsrc, color_source);

    //
    //  fuzz at the sink 
    //

    fuzz_smear(u_smr, psi,psi_fuzz, fuzz_width, j_decay) ;
    FermToProp(psi_fuzz, quark_propagator_Fsink_Fsrc, color_source);


  }  //color_source
    
      

  
  //
  // compute some simple baryon correlators
  //

  push(xml_out, "baryon_correlators");

  // describe the source
  string NN ;
  write(xml_out, "source_time", t_source);
  push(xml_out, "smearing_info");
  NN = "L" ; 
  write_smearing_info(NN, LOCAL_SRC,xml_out,fuzz_width) ;

  NN = "F" ; 
  write_smearing_info(NN,FUZZED_SRC,xml_out,fuzz_width) ;
    
  pop(xml_out);

  // write out the baryon correlators 
  string b_tag("srcLLL_sinkLLL_nucleon") ;  
  ks_compute_baryon(b_tag,quark_propagator_Lsink_Lsrc, 
		    quark_propagator_Lsink_Lsrc, 
		    quark_propagator_Lsink_Lsrc, 
		    xml_out, j_decay, 
		    input.param.nrow[3]) ;

  // single quark fuzzed

  b_tag = "srcLLL_sinkFLL_nucleon" ;
  ks_compute_baryon(b_tag,
		    quark_propagator_Fsink_Lsrc, 
		    quark_propagator_Lsink_Lsrc, 
		    quark_propagator_Lsink_Lsrc, 
		    xml_out, j_decay, 
		    input.param.nrow[3]) ;

  b_tag = "srcFLL_sinkLLL_nucleon" ;
  ks_compute_baryon(b_tag,
		    quark_propagator_Lsink_Fsrc, 
		    quark_propagator_Lsink_Lsrc, 
		    quark_propagator_Lsink_Lsrc, 
		    xml_out, j_decay, 
		    input.param.nrow[3]) ;

  b_tag = "srcFLL_sinkFLL_nucleon" ;
  ks_compute_baryon(b_tag,
		    quark_propagator_Fsink_Fsrc, 
		    quark_propagator_Lsink_Lsrc, 
		    quark_propagator_Lsink_Lsrc, 
		    xml_out, j_decay, 
		    input.param.nrow[3]) ;



  // double quark fuzzed

  b_tag = "srcLLL_sinkFFL_nucleon" ;
  ks_compute_baryon(b_tag,
		    quark_propagator_Fsink_Lsrc, 
		    quark_propagator_Fsink_Lsrc, 
		    quark_propagator_Lsink_Lsrc, 
		    xml_out, j_decay, 
		    input.param.nrow[3]) ;

  b_tag = "srcFFL_sinkLLL_nucleon" ;
  ks_compute_baryon(b_tag,
		    quark_propagator_Lsink_Fsrc, 
		    quark_propagator_Lsink_Fsrc, 
		    quark_propagator_Lsink_Lsrc, 
		    xml_out, j_decay, 
		    input.param.nrow[3]) ;

  b_tag = "srcFFL_sinkFFL_nucleon" ;
  ks_compute_baryon(b_tag,
		    quark_propagator_Fsink_Fsrc, 
		    quark_propagator_Fsink_Fsrc, 
		    quark_propagator_Lsink_Lsrc, 
		    xml_out, j_decay, 
		    input.param.nrow[3]) ;


  // treble quark fuzzed

  b_tag = "srcLLL_sinkFFF_nucleon" ;
  ks_compute_baryon(b_tag,
		    quark_propagator_Fsink_Lsrc, 
		    quark_propagator_Fsink_Lsrc, 
		    quark_propagator_Fsink_Lsrc, 
		    xml_out, j_decay, 
		    input.param.nrow[3]) ;

  b_tag = "srcFFF_sinkLLL_nucleon" ;
  ks_compute_baryon(b_tag,
		    quark_propagator_Lsink_Fsrc, 
		    quark_propagator_Lsink_Fsrc, 
		    quark_propagator_Lsink_Fsrc, 
		    xml_out, j_decay, 
		    input.param.nrow[3]) ;

  b_tag = "srcFFF_sinkFFF_nucleon" ;
  ks_compute_baryon(b_tag,
		    quark_propagator_Fsink_Fsrc, 
		    quark_propagator_Fsink_Fsrc, 
		    quark_propagator_Fsink_Fsrc, 
		    xml_out, j_decay, 
		    input.param.nrow[3]) ;



  pop(xml_out);  // baryon correlators

      
  pop(xml_out);
  xml_out.close();
  xml_in.close();

  // Time to bolt
  Chroma::finalize();
  QDPIO::cout << "CHROMA_RUN_COMPLETE " << endl;
  exit(0);
}


//
// wrapper routine for baryon operators
//

void ks_compute_baryon(string name,
		       LatticeStaggeredPropagator & quark_propagator, 
		       XMLFileWriter & xml_out, 
		       int j_decay, int tlength)
{
  int bc_spec = 0 ;
  multi1d<int> coord(Nd);
  coord[0]=0; coord[1] = 0; coord[2] = 0; coord[3] = 0;
  
  multi1d<Complex>  barprop(tlength) ;
  
  baryon_s(quark_propagator,barprop,
	   coord,j_decay, bc_spec) ;
  

  write(xml_out, name, barprop);

 

}



//
// wrapper routine for baryon operators
//

void ks_compute_baryon(string name,
		       LatticeStaggeredPropagator & quark_propagator_a, 
		       LatticeStaggeredPropagator & quark_propagator_b, 
		       LatticeStaggeredPropagator & quark_propagator_c, 
		       XMLFileWriter & xml_out, 
		       int j_decay, int tlength)
{
  int bc_spec = 0 ;
  multi1d<int> coord(Nd);
  coord[0]=0; coord[1] = 0; coord[2] = 0; coord[3] = 0;
  
  multi1d<Complex>  barprop(tlength) ;
  
  baryon_s(quark_propagator_a,quark_propagator_b,quark_propagator_c,
	   barprop,coord,j_decay, bc_spec) ;
  
  write(xml_out, name, barprop);

}




int ks_compute_quark_propagator(LatticeStaggeredFermion & psi,
				stag_src_type type_of_src, 
				int fuzz_width,
				multi1d<LatticeColorMatrix> & u , 
				multi1d<LatticeColorMatrix> & u_smr,
				Handle<const SystemSolver<LatticeStaggeredFermion> > & qprop,
				XMLFileWriter & xml_out,
				Real RsdCG, Real Mass, 
				int j_decay, 
				int src_ind, int color_source)
{
  LatticeStaggeredFermion q_source ;
  LatticeStaggeredFermion q_source_fuzz ; 
  int ncg_had = 0 ;

  QDPIO::cout << "Inversion for Color =  " << color_source << endl;
  q_source = zero ;

  if( type_of_src == LOCAL_SRC )
  {
    q_source = zero ;
    multi1d<int> coord(Nd);

    PropIndexTodelta(src_ind, coord) ; 
    srcfil(q_source, coord,color_source ) ;
  }
  else if( type_of_src == GAUGE_INVAR_LOCAL_SOURCE  )
  {
    q_source = zero ;
    multi1d<int> coord(Nd);
	      
    // start with local source 
    coord[0]=0; coord[1] = 0; coord[2] = 0; coord[3] = 0;
    srcfil(q_source, coord,color_source ) ;
	      
    // now do the shift
    PropIndexTodelta(src_ind, coord) ; 
    q_source_fuzz = q_source  ;
    q_source = shiftDeltaPropCov(coord,q_source_fuzz,u,
				 false); 

  }
  else if( type_of_src == FUZZED_SRC )
  {
    q_source = zero ;
    multi1d<int> coord(Nd);

    PropIndexTodelta(src_ind, coord) ; 
    srcfil(q_source, coord,color_source ) ;


    fuzz_smear(u_smr, q_source,q_source_fuzz, 
	       fuzz_width, j_decay) ; 

    q_source = q_source_fuzz  ;
  }



  // Use the last initial guess as the current guess

  // Compute the propagator for given source color/spin 
  // int n_count;

  StopWatch swatch;
  swatch.start();

  SystemSolverResults_t res n_count = (*qprop)(psi, q_source);
  swatch.stop();
  double time_in_sec  = swatch.getTimeInSeconds();

    
  ncg_had += res.n_count;

  // this is done for xmldif reasons
  if( src_ind == 0 )
  {
    push(xml_out,"Qprop");
    write(xml_out, "Staggered_src_tag" , src_ind);
    write(xml_out, "Mass" , Mass);
    write(xml_out, "RsdCG", RsdCG);
    write(xml_out, "n_count", n_count);
    write(xml_out, "time_in_sec",time_in_sec );
    pop(xml_out);
  }


  return ncg_had ;
}



void write_smearing_info(string name, stag_src_type type_of_src,
			 XMLFileWriter &xml_out, int fuzz_width )
{


  push(xml_out, "smearing_basis");
  write(xml_out, "element_tag", name);

  if( type_of_src == LOCAL_SRC )
  { write(xml_out, "source_type", "LOCAL_SRC"); }
  else if( type_of_src == GAUGE_INVAR_LOCAL_SOURCE  )
  { write(xml_out, "source_type", "GAUGE_INVAR_LOCAL_SOURCE"); }
  else if( type_of_src == FUZZED_SRC )
  { 
    write(xml_out, "source_type", "FUZZED_SRC"); 
    write(xml_out, "fuzzed_width", fuzz_width); 
  }

  pop(xml_out);

}
