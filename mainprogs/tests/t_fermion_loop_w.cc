// $Id: t_fermion_loop_w.cc,v 3.2 2006-07-03 15:26:11 edwards Exp $
/*! \file
 *  \brief Main code for  generation of disconnected 
 *         loops
 *   
 *   This version is for Wilson fermions.
 *   This code should work for su3 or su4.
 *
 *   I don't fully understand the conventions for kappa values
 *   and masses. I will figure this out when I compare against 
 *   the UKQCD code.
 */

#include <iostream>
#include <cstdio>

#define MAIN

// Include everything...
#include "chroma.h"


using namespace Chroma;


// copied from t_ritz.cc
enum GaugeStartType { HOT_START = 0, COLD_START = 1, FILE_START_NERSC = 2 };


void loops(const LatticeFermion &q_source,
           const LatticeFermion &psi,
           int length,
           XMLWriter& xml_gamma,
           const string& xml_tag) ;

void z2_src(LatticeFermion& a) ;

void z2_src(LatticeFermion& a, int slice, int mu) ;


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
 
  GaugeStartType  cfg_type;       // storage order for stored gauge configuration
  PropType prop_type;      // storage order for stored propagator

  SysSolverCGParams  invParam;

  Real GFAccu, OrPara;    // Gauge fixing tolerance and over-relaxation param
  int GFMax;              // Maximum gauge fixing iterations
  int number_sample ;    // number of z2_noise samples

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
      if (ferm_type_str == "WILSON") {
	input.param.FermTypeP = FERM_TYPE_WILSON;
      } 
    }

    // GTF NOTE: I'm going to switch on FermTypeP here because I want
    // to leave open the option of treating masses differently.
    switch (input.param.FermTypeP) {
    case FERM_TYPE_WILSON  :

      QDPIO::cout << "Compute fermion loops for Wilson fermions" << endl;

      read(paramtop, "Mass", input.param.Mass);
      read(paramtop, "u0" , input.param.u0);
      read(paramtop, "number_sample" , input.param.number_sample);

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
	input.param.cfg_type = FILE_START_NERSC  ;
      }
      else if (cfg_type_str == "HOT") {
	input.param.cfg_type = HOT_START;
      }       else if (cfg_type_str == "COLD") {
	input.param.cfg_type = COLD_START;
      } else {
	QDP_error_exit("Only know NERSC/HOT/COLD files yet");
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
//    input.param.invParam.invType = CG_INVERTER;   //need to fix this
    read(paramtop, "RsdCG", input.param.invParam.RsdCG);
    read(paramtop, "MaxCG", input.param.invParam.MaxCG);
    //    read(paramtop, "GFAccu", input.param.GFAccu);
    //    read(paramtop, "OrPara", input.param.OrPara);
    //    read(paramtop, "GFMax", input.param.GFMax);

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

//! Test fermion loops
/*! \defgroup t_fermion_loop Test fermion loops
 *  \ingroup testsmain
 */

int main(int argc, char **argv)
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);

  // Input parameter structure
  Propagator_input_t  input;

  // Instantiate xml reader for DATA
//  XMLReader xml_in("INPUT_W.xml");
  XMLReader xml_in(Chroma::getXMLInputFileName());

  // Read data
  read(xml_in, "/propagator", input);

  // Specify lattice size, shape, etc.
  Layout::setLattSize(input.param.nrow);
  Layout::create();

  // Read in the configuration along with relevant information.
  multi1d<LatticeColorMatrix> u(Nd);
  
  XMLReader gauge_xml;

  QDPIO::cout << "Calculation for SU(" << Nc << ")" << endl;
  switch (input.param.cfg_type) 
  {
  case FILE_START_NERSC :
    // su3 specific (at the moment)
    readArchiv(gauge_xml, u, input.cfg.cfg_file);
    break;
  case HOT_START :
    // create a hot configuration
    for(int dir = 0 ; dir < Nd ; ++dir)
      {
	gaussian(u[dir]);
	reunit(u[dir]) ; 
      }
    QDPIO::cout << "Hot/Random configuration created" <<  endl;
    break;
  default :
    QDP_error_exit("Configuration type is unsupported.");
  }

  // Check if the gauge field configuration is unitarized
  unitarityCheck(u);

  // Instantiate XML writer for XMLDAT
  XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();
  push(xml_out, "z2_loops");

  // Write out the input
  write(xml_out, "Input", xml_in);

  // Write out the config header
  write(xml_out, "Config_info", gauge_xml);

  // Write out the source header
  //  write(xml_out, "Source_info", source_xml);

  push(xml_out, "Output_version");
  write(xml_out, "out_version", 1);
  pop(xml_out);

  xml_out.flush();


  // Calculate some gauge invariant observables just for info.
  MesPlq(xml_out, "Observables", u);
  xml_out.flush();

  // 
  //  gauge invariance test
  //  

  // this parameter will be read from the input file
  bool do_gauge_transform ;
  do_gauge_transform = false ;
  //  do_gauge_transform = true ;


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


  int j_decay = Nd-1;

  
  // set up the calculation of quark propagators 
  // Typedefs to save typing
  typedef LatticeFermion               T;
  typedef multi1d<LatticeColorMatrix>  P;
  typedef multi1d<LatticeColorMatrix>  Q;


  // Create a fermion BC. Note, the handle is on an ABSTRACT type.
  Handle< CreateFermState<T,P,Q> >  cfs(new SimpleFermBC<T,P,Q>(input.param.boundary));

  //
  // Initialize fermion action
  //
  UnprecWilsonFermAct S_f(cfs,input.param.Mass);

  // Set up a state for the current u,
  Handle< FermState<T,P,Q> > state(S_f.createState(u));
  GroupXML_t inv_param;
  {
    XMLBufferWriter xml_buf;
    write(xml_buf, "InvertParam", input.param.invParam);
    XMLReader xml_in(xml_buf);
    inv_param = readXMLGroup(xml_in, "/InvertParam", "invType");
  }
  Handle< SystemSolver<LatticeFermion> > qprop(S_f.qprop(state,inv_param));

  //
  // Loop over the source color and spin , creating the source
  // and calling the relevant propagator routines. 
  // 
  //
  LatticePropagator quark_propagator;
  XMLBufferWriter xml_buf;
  int ncg_had = 0;
  int n_count;


  multi2d<DComplex> loops(16, input.param.nrow[3]);


  LatticeFermion q_source, psi;


  for(int sample = 0 ; sample < input.param.number_sample ; ++sample)
    {
      QDPIO::cout << "Inversion for sample =  " << sample << endl;
      
      q_source = zero ;
      //      gaussian(q_source);
      z2_src(q_source) ;
      //      z2_src(q_source,0,3) ;

      // DEBUG write out the source
      // write(xml_out, "q_source", q_source); 


      // initial guess is zero
      psi = zero;   // note this is ``zero'' and not 0

      // Compute the propagator for given source color/spin 
      // int n_count;
      
      SystemSolverResults_t res = (*qprop)(psi, q_source);         
      ncg_had += res.n_count;
      
      push(xml_out,"Qprop");
      write(xml_out, "Mass" , input.param.Mass);
      write(xml_out, "RsdCG", input.param.invParam.RsdCG);
      write(xml_out, "Sample" , sample);
      write(xml_out, "n_count", res.n_count);
      pop(xml_out);


      string  xml_tag = "dis_loops" ;
      ::loops(q_source,psi,input.param.nrow[3],xml_out,xml_tag) ;
      
    }  // numnber of samples
    
  
  // write out the loops
  push(xml_out,"loop_diagrams");

  pop(xml_out);


  pop(xml_out);
  xml_out.close();
  xml_in.close();

  // Time to bolt
  Chroma::finalize();

  exit(0);
}
