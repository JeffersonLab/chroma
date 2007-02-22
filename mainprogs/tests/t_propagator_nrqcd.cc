// $Id: t_propagator_nrqcd.cc,v 3.2 2007-02-22 21:11:50 bjoo Exp $
/*! \file
 *  \brief Main code for NRQCD propagator generation
 *   
 *
 *
 *   This is a template program to help the 
 *   Glasgow group impliment the 
 *   NRQCD action to whatever order they want.
 *
 *  See Thacker and Lepage , PRD 43, 1991 ,196
 *
 *
 *   To start this off I will just add in the
 *   additional code to this file. Eventally (but quickly), 
 *   the fermion operator will be put in the 
 *   appropriate place (whatever that is),
 *
 *  The NRQCD evolution equation is not a natural map
 *  to the chroma/qdp++ system. Apply the operators to the
 *  full lattice but only look at a specific time slice,
 *  using the set notation.
 *
 *   For NRQCD code need the field strength (and derivatives#
 *   of.
 * 
 *    This code needs to be converted from SZIN.
 *    ./actions/ferm/fermacts/prec_clover_fermact_w.cc
 *
 *
 *   LATTICE_FIELD_STRENGTH(f);
 *   To get at MesField (u, f);
 *
 *
 */

#include <iostream>
#include <cstdio>

#define MAIN

// Include everything...
#include "chroma.h"


using namespace Chroma;


// ----- start of prototypes for NRQCD -----
void time_evolve(LatticeFermion & Gplus, const LatticeFermion & Gnow, int t ) ;

void compute_nrqcd_prop(LatticeFermion & G, 
			const LatticeFermion & Gsource, 
			const multi1d<LatticeColorMatrix>& u, 
			const  Real Mass,
			int n, int nt ) ;

// ----- start of prototypes for NRQCD -----


// copied from t_ritz.cc
enum GaugeStartType { HOT_START = 0, COLD_START = 1, FILE_START_NERSC = 2 };


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

  GroupXML_t  invParam;

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
      if (ferm_type_str == "WILSON") {
	input.param.FermTypeP = FERM_TYPE_WILSON;
      } 
    }

    // GTF NOTE: I'm going to switch on FermTypeP here because I want
    // to leave open the option of treating masses differently.
    switch (input.param.FermTypeP) {
    case FERM_TYPE_WILSON  :

      QDPIO::cout << " PROPAGATOR: Propagator for NRQCD" << endl;

      read(paramtop, "Mass", input.param.Mass);
      read(paramtop, "u0" , input.param.u0);

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

  // Instantiate xml reader for DATA
  // XMLReader xml_in("INPUT_W.xml");
  XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();

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
    input.cfg.cfg_file = "t_nersc.cfg" ; 
    readArchiv(gauge_xml, u, input.cfg.cfg_file);
    break;
  case COLD_START:
    for(int j = 0; j < Nd; j++) {
      u(j) = Real(1);
    }
    QDPIO::cout << "COLD  unit configuration created" <<  endl;
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
  push(xml_out, "propagator");

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

  // Create a fermion BC. Note, the handle is on an ABSTRACT type.
  Handle< FermBC<LatticeFermion> >  fbc(new SimpleFermBC<LatticeFermion>(input.param.boundary));

  //
  // Initialize fermion action
  //
  UnprecWilsonFermAct S_f(fbc,input.param.Mass);

  // Set up a state for the current u,
  Handle<const ConnectState > state(S_f.createState(u));


  //
  // Loop over the source color and spin , creating the source
  // and calling the relevant propagator routines. 
  // 
  //
  LatticePropagator quark_propagator;
  XMLBufferWriter xml_buf;
  int ncg_had = 0;
  int n_count;

  LatticeFermion q_source, psi;

  int t_source = 0;
#ifdef BBBBBBBBBBBB
  QDPIO::cout << "Source time slice = " << t_source << endl;

  for(int color_source = 0; color_source < Nc; ++color_source) 
    for(int spin_source = 0 ; spin_source < Ns ; ++spin_source)
      {
	QDPIO::cout << "Inversion for Color =  " << color_source << endl;
	QDPIO::cout << "Inversion for Spin =  " << spin_source << endl;


	q_source = zero ;
	multi1d<int> coord(Nd);
	coord[0]=0; coord[1] = 0; coord[2] = 0; coord[3] = 0;
	srcfil(q_source, coord, color_source, spin_source);

	// initial guess is zero
	psi = zero;   // note this is ``zero'' and not 0


        // Compute the propagator for given source color/spin 
        // int n_count;

	const int n_nrqcd = 3 ; // read in from xml
	compute_nrqcd_prop(psi,q_source,u,n_nrqcd ,input.param.Mass,nt);

    
	//        ncg_had += n_count;
      
        /*
         * Move the solution to the appropriate components
         * of quark propagator.
        */
        FermToProp(psi, quark_propagator, color_source, spin_source);
      }  //spin / color_source
    
#endif  

      // compute the meson spectrum

      // source timeslice
      int t0 = 0 ;

      // create averaged Fourier phases with (mom)^2 <= 10
      SftMom phases(10, true, j_decay) ;

      mesons(quark_propagator,quark_propagator,
	     phases, t0, xml_out,
	     "Point_Point_NRQCD_Mesons") ;


  xml_out.close();
  xml_in.close();

  // Time to bolt
  Chroma::finalize();

  exit(0);
}


class TimeSliceFunc : public SetFunc
{
public:
  TimeSliceFunc(int dir): dir_decay(dir) {}
                                                                                
  int operator() (const multi1d<int>& coordinate) const {return coordinate[dir_decay];}
  int numSubsets() const {return Layout::lattSize()[dir_decay];}
                                                                                
  int dir_decay;
                                                                                
private:
  TimeSliceFunc() {}  // hide default constructor
};

/*

 G_[t+1] = U_t\dagger G[t]

 for time slice t


*/

void time_evolve(LatticeFermion & Gplus, 
		 const LatticeFermion & Gnow, 
		 const multi1d<LatticeColorMatrix>& u, 
		 int t )
{
  const int tdir = 3 ;

  Set timeslice;
  timeslice.make(TimeSliceFunc(tdir)) ;

  Subset this_time = timeslice[t+1] ;

  LatticeFermion tmp = shift(adj(u[tdir])*Gnow, BACKWARD, tdir);
  Gplus[this_time] = tmp ;
  //  Gplus = tmp ;  // DEBUG doesn't work

}


/*

  Apply the lowest order NRQCD Hamiltonian to
  a LatticeFermion on a specific timeslice. 


  H =  -1 \sum_i \Delta_i \Delta_-i 
       2M


    Gout = H Gin  on time slice t
*/

void apply_lowest_ke(LatticeFermion & Gout, 
		 const LatticeFermion & Gin, 
		 const multi1d<LatticeColorMatrix>& u, 
		 const  Real Mass,
		 int t )
{
  const int tdir = 3 ;

  Set timeslice;
  timeslice.make(TimeSliceFunc(tdir)) ;


  Subset this_time = timeslice[t] ;

  LatticeFermion  Gout_all ;

  Gout_all = Gin / Mass ;

  LatticeFermion tmp ;

  for(int dir = 0 ; dir < 3 ; ++dir)
    {
      // positive direction
      tmp = Gin ;
      displacement(u,tmp,1,dir);
      Gout_all -= tmp/(2.0*Mass) ; 

      // negative direction
      tmp = Gin ;
      displacement(u,tmp,-1,dir);
      Gout_all -= tmp/(2.0*Mass) ; 


    }


  // restrict to the timeslice
  Gout[this_time] = Gout_all ;

}





//
//  Compute the NRQCD propagator from a 
//  source in Gsource
//
// Impliment equation 28 in the Thacker/Lepage paper
//
//
//

void compute_nrqcd_prop(LatticeFermion & G, 
		 const LatticeFermion & Gsource, 
		 const multi1d<LatticeColorMatrix>& u, 
		 const  Real Mass,
		 int n, int nt )
{

  LatticeFermion  Gold ; 
  LatticeFermion  G_tmp ; 

  Gold = Gsource ;

  for(int t = 0 ; t < nt ; ++t)
    {

      for(int i = 0 ; i < n ; ++i)
	{
	  // apply ( 1 - H /n )
	  apply_lowest_ke(G_tmp,Gold,u,Mass,t);

	  G = G  - G_tmp / n ;
	  Gold =  G ; 
	}

      // evolve forward one time step
      time_evolve(G, Gold, u,t);
      Gold =  G ;

    } 

}
