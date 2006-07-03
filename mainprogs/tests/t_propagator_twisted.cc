// $Id: t_propagator_twisted.cc,v 3.2 2006-07-03 15:26:11 edwards Exp $
/*! \file
 *  \brief Main code for propagator generation of twisted mass QCD
 *   
 *   This version is for TWISTED Wilson fermions.
 *   This code should work for su3 or su4.
 *
 *   See hep-lat/0411001 for an introduction to twisted mass
 *   QCD. This code is currently being debugged. When the 
 *   code is debugged the input xml file will include 
 *   a test case with a mass term in the gamma_5 component.
 */

#include <iostream>
#include <cstdio>

#define MAIN

// Include everything...
#include "chroma.h"


using namespace Chroma;


/*
 * Input 
 */


// Parameters which must be determined from the XML input
// and written to the XML output
struct Param_t
{
  FermType     FermTypeP;
  Real         u0;        // Tadpole Factor
 
  UnprecParWilsonFermActParams mass_param  ; 

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

  read(inputtop, "prop_file", input.prop_file);
}



// Reader for input parameters
//  first called
//

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

    read(paramtop,"Twisted_mass",input.param.mass_param); 

//    read(paramtop, "invType", input.param.invType);
    read(paramtop, "RsdCG", input.param.invParam.RsdCG);
    read(paramtop, "MaxCG", input.param.invParam.MaxCG);

    read(paramtop, "nrow", input.param.nrow);
    read(paramtop, "boundary", input.param.boundary);
    read(paramtop, "t_srce", input.param.t_srce);

  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading data: " << e << endl;
    throw;
  }

  //
  //   outside <param>  </param>
  //


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
/*! \defgroup t_propagator_twisted Propagator generation
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
  // XMLReader xml_in("INPUT_t_propagator_twisted.xml");
  XMLReader xml_in(Chroma::getXMLInputFileName());

  // Read data
  read(xml_in, "/propagator", input);

  // Specify lattice size, shape, etc.
  Layout::setLattSize(input.param.nrow);
  Layout::create();

  // Read in the configuration along with relevant information.
  multi1d<LatticeColorMatrix> u(Nd);

  QDPIO::cout << "Computation of Meson Correlators for Twisted Mass QCD" << endl;
  QDPIO::cout << "Calculation for SU(" << Nc << ")" << endl;
  XMLReader gauge_file_xml, gauge_xml;
 
  // Start up the gauge field
  gaugeStartup(gauge_file_xml, gauge_xml, u, input.cfg);

  // Check if the gauge field configuration is unitarized
  unitarityCheck(u);

  // Instantiate XML writer for output
  // XMLFileWriter xml_out("t_propagator_twisted.xml");
  XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();
  push(xml_out, "propagator");

  // Write out the input
  write(xml_out, "Input", xml_in);

  // Write out the config header
  write(xml_out, "Config_info", gauge_xml);

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
  bool do_gauge_transform = false ;
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




  
  // -----set up the calculation of quark propagators ---------

  // Create the fermion boundary conditions.
  Handle< FermBC<LatticeFermion> >  fbc(new SimpleFermBC<LatticeFermion>(input.param.boundary));

  //
  // Initialize fermion action
  //
  UnprecParWilsonFermAct  S_f(fbc,input.param.mass_param) ;

  GroupXML_t inv_param;
  {
    XMLBufferWriter xml_buf;
    write(xml_buf, "InvertParam", input.param.invParam);
    XMLReader xml_in(xml_buf);
    inv_param = readXMLGroup(xml_in, "/InvertParam", "invType");
  }

  // Set up a state for the current u,
  Handle<const ConnectState > state(S_f.createState(u));
  Handle<const SystemSolver<LatticeFermion> > qprop(S_f.qprop(state,inv_param));

  LatticePropagator quark_propagator;
  XMLBufferWriter xml_buf;
  int ncg_had = 0;

  LatticeFermion q_source, psi;

  multi1d<int> coord(Nd);
  coord = input.param.t_srce ; 
  //  coord[0]=0; coord[1] = 0; coord[2] = 0; coord[3] = 0;
  int t_source = coord[Nd - 1] ;
  QDPIO::cout << "Source time slice = " << t_source << endl;

  //
  // Loop over the source color and spin , creating the source
  // and calling the relevant propagator routines. 
  //

  for(int color_source = 0; color_source < Nc; ++color_source) 
    for(int spin_source = 0 ; spin_source < Ns ; ++spin_source)
      {
	QDPIO::cout << "Inversion for Color =  " << color_source ; 
	QDPIO::cout << " Spin =  " << spin_source << endl;

	q_source = zero ;
	srcfil(q_source, coord, color_source, spin_source);

	// initial guess is zero
	psi = zero;   // note this is ``zero'' and not 0


        // Compute the propagator for given source color/spin 
	SystemSolverResults_t res = (*qprop)(psi, q_source);     
        ncg_had += res.n_count;
      
        push(xml_out,"Qprop");
        write(xml_out, "Mass" , input.param.mass_param.Mass);
        write(xml_out, "gamma5_mass" , input.param.mass_param.H);
        write(xml_out, "RsdCG", input.param.invParam.RsdCG);
        write(xml_out, "n_count", res.n_count);
        pop(xml_out);

        /*
         * Move the solution to the appropriate components
         * of quark propagator.
        */
        FermToProp(psi, quark_propagator, color_source, spin_source);
      }  //spin / color_source
    
  
  // compute the meson spectrum

  // create averaged Fourier phases with (mom)^2 <= 10
  int j_decay = Nd-1;
  SftMom phases(10, true, j_decay) ;

  mesons(quark_propagator,quark_propagator,
	 phases, t_source, xml_out,
         "Point_Point_Twisted_Wilson_Mesons") ;

  pop(xml_out);

  xml_out.close();
  xml_in.close();

  // Time to bolt
  Chroma::finalize();

  exit(0);
}
