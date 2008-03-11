// $Id: t_propagator_s.cc,v 3.10 2008-03-11 21:37:51 mcneile Exp $
/*! \file
 *  \brief Main code for propagator generation
 *
 *  I using this code to DEBUG the HISQ inverter.
 *
 */

#include <iostream>
#include <cstdio>

#define MAIN

// Include everything...
#include "chroma.h"

#include "meas/hadron/stag_propShift_s.h"
#include "meas/hadron/pion_local_s.h"

#include "io/xml_group_reader.h"

/*
 *  Here we have various temporary definitions
 */

using namespace Chroma;


/*
 * Input 
 */


// Parameters which must be determined from the XML input
// and written to the XML output
struct Param_t
{
  Real         Mass;      // Staggered mass
  Real         u0;        // Tadpole Factor
 
  CfgType  cfg_type;       // storage order for stored gauge configuration

  SysSolverCGParams  invParam;
  
   GroupXML_t fermact;

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

  QDPIO::cout << "Propagator for Staggered fermions" << endl;
  //    Read the common bits
   try 
   {
     XMLReader paramtop(inputtop, "param");  // push into 'param' group
       
     // This is a Group XML structure which has 2 members 
     // One is the 'id' which will  be found in the "FermAct" sub-tag
     // One is the root tag of the Group in this case "FermionAction"
     input.param.fermact = readXMLGroup(paramtop, "FermionAction", "FermAct");

     read(paramtop, "RsdCG", input.param.invParam.RsdCG);
     read(paramtop, "MaxCG", input.param.invParam.MaxCG);
     read(paramtop, "GFAccu", input.param.GFAccu);
     read(paramtop, "OrPara", input.param.OrPara);
     read(paramtop, "GFMax", input.param.GFMax);

     read(paramtop, "nrow", input.param.nrow);
     read(paramtop, "t_srce", input.param.t_srce);
   }
   catch (const string& e) 
   {
     QDPIO::cerr << "Error reading data: " << e << endl;
     throw;
   }


   //    Read in the gauge configuration file name
   try
   {
     read(inputtop, "Cfg", input.cfg);
   }
   catch (const string& e) 
   {
     QDPIO::cerr << "Error reading data: " << e << endl;
     throw;
   }
 }


 bool linkageHack(void)
{
  bool foo = true;

  // Inline Measurements
  foo &= StaggeredTypeFermActsEnv::registerAll();
  foo &= InlineAggregateEnv::registerAll();
  foo &= GaugeInitEnv::registerAll();

  return foo;
}
// ! Propagator generation
 /*! \defgroup t_propagator_s Propagator generation
  *  \ingroup testsmain
  *
  * Main program for propagator generation. 
  */

 int main(int argc, char **argv)
 {
   // Put the machine into a known state
   Chroma::initialize(&argc, &argv);

   linkageHack();
   // Input parameter structure
   Propagator_input_t  input;

   // Instantiate xml reader for DATA
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


   // Read the input file
  try
  {
  read(xml_in, "/propagator", input);
  }
    catch (...) 
  {
    QDPIO::cerr << "Error parsing the input file " << in_name << endl;
    throw;
  }
  // Specify lattice size, shape, etc.
  Layout::setLattSize(input.param.nrow);
  Layout::create();

  // Read in the configuration along with relevant information.
  multi1d<LatticeColorMatrix> u(Nd);

  XMLReader  gauge_file_xml,  gauge_xml;
  gaugeStartup(gauge_file_xml, gauge_xml, u, input.cfg);

  // Instantiate XML writer for the output file
  XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();
  push(xml_out, "hadron_corr");

  // Write out the input
  write(xml_out, "Input", xml_in);

  // Write out the config header
  write(xml_out, "Config_info", gauge_xml);

  // Write out the source header

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

  // Calcluate plaq on the gauge fixed field
  MesPlq(xml_out, "Is_this_gauge_invariant", u);
  xml_out.flush();

  // Typedefs to save typing
  typedef LatticeStaggeredFermion      T;
  typedef multi1d<LatticeColorMatrix>  P;
  typedef multi1d<LatticeColorMatrix>  Q;

  //
  // Initialize fermion action
  //

  XMLReader fermact_reader ;
  // Make a memory 'input stream' out of the XML, so we can open an
  // XML Reader on it.
  try{
  std::istringstream is(input.param.fermact.xml);

  // Open a reader on the memory stream.
  //  XMLReader fermact_reader(is);
  fermact_reader.open(is);
  }
  catch (...)
    {
    QDPIO::cerr << "Error reading action name " << endl;
    throw;
    }


  Handle< StaggeredTypeFermAct< T,P,Q> > fermact(TheStagTypeFermActFactory::Instance().createObject(input.param.fermact.id, fermact_reader, input.param.fermact.path));
  // Cast of a pointer to a reference?
  StaggeredTypeFermAct<T,P,Q>& S_f= *(fermact);


  // Set up a state for the current u,
  // (compute fat & triple links)
  // Use S_f.createState so that S_f can pass in u0

  Handle< FermState<T,P,Q> > state(S_f.createState(u));

  //
  // Loop over the source color, creating the source
  // and calling the relevant propagator routines. 
  // 
  //

  LatticeStaggeredPropagator quark_propagator;
  XMLBufferWriter xml_buf;
  int ncg_had = 0;
  int n_count;
  int t_length = input.param.nrow[3];

  LatticeStaggeredFermion q_source, psi;

  push(xml_out, "point_source");
  push(xml_out,"Inverter_properties");
  write(xml_out, "Mass" , input.param.Mass);
  write(xml_out, "RsdCG", input.param.invParam.RsdCG);
  pop(xml_out);

  GroupXML_t inv_param;
  {
    XMLBufferWriter xml_buf;
    write(xml_buf, "InvertParam", input.param.invParam);
    XMLReader xml_in(xml_buf);
    inv_param = readXMLGroup(xml_in, "/InvertParam", "invType");
  }
  Handle< SystemSolver<LatticeStaggeredFermion> > qprop(S_f.qprop(state,inv_param));

  /** do inversions **************************/

  for(int t_source = 0; t_source < 1; t_source += 2) {
    QDPIO::cout << "Source time slice = " << t_source << endl;

      psi = zero;   // note this is ``zero'' and not 0
      multi1d<int>  coord(4)   ;  
      coord[0] = 0 ; 
      coord[1] = 0 ; 
      coord[2] = 0 ; 
      coord[3] =  t_source ; 

      for(int color_source = 0; color_source < Nc; ++color_source) {
        QDPIO::cout << "Inversion for Color =  " << color_source << endl;

        q_source = zero ;
        
	//  Use a point source
	srcfil(q_source, coord, color_source);

        // Use the last initial guess as the current guess

        // Compute the propagator for given source color/spin 
	SystemSolverResults_t res = (*qprop)(psi, q_source);
        ncg_had += res.n_count;
      
	push(xml_out,"Inverter_performance");
        write(xml_out, "color", color_source);
        write(xml_out, "iterations", res.n_count);
        pop(xml_out);

        /*
         * Move the solution to the appropriate components
         * of quark propagator.
        */
        FermToProp(psi, quark_propagator, color_source);
      }  //color_source
    
    int t_eff;

    push(xml_out, "Hadrons_from_time_source");
    write(xml_out, "source_time", t_source);

    staggered_local_pion pion(t_length,u) ;
    pion.compute(quark_propagator,quark_propagator,j_decay) ; 
    pion.dump(t_source,xml_out) ; 

    pop(xml_out);

  } //t_source;
  pop(xml_out);

  pop(xml_out);
  xml_out.close();
  xml_in.close();


  // Time to bolt
  Chroma::finalize();

  exit(0);
}
