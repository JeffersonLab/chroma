// $Id: t_propagator_s.cc,v 3.6 2008-01-04 15:09:52 mcneile Exp $
/*! \file
 *  \brief Main code for propagator generation
 *
 *  I using this code to DEBUG the HISQ inverter.
 *
 */

#define HISQ_OPTION
#undef  ASQTAD_OPTION

#include <iostream>
#include <cstdio>

#define MAIN

// Include everything...
#include "chroma.h"

#include "meas/hadron/stag_propShift_s.h"
#include "meas/hadron/pion_local_s.h"

#include "actions/ferm/fermacts/hisq_fermact_params_s.h"
#include "actions/ferm/fermacts/hisq_fermact_s.h"
#include "actions/ferm/qprop/hisq_qprop.h"

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
 
  CfgType  cfg_type;       // storage order for stored gauge configuration
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


  //    Read the common bits
   try 
   {
     XMLReader paramtop(inputtop, "param");  // push into 'param' group

     {
       string ferm_type_str;
       read(paramtop, "FermTypeP", ferm_type_str);
       if (ferm_type_str == "STAGGERED") {
 	input.param.FermTypeP = FERM_TYPE_STAGGERED;
       } 
     }

     //      GTF NOTE: I'm going to switch on FermTypeP here because I want
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

/*      This may need changing in the future but since the plan is
	not to write out propagators we'll leave alone for now.
*/
   {
       string prop_type_str;
       read(paramtop, "prop_type", prop_type_str);
       if (prop_type_str == "SZIN") {
 	input.param.prop_type = PROP_TYPE_SZIN;
       } else {
 	QDP_error_exit("Dont know non SZIN files yet");
       }
     }

//     read(paramtop, "invType", input.param.invType);
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


   //    Read in the gauge configuration file name
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

  // Create a fermion state
  Handle< CreateFermState<T,P,Q> > cfs(new CreateSimpleFermState<T,P,Q>(input.param.boundary));

  //
  // Initialize fermion action
  //

  // I know this is an ugly hack from my youth
#if defined(ASQTAD_OPTION)
  QDPIO::cout << "Using ASQTAD inverter " <<  endl;
  write(xml_out, "Action", "ASQTAD");
  AsqtadFermActParams asq_param;
  asq_param.Mass = input.param.Mass;
  asq_param.u0   = input.param.u0;
  AsqtadFermAct S_f(cfs, asq_param);
#else

#if defined(HISQ_OPTION)
  QDPIO::cout << "Using HISQ inverter " <<  endl;
  write(xml_out, "Action", "HISQ");
  HisqFermActParams hisq_param;
  hisq_param.Mass = input.param.Mass;
  hisq_param.u0   = input.param.u0;
  HisqFermAct S_f(cfs, hisq_param);
#else
#error Need to chose HISQ or ASQTAD
#endif

#endif

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
