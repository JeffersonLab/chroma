// $Id: t_stagg_baryon.cc,v 1.1 2009-03-18 16:00:13 mcneile Exp $
/*! \file
 *  \brief Main code for staggered charmed baryons
 *
 *  This code is a wrapper for the calculation of staggered charmed
 *  baryons.
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
#include "meas/hadron/baryon_s.h"


enum staggered_src_type { local_s, wall_s , wall_o , wall_q , wall_o_and_q,
 local_o, local_q , local_o_and_q   };


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

  CfgType  cfg_type;       // storage order for stored gauge configuration

  SysSolverCGParams  invParam;

  GroupXML_t fermact;
  GroupXML_t fermact_A;

  bool  gauge_trans ; 
  staggered_src_type src_type ; 
  std::string src_type_str ; 

  multi1d<int> nrow;
  multi1d<int> boundary;
  int t_srce;
};



struct Propagator_input_t
{
  Param_t          param;
  Cfg_t            cfg;
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
     input.param.fermact_A = readXMLGroup(paramtop, "FermionAction_A", "FermAct_A");

     read(paramtop, "RsdCG", input.param.invParam.RsdCG);
     read(paramtop, "MaxCG", input.param.invParam.MaxCG);
     read(paramtop, "MinCG", input.param.invParam.MinCG);

     read(paramtop, "gauge_trans", input.param.gauge_trans);
     read(paramtop, "src_type", input.param.src_type_str);

     if( input.param.src_type_str.compare("LOCAL") == 0 )
       {
	 input.param.src_type = local_s ; 
       }
     else if( input.param.src_type_str.compare("WALL") == 0 )
       {
	 input.param.src_type = wall_s ; 
       }
     else if( input.param.src_type_str.compare("WALL_Q") == 0 )
       {
	 input.param.src_type = wall_q ; 
       }
     else if( input.param.src_type_str.compare("LOCAL_Q") == 0 )
       {
	 input.param.src_type = local_q ; 
       }
     else 
       {
       QDPIO::cerr << "Input src = " << input.param.src_type_str << " unknown." << endl;
       QDP_abort(1);
       }



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



#include "qdp_util.h"

void dump_text_src(LatticeStaggeredFermion psi, multi1d<int> nrow) 
{

  QDPIO::cout << "Dump of source " << endl;
  StaggeredFermion psi_local ; 
  ColorVector psi_color ;
  multi1d<int> coords ;
  Complex z ; 

  for(int site=0; site < Layout::vol(); ++site)
    {
      multi1d<int> coord = crtesn(site, Layout::lattSize());
      printf("[%d,%d,%d,%d] \n",
	     coord[0],coord[1],coord[2],coord[3] )  ;
      psi_local = peekSite(psi,coord) ;
      psi_color = peekSpin(psi_local,0) ; 
      for(int i=0 ; i < Nc ; ++i)
	  {
	    z = peekColor(psi_color,i) ;
	    QDPIO::cout <<  i << " " <<  z << endl;
	  }
      
    }



}


//
//  Staggered wall source from Gupta paper
//  This is known as an even source
//  This is a copy of the routine in the MILC code.
//

void walfil_q(LatticeStaggeredFermion& a, int slice, int mu, int color_index)
{
  START_CODE();

  if (color_index >= Nc || color_index < 0)
    QDP_error_exit("invalid color index", color_index);
  
  int spin_index = 0 ;

  // Write ONE to all field
  Real nnn = -1/4.0 ;
  Complex sitecomp = cmplx(nnn,0);
  ColorVector sitecolor = zero;
  StaggeredFermion sitefield = zero;
  
  pokeSpin(sitefield,
	   pokeColor(sitecolor,sitecomp,color_index),
	   spin_index);
  
  // Narrow the context to the desired slice.
  LatticeStaggeredFermion tmp;
  tmp = sitefield;  // QDP (not installed version) now supports   construct OLattice = OScalar
  
  // Auxiliary: Coordinates to use in "where" clauses
  multi1d<LatticeInteger> x(Nd);
  // Fill x with lattice coordinates
  for( int sigma = 0; sigma < Nd; sigma++) {
    x[ sigma ] = Layout::latticeCoordinate(sigma);
  }

  a = where(Layout::latticeCoordinate(mu) == slice && 
	    (x[0]+x[1]+x[2] ) % 2 ==  0,
	    tmp, LatticeStaggeredFermion(zero));

  
  END_CODE();
}



//
//  Staggered wall source from Gupta paper
//  This is known as an odd source
//
void walfil_o(LatticeStaggeredFermion& a, int slice, int mu, int color_index)
{
  START_CODE();
  
  if (color_index >= Nc || color_index < 0)
    QDP_error_exit("invalid color index", color_index);
  
  int spin_index = 0 ;

  // Write ONE to all field
  Real nnn = -1/4.0 ;
  Complex sitecomp = cmplx(nnn,0);
  ColorVector sitecolor = zero;
  StaggeredFermion sitefield = zero;
  
  pokeSpin(sitefield,
	   pokeColor(sitecolor,sitecomp,color_index),
	   spin_index);
  
  // Narrow the context to the desired slice.
  LatticeStaggeredFermion tmp;
  tmp = sitefield;  // QDP (not installed version) now supports   construct OLattice = OScalar
  
  // Auxiliary: Coordinates to use in "where" clauses
  multi1d<LatticeInteger> x(Nd);
  // Fill x with lattice coordinates
  for( int sigma = 0; sigma < Nd; sigma++) {
    x[ sigma ] = Layout::latticeCoordinate(sigma);
  }

  a = where(Layout::latticeCoordinate(mu) == slice && 
	    (x[0]+x[1]+x[2] ) % 2 ==  1 ,
	    tmp, LatticeStaggeredFermion(zero));
  
  END_CODE();
}

//
//  Putt all 1's in the main hypercube
//

void srcfil_local_o(LatticeStaggeredFermion& a, int slice, int mu, 
		    int color_index)
{

    Real one = -1;
    Complex sitecomp = cmplx(one,0);
    ColorVector sitecolor = zero;
    StaggeredFermion sitefield = zero;

    a = zero ; 

    const int spin_index = 0 ; 
    multi1d<int> coord(4) ; 

    coord[3] = slice ; 


    coord[0] = 1 ; coord[1] = 0 ; coord[2] = 0 ; 
    sitecolor = zero; sitefield = zero;
    pokeSite(a, pokeSpin(sitefield,
			 pokeColor(sitecolor,sitecomp,color_index),
			 spin_index), coord);

    coord[0] = 0 ; coord[1] = 1 ; coord[2] = 0 ; 
    sitecolor = zero; sitefield = zero;
    pokeSite(a, pokeSpin(sitefield,
			 pokeColor(sitecolor,sitecomp,color_index),
			 spin_index), coord);

    coord[0] = 0 ; coord[1] = 0 ; coord[2] = 1 ; 
    sitecolor = zero; sitefield = zero;
    pokeSite(a, pokeSpin(sitefield,
			 pokeColor(sitecolor,sitecomp,color_index),
			 spin_index), coord);

    coord[0] = 1 ; coord[1] = 1 ; coord[2] = 1 ; 
    sitecolor = zero; sitefield = zero;
    pokeSite(a, pokeSpin(sitefield,
			 pokeColor(sitecolor,sitecomp,color_index),
			 spin_index), coord);

}





//
//  Local even source
//

void srcfil_local_q(LatticeStaggeredFermion& a, int slice, int mu, 
		    int color_index)
{
    Real one = -1;
    Complex sitecomp = cmplx(one,0);

    a = zero ; 

    ColorVector sitecolor = zero;
    StaggeredFermion sitefield = zero;

    const int spin_index = 0 ; 
    multi1d<int> coord(4) ; 

    coord[3] = slice ; 

    coord[0] = 0 ; coord[1] = 0 ; coord[2] = 0 ; 
    sitecolor = zero; sitefield = zero;
    pokeSite(a, pokeSpin(sitefield,
			 pokeColor(sitecolor,sitecomp,color_index),
			 spin_index), coord);


    coord[0] = 1 ; coord[1] = 1 ; coord[2] = 0 ; 
    sitecolor = zero; sitefield = zero;
    pokeSite(a, pokeSpin(sitefield,
			 pokeColor(sitecolor,sitecomp,color_index),
			 spin_index), coord);


    coord[0] = 1 ; coord[1] = 0 ; coord[2] = 1 ; 
    sitecolor = zero; sitefield = zero;
    pokeSite(a, pokeSpin(sitefield,
			 pokeColor(sitecolor,sitecomp,color_index),
			 spin_index), coord);

    coord[0] = 0 ; coord[1] = 1 ; coord[2] = 1 ; 
    sitecolor = zero; sitefield = zero;
    pokeSite(a, pokeSpin(sitefield,
			 pokeColor(sitecolor,sitecomp,color_index),
			 spin_index), coord);

}


void create_stagg_source(LatticeStaggeredFermion & q_source,  
			 staggered_src_type src_type,
			 int color_source, int j_decay, int t_source)
{

  multi1d<int>  coord(4)   ;
  coord[0] = 0 ;
  coord[1] = 0 ;
  coord[2] = 0 ;
  coord[j_decay] = t_source ;

  q_source = zero ;
  if( src_type == local_s  )
    {
      //  Use a point source
      srcfil(q_source, coord, color_source);
    }    
  else if( src_type == wall_s  )
    {
      int  src_index = 0 ; 
      walfil(q_source, t_source, j_decay, color_source, src_index);
    }    
  else if( src_type ==  wall_q )
    {
      LatticeStaggeredFermion q_source_a , q_source_b  ;
      walfil_o(q_source_a, t_source, j_decay, color_source);
      walfil_q(q_source_b, t_source, j_decay, color_source);
      q_source = q_source_a + q_source_b ;
    }    
  else if( src_type ==  local_q )
    {
      LatticeStaggeredFermion q_source_a , q_source_b  ;
      srcfil_local_q(q_source_a, t_source, j_decay, color_source);
      srcfil_local_o(q_source_b, t_source, j_decay, color_source);
      q_source = q_source_a + q_source_b ;
    }    
  else
    {
      QDPIO::cerr << "Error,  src_type  " <<  src_type << " out of range" << endl;
      exit(0) ;
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

  if( input.param.gauge_trans )
    {
      QDPIO::cout << "DOING GAUGE TRANSFORM "  << endl;
      rgauge(u); 
      QDPIO::cout << "Random gauge transform on gauge fields done." << endl;


    }


  // Instantiate XML writer for the output file
  XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();
  push(xml_out, "hadron_corr");

  // Write out the input
  write(xml_out, "Input", xml_in);

  // Write out the config header
  write(xml_out, "Config_info", gauge_xml);

  // Check if the gauge field configuration is unitarized
  unitarityCheck(u);

  // Calculate some gauge invariant observables just for info.
  MesPlq(xml_out, "Gauge_Observables", u);
  xml_out.flush();

  // Fix to the coulomb gauge
  int n_gf;
  int j_decay = Nd-1;

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
      QDPIO::cerr << "Error reading first action name " << endl;
      throw;
    }


  XMLReader fermact_reader_A ;
  // Make a memory 'input stream' out of the XML, so we can open an
  // XML Reader on it.
  try{
  std::istringstream is(input.param.fermact_A.xml);

  // Open a reader on the memory stream.
  //  XMLReader fermact_reader(is);
  fermact_reader_A.open(is);
  }
  catch (...)
    {
      QDPIO::cerr << "Error reading first action name " << endl;
      throw;
    }


  Handle< StaggeredTypeFermAct< T,P,Q> > fermact(TheStagTypeFermActFactory::Instance().createObject(input.param.fermact.id, fermact_reader, input.param.fermact.path));

  Handle< StaggeredTypeFermAct< T,P,Q> > fermact_A(TheStagTypeFermActFactory::Instance().createObject(input.param.fermact_A.id, fermact_reader_A, input.param.fermact_A.path));

  // Cast of a pointer to a reference?
  StaggeredTypeFermAct<T,P,Q>& S_f   = *(fermact);
  StaggeredTypeFermAct<T,P,Q>& S_f_A = *(fermact_A);

  Real Mass = S_f.getQuarkMass() ;
  Real Mass_A = S_f_A.getQuarkMass() ;

  LatticeStaggeredPropagator quark_propagator;
  LatticeStaggeredPropagator quark_propagator_A;

  XMLBufferWriter xml_buf;
  int ncg_had = 0;
  int n_count;
  int t_length = input.param.nrow[3];

  LatticeStaggeredFermion q_source, psi;


  GroupXML_t inv_param;
  {
    XMLBufferWriter xml_buf;
    write(xml_buf, "InvertParam", input.param.invParam);
    XMLReader xml_in(xml_buf);
    inv_param = readXMLGroup(xml_in, "/InvertParam", "invType");
  }

  Handle< FermState<T,P,Q> > state(S_f.createState(u));
  Handle< SystemSolver<LatticeStaggeredFermion> > qprop(S_f.qprop(state,inv_param));

  Handle< FermState<T,P,Q> > state_A(S_f_A.createState(u));
  Handle< SystemSolver<LatticeStaggeredFermion> > qprop_A(S_f_A.qprop(state,inv_param));


  /** do inversions **************************/
  multi1d<int>  coord(4)   ;
  coord[0] = 0 ;
  coord[1] = 0 ;
  coord[2] = 0 ;
  coord[j_decay] = input.param.t_srce ;
  int t_source = input.param.t_srce ;
  QDPIO::cout << "Source time slice = " << input.param.t_srce << endl;

  psi = zero;   // note this is ``zero'' and not 0

  //
  // Loop over the source color, creating the source
  // and calling the relevant propagator routines.
  //

  push(xml_out, "first_inversion");

  push(xml_out,"Inverter_properties");
  write(xml_out, "Type_of_source",input.param.src_type_str);
  write(xml_out, "Mass" , Mass);
  write(xml_out, "RsdCG", input.param.invParam.RsdCG);
  write(xml_out, "SrcType", input.param.src_type_str);
  pop(xml_out);
  
  for(int color_source = 0; color_source < Nc; ++color_source) {
    QDPIO::cout << "Inversion[A] for Color =  " << color_source << endl;
    
    create_stagg_source(q_source,input.param.src_type,
			color_source, j_decay, t_source) ;

    // Use the last initial guess as the current guess

    // Compute the propagator for given source color/spin
    SystemSolverResults_t res = (*qprop)(psi, q_source);
    ncg_had += res.n_count;

    push(xml_out,"Inverter_performance");
    write(xml_out, "color", color_source);
    write(xml_out, "iterations", res.n_count);
    write(xml_out, "Final_RsdCG", res.resid);
    pop(xml_out);

    /*
     * Move the solution to the appropriate components
     * of quark propagator.
     */
    FermToProp(psi, quark_propagator, color_source);
  }  //color_source, first inversion

  pop(xml_out); //  inversion information

  push(xml_out, "second_inversion");

  push(xml_out,"Inverter_properties");
  write(xml_out, "Type_of_source", input.param.src_type_str);
  write(xml_out, "Mass" , Mass);
  write(xml_out, "RsdCG", input.param.invParam.RsdCG);
  write(xml_out, "SrcType", input.param.src_type_str);
  pop(xml_out);
  
  for(int color_source = 0; color_source < Nc; ++color_source) {
    QDPIO::cout << "Inversion[B] for Color =  " << color_source << endl;
    
    create_stagg_source(q_source,input.param.src_type,
			color_source, j_decay, t_source) ;


    /**DEBUG ** dump_text_src(q_source,input.param.nrow) ; exit(0) ; **/

    // Use the last initial guess as the current guess

    // Compute the propagator for given source color/spin
    SystemSolverResults_t res = (*qprop_A)(psi, q_source);
    ncg_had += res.n_count;

    push(xml_out,"Inverter_performance");
    write(xml_out, "color", color_source);
    write(xml_out, "iterations", res.n_count);
    write(xml_out, "Final_RsdCG", res.resid);
    pop(xml_out);

    /*
     * Move the solution to the appropriate components
     * of quark propagator.
     */
    FermToProp(psi, quark_propagator_A, color_source);
  }  //color_source, first inversion

  pop(xml_out); //  second inversion information

  int Bc_spec = 0; // perhaps 1 or -1




  {
    push(xml_out, "Hadron_A") ;
    write(xml_out, "source_time", t_source);
    write(xml_out, "Mass_A" , Mass_A);

    staggered_local_pion pion(t_length,u) ;
    pion.compute(quark_propagator_A,quark_propagator_A,j_decay) ;
    pion.dump(t_source,xml_out) ;

    multi1d<Complex>  barprop(t_length);

    // call the baryon_s function in baryon_s.cc
    baryon_s(quark_propagator_A,quark_propagator_A,
	     quark_propagator_A,
	     barprop,coord, j_decay,Bc_spec) ;
    write(xml_out, "baryon_corr", barprop) ; 

    pop(xml_out);  // Hadrons
  }


  {
    push(xml_out, "Hadron_B") ;
    write(xml_out, "source_time", t_source);
    write(xml_out, "Mass" , Mass);

    staggered_local_pion pion(t_length,u) ;
    pion.compute(quark_propagator,quark_propagator,j_decay) ;
    pion.dump(t_source,xml_out) ;

    multi1d<Complex>  barprop(t_length);

    // call the baryon_s function in baryon_s.cc
    baryon_s(quark_propagator,quark_propagator,quark_propagator,
	   barprop,coord, j_decay,Bc_spec) ;
    write(xml_out, "baryon_corr", barprop) ; 

    pop(xml_out);  // Hadrons
  }


  {
    push(xml_out, "Hadron_ABB") ;
    write(xml_out, "source_time", t_source);
    write(xml_out, "Mass_A" , Mass_A);
    write(xml_out, "Mass_B" , Mass);

    multi1d<Complex>  barprop(t_length);

    // call the baryon_s function in baryon_s.cc
    baryon_s(quark_propagator_A,quark_propagator,quark_propagator,
	   barprop,coord, j_decay,Bc_spec) ;
    write(xml_out, "baryon_corr", barprop) ; 

    pop(xml_out);  // Hadrons
  }


  {
    push(xml_out, "Hadron_AAB") ;
    write(xml_out, "source_time", t_source);
    write(xml_out, "Mass_A" , Mass_A);
    write(xml_out, "Mass_B" , Mass);

    multi1d<Complex>  barprop(t_length);

    // call the baryon_s function in baryon_s.cc
    baryon_s(quark_propagator_A,quark_propagator_A,
	     quark_propagator,barprop,coord, j_decay,Bc_spec) ;
    write(xml_out, "baryon_corr", barprop) ; 

    pop(xml_out);  // Hadrons
  }

  //
  // spin 3/2 baryons
  //

  bool do_spin3by_2 = false ; 
  if( input.param.src_type == local_q || 
      input.param.src_type == wall_q )
    { do_spin3by_2 = true ; } 



  // **************************************************
  //            Class 4 operator
  // **************************************************


  if( do_spin3by_2 ) 
  {
    push(xml_out, "Hadron_class4_AAB") ;
    write(xml_out, "source_time", t_source);
    write(xml_out, "Mass_A" , Mass_A);
    write(xml_out, "Mass_B" , Mass);

    multi1d<Complex>  barprop(t_length);

    // call the baryon_s function in baryon_s.cc
    baryon_class4_s(quark_propagator_A,
		    quark_propagator_A,
		    quark_propagator,
		    barprop,coord, j_decay,Bc_spec) ;
    write(xml_out, "baryon_corr", barprop) ; 

    pop(xml_out);  // Hadrons
  }



  if( do_spin3by_2 ) 
  {
    push(xml_out, "Hadron_class4_ABB") ;
    write(xml_out, "source_time", t_source);
    write(xml_out, "Mass_A" , Mass_A);
    write(xml_out, "Mass_B" , Mass);

    multi1d<Complex>  barprop(t_length);

    // call the baryon_s function in baryon_s.cc
    baryon_class4_s(quark_propagator_A,
		    quark_propagator,
		    quark_propagator,
		    barprop,coord, j_decay,Bc_spec) ;
    write(xml_out, "baryon_corr", barprop) ; 

    pop(xml_out);  // Hadrons
  }



  if( do_spin3by_2 ) 
  {
    push(xml_out, "Hadron_class4_BBB") ;
    write(xml_out, "source_time", t_source);
    write(xml_out, "Mass_B" , Mass);

    multi1d<Complex>  barprop(t_length);

    // call the baryon_s function in baryon_s.cc
    baryon_class4_s(quark_propagator,
		    quark_propagator,
		    quark_propagator,
		    barprop,coord, j_decay,Bc_spec) ;
    write(xml_out, "baryon_corr", barprop) ; 

    pop(xml_out);  // Hadrons
  }





  if( do_spin3by_2 ) 
  {
    push(xml_out, "Hadron_class4_AAA") ;
    write(xml_out, "source_time", t_source);
    write(xml_out, "Mass_A" , Mass_A);

    multi1d<Complex>  barprop(t_length);

    // call the baryon_s function in baryon_s.cc
    baryon_class4_s(quark_propagator_A,quark_propagator_A,
	     quark_propagator_A,barprop,coord, j_decay,Bc_spec) ;
    write(xml_out, "baryon_corr", barprop) ; 

    pop(xml_out);  // Hadrons
  }



  // **************************************************
  //            Class 7 operator
  // **************************************************


  if( do_spin3by_2 ) 
  {
    push(xml_out, "Hadron_class7_AAB") ;
    write(xml_out, "source_time", t_source);
    write(xml_out, "Mass_A" , Mass_A);
    write(xml_out, "Mass_B" , Mass);

    multi1d<Complex>  barprop(t_length);

    // call the baryon_s function in baryon_s.cc
    baryon_class7_s(quark_propagator_A,quark_propagator_A,
	     quark_propagator,barprop,coord, j_decay,Bc_spec) ;
    write(xml_out, "baryon_corr", barprop) ; 

    pop(xml_out);  // Hadrons
  }


  if( do_spin3by_2 ) 
  {
    push(xml_out, "Hadron_class7_ABB") ;
    write(xml_out, "source_time", t_source);
    write(xml_out, "Mass_A" , Mass_A);
    write(xml_out, "Mass_B" , Mass);

    multi1d<Complex>  barprop(t_length);

    // call the baryon_s function in baryon_s.cc
    baryon_class7_s(quark_propagator_A,quark_propagator,
	     quark_propagator,barprop,coord, j_decay,Bc_spec) ;
    write(xml_out, "baryon_corr", barprop) ; 

    pop(xml_out);  // Hadrons
  }




  if( do_spin3by_2 ) 
  {
    push(xml_out, "Hadron_class7_AAA") ;
    write(xml_out, "source_time", t_source);
    write(xml_out, "Mass_A" , Mass_A);

    multi1d<Complex>  barprop(t_length);

    // call the baryon_s function in baryon_s.cc
    baryon_class7_s(quark_propagator_A,
		    quark_propagator_A,
		    quark_propagator_A,
		    barprop,coord, j_decay,Bc_spec) ;
    write(xml_out, "baryon_corr", barprop) ; 

    pop(xml_out);  // Hadrons
  }




  if( do_spin3by_2 ) 
  {
    push(xml_out, "Hadron_class7_BBB") ;
    write(xml_out, "source_time", t_source);
    write(xml_out, "Mass_B" , Mass);

    multi1d<Complex>  barprop(t_length);

    // call the baryon_s function in baryon_s.cc
    baryon_class7_s(quark_propagator,
		    quark_propagator,
		    quark_propagator,
		    barprop,coord, j_decay,Bc_spec) ;
    write(xml_out, "baryon_corr", barprop) ; 

    pop(xml_out);  // Hadrons
  }





  // **************************************************
  //            Class 7 NLT operator
  // **************************************************


  if( do_spin3by_2 ) 
  {
    push(xml_out, "Hadron_class7_NLT_AAB") ;
    write(xml_out, "source_time", t_source);
    write(xml_out, "Mass_A" , Mass_A);
    write(xml_out, "Mass_B" , Mass);

    multi1d<Complex>  barprop(t_length);

    // call the baryon_s function in baryon_s.cc
    baryon_class7_NLT_s(quark_propagator_A,quark_propagator_A,
	     quark_propagator,u,barprop,coord, j_decay,Bc_spec) ;
    write(xml_out, "baryon_corr", barprop) ; 

    pop(xml_out);  // Hadrons
  }


  if( do_spin3by_2 ) 
  {
    push(xml_out, "Hadron_class7_NLT_ABB") ;
    write(xml_out, "source_time", t_source);
    write(xml_out, "Mass_A" , Mass_A);
    write(xml_out, "Mass_B" , Mass);

    multi1d<Complex>  barprop(t_length);

    // call the baryon_s function in baryon_s.cc
    baryon_class7_NLT_s(quark_propagator_A,quark_propagator,
	     quark_propagator,u,barprop,coord, j_decay,Bc_spec) ;
    write(xml_out, "baryon_corr", barprop) ; 

    pop(xml_out);  // Hadrons
  }




  if( do_spin3by_2 ) 
  {
    push(xml_out, "Hadron_class7_NLT_AAA") ;
    write(xml_out, "source_time", t_source);
    write(xml_out, "Mass_A" , Mass_A);

    multi1d<Complex>  barprop(t_length);

    // call the baryon_s function in baryon_s.cc
    baryon_class7_NLT_s(quark_propagator_A,
		    quark_propagator_A,
			quark_propagator_A,u,
		    barprop,coord, j_decay,Bc_spec) ;
    write(xml_out, "baryon_corr", barprop) ; 

    pop(xml_out);  // Hadrons
  }




  if( do_spin3by_2 ) 
  {
    push(xml_out, "Hadron_class7_NLT_BBB") ;
    write(xml_out, "source_time", t_source);
    write(xml_out, "Mass_B" , Mass);

    multi1d<Complex>  barprop(t_length);

    // call the baryon_s function in baryon_s.cc
    baryon_class7_NLT_s(quark_propagator,
			quark_propagator,
			quark_propagator,u,
			barprop,coord, j_decay,Bc_spec) ;
    write(xml_out, "baryon_corr", barprop) ; 

    pop(xml_out);  // Hadrons
  }





  // ***** ------- ********** ----- ****** ----------------
  pop(xml_out);

  xml_in.close();
  xml_out.close();

  // Time to bolt
  Chroma::finalize();

  exit(0);
}
