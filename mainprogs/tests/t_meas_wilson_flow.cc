// $Id: wilson_flow_w.cc,v 1.12 2011/12/11 17:23:22 cmcneile Exp cmcneile $
/*! \file
 *  \brief Wrapper code to the measure the Wilson flow
 *   
 *
 *
 */

#include <iostream>
#include <cstdio>

#define MAIN

// Include everything...
#include "chroma.h"

#include "meas/glue/wilson_flow_w.h"

using namespace Chroma;


 bool linkageHack(void)
{
  bool foo = true;

  // Inline Measurements
  foo &= InlineAggregateEnv::registerAll();
  foo &= GaugeInitEnv::registerAll();

  return foo;
}





/*
 * Input 
 */


// Parameters which must be determined from the XML input
// and written to the XML output
struct Param_t
{

  Real GFAccu, OrPara;    // Gauge fixing tolerance and over-relaxation param
  int GFMax;              // Maximum gauge fixing iterations

  multi1d<int> nrow;
  multi1d<int> boundary;
  int nstep ; 
  bool do_gauge_transform ;

  Real wflow_eps ; 

};



struct Propagator_input_t
{
  IO_version_t     io_version;
  Param_t          param;
  Cfg_t            cfg;
};


//
// Reader for input parameters
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
    QDPIO::cerr << "Error reading IO version: " << e << endl;
    throw;
  }

   XMLReader paramtop(inputtop, "param"); // push into 'param' group

   try 
     {

       read(paramtop, "nrow", input.param.nrow);
       read(paramtop, "wflow_eps", input.param.wflow_eps);
       read(paramtop, "nstep", input.param.nstep);
       read(paramtop, "gauge_trans", input.param.do_gauge_transform);

     }
   catch (const string& e) 
     {
       QDPIO::cerr << "Error parameters: " << e << endl;
       throw;
     }


  //
  //   outside <param>  </param>
  //


  // Read in the gauge configuration file name
  try
  {
    read(inputtop, "Cfg", input.cfg);
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error gauge name data: " << e << endl;
    throw;
  }
}




//! Measure Wilson flow 
/*! \defgroup calculation of Wilson flow 
 *  \ingroup 
 *
 * Main program for calculation of Wilson flow
 */

int main(int argc, char **argv)
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);

  QDPIO::cout << "Measure the Wilson flow " << endl;
  QDPIO::cout << "Calculation for SU(" << Nc << ")" << endl;
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


  // Read data
  read(xml_in, "/WilsonFlow", input);


  // Specify lattice size, shape, etc.
  Layout::setLattSize(input.param.nrow);
  Layout::create();

  // Read in the configuration along with relevant information.
  multi1d<LatticeColorMatrix> u(Nd);



  XMLReader gauge_file_xml, gauge_xml;
  // Start up the gauge field
  gaugeStartup(gauge_file_xml, gauge_xml, u, input.cfg);


  // Check if the gauge field configuration is unitarized
  unitarityCheck(u);

  // Instantiate XML writer for XMLDAT
  XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();
  push(xml_out, "wilsonFlow");

  // Write out the input
  write(xml_out, "Input", xml_in);


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

  if( input.param.do_gauge_transform )
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
  else
    {
       QDPIO::cout << "NO RANDOM GAUGE TRANSFORM" << endl;
    }



  int j_decay = Nd-1;
  int jomit=3 ; 
  wilson_flow(xml_out, u, input.param.nstep,input.param.wflow_eps ,jomit) ;


  pop(xml_out);

  xml_out.close();
  xml_in.close();


  // Time to bolt
  QDPIO::cout << "End of measurements " << endl;
  Chroma::finalize();

  exit(0);
}
