/*
 *  $Id: hypsmear.cc,v 1.20 2005-02-21 19:28:59 edwards Exp $
 *
 *  This is the top-level routine for HYP smearing.
 *  It is a wrapper for Urs' and Robert's implmenetation of the HYP
 *  smearing algorithm
 */


#include <iostream>
#include <cstdio>

#include "chroma.h"

#include <sys/time.h>   // for timings

using namespace Chroma;

/*
 *  This is the reader for the input parameters
 */

/*
 * Input 
 */

// Parameters which must be determined from the XML input
// and written to the XML output
struct Param_t
{
  Real alpha1;			// Smearing parameters
  Real alpha2;
  Real alpha3;

  int num_smear;                // Number of smearing iterations

  multi1d<int> nrow;		// Lattice dimension
  int j_decay;			// Direction corresponding to time

  /*
   *  Now some various rules for truncating the configuration
   */
  int trunc;			// Whether to truncate the output
  int t_start;			// Starting time slice
  int t_end;			// Ending time slice
};

struct Hyp_t
{
  CfgType  cfg_type;       // storage order for stored gauge configuration
  string   hyp_file;
};

struct Hypsmear_input_t
{
  Param_t          param;
  Cfg_t            cfg;
  Hyp_t            hyp;
};



//! Target file
void read(XMLReader& xml, const string& path, Hyp_t& input)
{
  XMLReader inputtop(xml, path);

  read(inputtop, "cfg_type", input.cfg_type);
  read(inputtop, "hyp_file", input.hyp_file);
}


//! Parameters for running code
void read(XMLReader& xml, const string& path, Param_t& param)
{
  XMLReader paramtop(xml, path);

  int version;
  read(paramtop, "version", version);

  switch (version) 
  {
  case 3:
    /**************************************************************************/

    /* this version allows a variable num_smear */
    read(paramtop, "num_smear", param.num_smear);

    if( param.num_smear < 0 )
    {
      QDP_error_exit( "hypsmear.cc: invalid number of hyp smearing iterations, num_smear = %d", param.num_smear );
    }

    break;

  case 2:
    /**************************************************************************/

    /* this version only allows num_smear = 1 */
    param.num_smear = 1;

    break;

  default :
    /**************************************************************************/

    QDPIO::cerr << "Input parameter version " << version << " unsupported." << endl;
    QDP_abort(1);
  }

  read(paramtop, "alpha1", param.alpha1);
  read(paramtop, "alpha2", param.alpha2);
  read(paramtop, "alpha3", param.alpha3);

  read(paramtop, "nrow", param.nrow);

  /*
   *  Now information about whether to truncate the configuration
   */
  read(paramtop, "trunc", param.trunc);
  switch(param.trunc){
  case 1:
    read(paramtop, "t_start", param.t_start);
    read(paramtop, "t_end", param.t_end);
    read(paramtop, "j_decay", param.j_decay);
    break;
  default:
    break;
  }
}


// Reader for input parameters
void read(XMLReader& xml, const string& path, Hypsmear_input_t& input)
{
  XMLReader inputtop(xml, path);

  // Read all the input groups
  try
  {
    // Read program parameters
    read(inputtop, "Param", input.param);

    // Read in the gauge configuration info
    read(inputtop, "Cfg", input.cfg);

    // Read in the hyp file info
    read(inputtop, "Hyp", input.hyp);
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading hypsmear data: " << e << endl;
    throw;
  }
}



int main(int argc, char *argv[])
{
  // Put the machine into a known state
  ChromaInitialize(&argc, &argv);

  START_CODE();

  // Parameter structure for the input
  Hypsmear_input_t input;

  // Instantiate xml reader for DATA
  XMLReader xml_in("./DATA");

  // Read data
  read(xml_in, "/hypsmear", input);

  Layout::setLattSize(input.param.nrow);
  Layout::create();

  QDPIO::cout << " HYPSMEAR: HYP smearing of gauge config" << endl;

  // Read in the configuration along with relevant information.
  multi1d<LatticeColorMatrix> u(Nd);
  XMLReader gauge_file_xml, gauge_xml;

  clock_t t1 = clock();

  // Startup gauge
  gaugeStartup(gauge_file_xml, gauge_xml, u, input.cfg);

  clock_t t2 = clock();
  QDPIO::cout << "Gauge read took " << (double)((int)(t2)-(int)(t1))/(double)(CLOCKS_PER_SEC) << " secs" << endl;


  // Instantiate XML writer for XMLDAT
  XMLFileWriter& xml_out = TheXMLOutputWriter::Instance();
  push(xml_out, "hypsmear");

  proginfo(xml_out);    // Print out basic program info

  // Write out the input
  write(xml_out, "Input", xml_in);

  // Write out the config header
  write(xml_out, "Config_info", gauge_xml);

  push(xml_out, "Output_version");
  write(xml_out, "out_version", 1);
  pop(xml_out);

  xml_out.flush();


  // Check if the gauge field configuration is unitarized
  t1 = clock();
  unitarityCheck(u);
  t2 = clock();
  QDPIO::cout << "Unitarity took " << (double)((int)(t2)-(int)(t1))/(double)(CLOCKS_PER_SEC) << " secs" << endl;
  

  // Calculate some gauge invariant observables just for info.
  Double w_plaq, s_plaq, t_plaq, link;
  t1 = clock();
  MesPlq(u, w_plaq, s_plaq, t_plaq, link);
  t2 = clock();
  QDPIO::cout << "Plaquette took " << (double)((int)(t2)-(int)(t1))/(double)(CLOCKS_PER_SEC) << " secs" << endl;

  push(xml_out, "Observables");
  write(xml_out, "w_plaq",w_plaq);
  write(xml_out, "s_plaq", s_plaq);
  write(xml_out, "t_plaq", t_plaq);
  write(xml_out, "link", link);
  pop(xml_out);

  xml_out.flush();


  // Now hyp smear
  multi1d<LatticeColorMatrix> u_hyp(Nd);

  Real BlkAccu = 1.0e-5;
  int BlkMax = 100;


  t1 = clock();
  for( int n = 0; n < input.param.num_smear; n ++ )
  {
    Hyp_Smear(u, u_hyp, 
	      input.param.alpha1, input.param.alpha2, input.param.alpha3, 
	      BlkAccu, BlkMax);
    u = u_hyp;
  }
  if( input.param.num_smear == 0 )
  {
    u_hyp = u;
  }
  t2 = clock();
  QDPIO::cout << "Hypsmear took " << (double)((int)(t2)-(int)(t1))/(double)(CLOCKS_PER_SEC) << " secs" << endl;

  // Calculate some gauge invariant observables just for info.

  MesPlq(u_hyp, w_plaq, s_plaq, t_plaq, link);

  push(xml_out, "HYP_observables");
  write(xml_out, "w_plaq", w_plaq);
  write(xml_out, "s_plaq", s_plaq);
  write(xml_out, "t_plaq", t_plaq);
  write(xml_out, "link", link);
  pop(xml_out);

  // Now write the configuration to disk

  QDPIO::cout << "call szin trunc-er" << endl;

  t1 = clock();

  switch (input.hyp.cfg_type) 
  {
  case CFG_TYPE_SZIN :
  {
    SzinGauge_t szin_out;

    switch(input.param.trunc)
    {
    case 1:
      QDPIO::cout << "Call writeSzinTrunc" << endl;
      writeSzinTrunc(szin_out, u_hyp, input.param.j_decay,
		     input.param.t_start, input.param.t_end,
		     input.hyp.hyp_file);
      break;
    default:
      QDPIO::cout << "Call writeSzin" << endl;
      writeSzin(szin_out, u_hyp,
		input.hyp.hyp_file);
      break;
    }
    break;
  }
  default :
    QDP_error_exit("Configuration type is unsupported.");
  }

  t2 = clock();
  QDPIO::cout << "Gauge write took " << (double)((int)(t2)-(int)(t1))/(double)(CLOCKS_PER_SEC) << " secs" << endl;

  pop(xml_out);

  END_CODE();

  // Time to bolt
  ChromaFinalize();

  exit(0);
}
