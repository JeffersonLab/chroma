/*
 *  $Id: hypsmear4d.cc,v 1.2 2005-02-21 19:28:59 edwards Exp $
 *
 *  This is the top-level routine for HYP smearing.
 *  Reading in gauge fields and writing out hyp-smear gauge fields
 *  in a generic format. The program hypsmear.cc allows truncating
 *  the config size, but only supports szin format. This program
 *  writes only SCIDAC format.
 */


#include <iostream>
#include <cstdio>

#include "chroma.h"

#include <sys/time.h>   // for timings

using namespace Chroma;

struct Param_t
{
  Real alpha1;			// Smearing parameters
  Real alpha2;
  Real alpha3;

  int link_smear_num;           // Number of smearing iterations

  multi1d<int> nrow;		// Lattice dimension

};

struct Hyp_t
{
  QDP_volfmt_t  volfmt;	   // single or multi file volume format
  CfgType  hyp_type;       // storage type for hyp config
  string   hyp_file;	   // storage for hyp config
};

struct Hypsmear_input_t
{
  Param_t          param;
  Cfg_t            cfg;
  Hyp_t            hyp;
};



// Reader for out gauge file
void read(XMLReader& xml, const string& path, Hyp_t& input)
{
  XMLReader inputtop(xml, path);
  read(inputtop, "volfmt", input.volfmt);
  read(inputtop, "hyp_type", input.hyp_type);
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

    read(paramtop, "link_smear_num", param.link_smear_num);

    if( param.link_smear_num < 0 )
    {
      QDP_error_exit( "hypsmear.cc: invalid number of hyp smearing iterations, link_smear_num = %d", param.link_smear_num );
    }
    break;

  case 2:

    param.link_smear_num = 1;
    break;

  default :

    QDPIO::cerr << "Input version " << version << " unsupported." << endl;
    QDP_abort(1);
  }

  read(paramtop, "alpha1", param.alpha1);
  read(paramtop, "alpha2", param.alpha2);
  read(paramtop, "alpha3", param.alpha3);

  read(paramtop, "nrow", param.nrow);
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

    // Read in the hyp outfile info
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

  SzinGauge_t  szin_gauge_header;

  // Startup gauge
  gaugeStartup(gauge_file_xml, gauge_xml, u, input.cfg);

  read(gauge_xml, "/szin", szin_gauge_header);

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
  unitarityCheck(u);

  // Calculate some gauge invariant observables
  Double w_plaq, s_plaq, t_plaq, link;
  MesPlq(u, w_plaq, s_plaq, t_plaq, link);

  push(xml_out, "Observables");
  write(xml_out, "w_plaq",w_plaq);
  write(xml_out, "s_plaq", s_plaq);
  write(xml_out, "t_plaq", t_plaq);
  write(xml_out, "link", link);

  pop(xml_out);
  xml_out.flush();


  // Now hyp smear
  multi1d<LatticeColorMatrix> u_hyp(Nd);
  u_hyp = u;

  Real BlkAccu = 1.0e-5;
  int BlkMax = 100;

  if (input.param.link_smear_num > 0)
  {

    QDPIO::cout << "HYP Smear gauge field" << endl;

    for(int i=0; i < input.param.link_smear_num; i++)
    {
    Hyp_Smear(u, u_hyp, 
	      input.param.alpha1, input.param.alpha2, input.param.alpha3, 
	      BlkAccu, BlkMax);
    u = u_hyp;
    }
    QDPIO::cout << "Gauge field HYP-smeared!" << endl;

    // Write out what is done
    push(xml_out, "Gauge field HYP-smeared!");
    push(xml_out,"Smearing_parameters");
    write(xml_out, "link_smear_num",input.param.link_smear_num);
    write(xml_out, "alpha1",input.param.alpha1);
    write(xml_out, "alpha2",input.param.alpha2);
    write(xml_out, "alpha3",input.param.alpha3);
    pop(xml_out);
    xml_out.flush();
  }

  // Again calculate some gauge invariant observables
  MesPlq(u_hyp, w_plaq, s_plaq, t_plaq, link);

  push(xml_out, "HYP_observables");
  write(xml_out, "w_plaq", w_plaq);
  write(xml_out, "s_plaq", s_plaq);
  write(xml_out, "t_plaq", t_plaq);
  write(xml_out, "link", link);

  pop(xml_out);
  xml_out.flush();

  // Now write the configuration to disk
  XMLBufferWriter gauge_file_xml_out, gauge_xml_out;
  push(gauge_file_xml_out, "gauge");
  write(gauge_file_xml_out, "id", int(0));
  pop(gauge_file_xml_out);
  write(gauge_xml_out, "szin", szin_gauge_header);
  writeGauge(gauge_file_xml_out, gauge_xml_out, u_hyp,
             input.hyp.hyp_file,input.hyp.volfmt, QDPIO_SERIAL);
  // writeSzin(szin_gauge_header, u_hyp, input.hyp.hyp_file);

  END_CODE();

  // Time to bolt
  ChromaFinalize();

  exit(0);
}
