/*
 *  This is the top-level routine for APE smearing.
 *  Reading in gauge fields and writing out ape-smeared gauge fields
 */

#include <iostream>
#include <cstdio>

#include "chroma.h"

using namespace Chroma;


struct Param_t
{
  int link_smear_num;
  Real link_smear_fact;		// Smearing parameters

  int j_decay;			// Decay direction
  multi1d<int> nrow;		// Lattice dimension
  
};

struct Ape_t
{
  QDP_volfmt_t  volfmt;
  CfgType	ape_type;
  string 	ape_file;	// storage for ape config
};

struct Apesmear_input_t
{
  Param_t          param;
  Cfg_t            cfg;
  Ape_t	           ape;
};

//! Parameters for running code
void read(XMLReader& xml, const string& path, Param_t& param)
{
  XMLReader paramtop(xml, path);

  int version;
  read(paramtop, "version", version);

  switch (version) 
  {
  case 2:

    read(paramtop, "link_smear_num", param.link_smear_num);

    if( param.link_smear_num < 0 )
    {
     QDP_error_exit( "apesmear.cc: invalid number of ape smearing iterations, link_smear_num = %d", param.link_smear_num );
    }

    break;

  default :

    QDPIO::cerr << "Input version " << version << " unsupported." << endl;
    QDP_abort(1);
  }

  read(paramtop, "link_smear_fact", param.link_smear_fact);
  read(paramtop, "j_decay", param.j_decay);
  read(paramtop, "nrow", param.nrow);

}

// Reader for out gauge file
void read(XMLReader& xml, const string& path, Ape_t& input)
{
  XMLReader inputtop(xml, path);
  read(inputtop, "volfmt", input.volfmt);
  read(inputtop, "ape_type", input.ape_type);
  read(inputtop, "ape_file", input.ape_file);
}

// Reader for input parameters
void read(XMLReader& xml, const string& path, Apesmear_input_t& input)
{
  XMLReader inputtop(xml, path);

  // Read all the input groups
  try
  {
    // Read program parameters
    read(inputtop, "Param", input.param);

    // Read in the gauge configuration info
    read(inputtop, "Cfg", input.cfg);

    // Read in the ape outfile
    read(inputtop, "Ape", input.ape);
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading apesmear data: " << e << endl;
    throw;
  }
}


int main(int argc, char *argv[])
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  START_CODE();

  // Parameter structure for the input
  Apesmear_input_t input;

  // Instantiate xml reader for DATA
  XMLReader xml_in("DATA");

  // Read data
  read(xml_in, "/apesmear", input);

  Layout::setLattSize(input.param.nrow);
  Layout::create();

  QDPIO::cout << " APESMEAR: APE smearing of gauge config" << endl;

  // Read in the configuration along with relevant information.
  multi1d<LatticeColorMatrix> u(Nd);
  XMLReader gauge_file_xml, gauge_xml;

  SzinGauge_t  szin_gauge_header;
  initHeader(szin_gauge_header);


  // Startup gauge
  gaugeStartup(gauge_file_xml, gauge_xml, u, input.cfg);

  read(gauge_xml, "/szin", szin_gauge_header);

  // Instantiate XML writer
  XMLFileWriter xml_out("XMLDAT");
  push(xml_out, "apesmear");
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
  // unitarityCheck(u);

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


  // Now ape smear
  multi1d<LatticeColorMatrix> u_ape(Nd);
  u_ape = u;

  if (input.param.link_smear_num > 0)
  {

    QDPIO::cout << "APE Smear gauge field" << endl;

    int BlkMax = 100;
    Real BlkAccu = 1.0e-5;

    for(int i=0; i < input.param.link_smear_num; ++i)
    {
      multi1d<LatticeColorMatrix> u_tmp(Nd);

      for(int mu = 0; mu < Nd; ++mu)
        if ( mu != input.param.j_decay )
          APE_Smear(u_ape, u_tmp[mu], mu, 0,
                    input.param.link_smear_fact, BlkAccu, BlkMax,
                    input.param.j_decay);
        else
          u_tmp[mu] = u_ape[mu];

      u_ape = u_tmp;
    }
    QDPIO::cout << "Gauge field APE-smeared!" << endl;
  }

  // Write out what is done
  push(xml_out, "Gauge field APE-smeared!");
  push(xml_out,"Smearing_parameters");
  write(xml_out, "link_smear_num",input.param.link_smear_num);
  write(xml_out, "link_smear_fact",input.param.link_smear_fact);
  write(xml_out, "j_decay",input.param.j_decay);
  pop(xml_out);
  xml_out.flush();
  
   // Check if the smeared gauge field is unitary
  unitarityCheck(u_ape);
  
  // Again calculate some gauge invariant observables
  MesPlq(u_ape, w_plaq, s_plaq, t_plaq, link);

  push(xml_out, "APE_observables");
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
  writeGauge(gauge_file_xml_out, gauge_xml_out, u_ape, input.ape.ape_file,
             input.ape.volfmt, QDPIO_SERIAL);
  // writeSzin(szin_gauge_header, u_ape, input.ape.ape_file);

  END_CODE();

  // Time to bolt
  QDP_finalize();

  exit(0);
}
