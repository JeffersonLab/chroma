/*
 *  This is the top-level routine for Stout-link smearing.
 *  Reading in gauge fields and writing out stout-smeared gauge fields
 */

#include <iostream>
#include <cstdio>

#include "chroma.h"

#include <sys/time.h>   // for timings

using namespace Chroma;

struct Param_t
{
  int link_smear_num;
  Real link_smear_fact;		// Smearing parameters

  int j_decay;			// Decay direction
  multi1d<int> nrow;		// Lattice dimension
  
};

struct Stout_t
{
  QDP_volfmt_t  volfmt;
  CfgType	stout_type;
  string 	stout_file;	// storage for stout config
};

struct Stoutsmear_input_t
{
  Param_t          param;
  Cfg_t            cfg;
  Stout_t	   stout;
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
     QDP_error_exit( "stoutsmear.cc: invalid number of stout smearing iterations, link_smear_num = %d", param.link_smear_num );
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
void read(XMLReader& xml, const string& path, Stout_t& input)
{
  XMLReader inputtop(xml, path);
  read(inputtop, "volfmt", input.volfmt);
  read(inputtop, "stout_type", input.stout_type);
  read(inputtop, "stout_file", input.stout_file);
}

// Reader for input parameters
void read(XMLReader& xml, const string& path, Stoutsmear_input_t& input)
{
  XMLReader inputtop(xml, path);

  // Read all the input groups
  try
  {
    // Read program parameters
    read(inputtop, "Param", input.param);

    // Read in the gauge configuration info
    read(inputtop, "Cfg", input.cfg);

    // Read in the stout outfile
    read(inputtop, "Stout", input.stout);
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading stoutsmear data: " << e << endl;
    throw;
  }
}


int main(int argc, char *argv[])
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  START_CODE();

  // Parameter structure for the input
  Stoutsmear_input_t input;

  // Instantiate xml reader for DATA
  XMLReader xml_in("DATA");

  // Read data
  read(xml_in, "/stoutsmear", input);

  Layout::setLattSize(input.param.nrow);
  Layout::create();

  QDPIO::cout << " STOUTSMEAR: STOUT smearing of gauge config" << endl;

  // Read in the configuration along with relevant information.
  multi1d<LatticeColorMatrix> u(Nd);
  XMLReader gauge_file_xml, gauge_xml;

  SzinGauge_t  szin_gauge_header;

  // Startup gauge
  gaugeStartup(gauge_file_xml, gauge_xml, u, input.cfg);

  read(gauge_xml, "/szin", szin_gauge_header);

  // Instantiate XML writer
  XMLFileWriter xml_out("XMLDAT");
  push(xml_out, "stoutsmear");
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
  MesPlq(xml_out, "Observables", u);
  xml_out.flush();

  // Now stout smear
  multi1d<LatticeColorMatrix> u_stout(Nd);
  u_stout = u;

  if (input.param.link_smear_num > 0)
  {

    QDPIO::cout << "STOUT Smear gauge field" << endl;

    for(int i=0; i < input.param.link_smear_num; ++i)
    {
      multi1d<LatticeColorMatrix> u_tmp(Nd);

      for(int mu = 0; mu < Nd; ++mu)
        if ( mu != input.param.j_decay )
          stout_smear(u_tmp[mu], u_stout, mu,
                      input.param.link_smear_fact,input.param.j_decay);
        else
          u_tmp[mu] = u_stout[mu];

      u_stout = u_tmp;
    }
    QDPIO::cout << "Gauge field STOUT-smeared!" << endl;
    
    // Write out what is done
    push(xml_out, "Gauge field STOUT-smeared!");
    push(xml_out,"Smearing_parameters");
    write(xml_out, "link_smear_num",input.param.link_smear_num);
    write(xml_out, "link_smear_fact",input.param.link_smear_fact);
    write(xml_out, "j_decay",input.param.j_decay);
    pop(xml_out);
    xml_out.flush();
  }

  // Check if the smeared gauge field is unitary
  unitarityCheck(u_stout);
  
  // Again calculate some gauge invariant observables
  MesPlq(xml_out, "STOUT_observables", u_stout);
  xml_out.flush();

  // Now write the configuration to disk
  XMLBufferWriter gauge_file_xml_out, gauge_xml_out;
  push(gauge_file_xml_out, "gauge");
  write(gauge_file_xml_out, "id", int(0));
  pop(gauge_file_xml_out);
  write(gauge_xml_out, "szin", szin_gauge_header);
  writeGauge(gauge_file_xml_out, gauge_xml_out, u_stout,
             input.stout.stout_file,input.stout.volfmt, QDPIO_SERIAL);
  // writeSzin(szin_gauge_header, u_stout, input.stout.stout_file);

  END_CODE();

  // Time to bolt
  QDP_finalize();

  exit(0);
}
