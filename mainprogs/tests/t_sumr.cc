// $Id: t_sumr.cc,v 1.5 2005-01-14 20:13:10 edwards Exp $

#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>

#include <cstdio>

#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

#include "chroma.h"

using namespace Chroma;

struct App_input_t {
  ChromaProp_t param;
  Cfg_t        cfg;
};

// Reader for input parameters
void read(XMLReader& xml, const string& path, App_input_t& input)
{
  XMLReader inputtop(xml, path);

  // Read the input
  try
  {
    // The parameters holds the version number
    read(inputtop, "Param", input.param);

    // Read in the gauge configuration info
    read(inputtop, "Cfg", input.cfg);
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading data: " << e << endl;
    throw;
  }
}


int main(int argc, char **argv)
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);



  App_input_t input;
  XMLReader xml_in("DATA");

  try {
    read(xml_in, "/ovlapTest", input);
  }
   catch( const string& e) { 
    QDPIO::cerr << "Caught Exception : " << e << endl;
    QDP_abort(1);
  }


  // Setup the lattice
  Layout::setLattSize(input.param.nrow);
  Layout::create();

  multi1d<LatticeColorMatrix> u(Nd);
  XMLReader gauge_file_xml, gauge_xml;
  gaugeStartup(gauge_file_xml, gauge_xml, u, input.cfg);

  XMLFileWriter xml_out("XMLDAT");
  push(xml_out,"t_g5eps_bj");


  // Measure the plaquette on the gauge
  Double w_plaq, s_plaq, t_plaq, link;
  MesPlq(u, w_plaq, s_plaq, t_plaq, link);
  push(xml_out, "plaquette");
  write(xml_out, "w_plaq", w_plaq);
  write(xml_out, "s_plaq", s_plaq);
  write(xml_out, "t_plaq", t_plaq);
  write(xml_out, "link", link);
  pop(xml_out);

  // Create a FermBC
  Handle<FermBC<LatticeFermion> >  fbc(new SimpleFermBC<LatticeFermion>(input.param.boundary));


  QDPIO::cout << "FERM_ACT_ZOLOTAREV_4D" << endl;
  const Zolotarev4DFermActParams& zolo4d = dynamic_cast<const Zolotarev4DFermActParams& > (*(input.param.FermActHandle));
      
  // Construct Fermact -- now uses constructor from the zolo4d params
  // struct
  Zolotarev4DFermAct S(fbc, zolo4d, xml_out);

  Handle<const ConnectState> connect_state(S.createState(u, zolo4d.StateInfo, xml_out,zolo4d.AuxFermActHandle->getMass()));


  int G5 = Ns*Ns - 1;
  LatticeFermion psi,chi;
  int n_count;


  // Solve on chiral point source
  multi1d<int> coord(4);
  coord[0] = 0;
  coord[1] = 0;
  coord[2] = 0;
  coord[3] = 0;
  QDP::StopWatch swatch;

  chi=zero;
  srcfil(chi, coord, 0, 0);
  
  psi = zero;
  swatch.reset();
  swatch.start();
  S.qprop(psi,
	  connect_state,
	  chi,
	  SUMR_INVERTER,
	  input.param.invParam.RsdCG,
	  input.param.invParam.MaxCG,
	  n_count);

  swatch.stop();

  double t = swatch.getTimeInSeconds();

  QDPIO::cout << "CG on point source took: " << n_count << " iterations" << endl;
  QDPIO::cout << "Wall clock time : " << t << " seconds" << endl;
  push(xml_out, "CGPointSource");
  write(xml_out, "n_count", n_count);
  write(xml_out, "t",       t);
  pop(xml_out);

  psi = zero;
  swatch.reset();
  swatch.start();
  S.qprop(psi,
	  connect_state,
	  chi,
	  REL_SUMR_INVERTER,
	  input.param.invParam.RsdCG,
	  input.param.invParam.MaxCG,
	  n_count);
  swatch.stop();
  t = swatch.getTimeInSeconds();
  QDPIO::cout << "SUMR on point source took: " << n_count << " iterations" << endl;
  QDPIO::cout << "Wall clock time : " << t << " seconds" << endl;
  push(xml_out, "SUMRPointSource");
  write(xml_out, "n_count", n_count);
  write(xml_out, "t",       t);
  pop(xml_out);

#if 0
  // Solve on non chiral sources
  gaussian(chi);
  chi /= sqrt(norm2(chi));

  psi = zero;
  swatch.reset();
  swatch.start();
  S.qprop(psi,
	  connect_state,
	  chi,
	  CG_INVERTER,
	  input.param.invParam.RsdCG,
	  input.param.invParam.MaxCG,
	  n_count);
  swatch.stop();
  t = swatch.getTimeInSeconds();

  QDPIO::cout << "CG on gaussian source took: " << n_count << " iterations" << endl;
  QDPIO::cout << "Wall clock time : " << t << " seconds" << endl;
  push(xml_out, "CGGaussianSource");
  write(xml_out, "n_count", n_count);
  write(xml_out, "t",       t);
  pop(xml_out);

  psi = zero;
  swatch.reset();
  swatch.start();
  S.qprop(psi,
	  connect_state,
	  chi,
	  SUMR_INVERTER,
	  input.param.invParam.RsdCG,
	  input.param.invParam.MaxCG,
	  n_count);
  swatch.stop();
  t = swatch.getTimeInSeconds();

  QDPIO::cout << "SUMR on gaussian source took: " << n_count << " iterations" << endl;
  QDPIO::cout << "Wall clock time : " << t << " seconds" << endl;
  push(xml_out, "SUMRGaussianSource");
  write(xml_out, "n_count", n_count);
  write(xml_out, "t",       t);
  pop(xml_out);
#endif

  pop(xml_out);
  QDP_finalize();
    
  exit(0);
}
