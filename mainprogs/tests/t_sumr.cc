// $Id: t_sumr.cc,v 1.1 2004-05-12 15:46:24 bjoo Exp $

#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>

#include <cstdio>

#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

#include "chroma.h"
#include "actions/ferm/invert/invsumr.h"

using namespace QDP;
using namespace std;

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

  // Make me a linop (this callls the initialise function)
  Handle<const LinearOperator<LatticeFermion> > D_op(S.linOp(connect_state));

  // Make me an epsilon
  Handle<const LinearOperator<LatticeFermion> > g5eps(S.lgamma5epsH(connect_state));

  int G5 = Ns*Ns - 1;

  LatticeFermion psi;
  gaussian(psi);
  psi /= sqrt(norm2(psi));

  LatticeFermion chi=zero;
 
  Real mu = input.param.FermActHandle->getMass();
  Real fact = (1 + mu) / (1 - mu);

  Complex zeta = fact;

  int n_count;

  InvSUMR<LatticeFermion>(*g5eps, psi, chi, zeta, Real(1), Real(1.0e-6), 10000, n_count);

  // Recorrect normalisation
  Real ftmp = Real(2)/(1 - mu);
  chi *= ftmp;


  LatticeFermion D_chi;

  (*D_op)(D_chi, chi, PLUS);
  LatticeFermion r;
  
  r = psi - D_chi;

  Double norm_r;

  norm_r = sqrt(norm2(r));

  QDPIO::cout << " || psi - D_chi || = " << norm_r << endl;

  QDPIO::cout << " || psi - D_chi || / || psi || " << norm_r / sqrt(norm2(psi)) << endl;
  pop(xml_out);
  QDP_finalize();
    
  exit(0);
}
