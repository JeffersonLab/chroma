// $Id: t_msumr.cc,v 1.1 2004-05-13 13:34:49 bjoo Exp $

#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>

#include <cstdio>

#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

#include "chroma.h"
#include "actions/ferm/invert/minvsumr.h"

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




  Handle<const LinearOperator<LatticeFermion> > g5eps(S.lgamma5epsH(connect_state));

  int numroot=3;
  int n_count;

  multi1d<Complex> zeta(numroot);
  multi1d<Real> rho(numroot);
  multi1d<LatticeFermion> x(numroot);
  multi1d<Real> epsilon(numroot);
  LatticeFermion b;

  zeta[0] = Real(1.01);
  rho[0]  = Real(1);
  epsilon[0] = Real(1.0e-4);
  for(int i=1; i < numroot; i++) { 
    zeta[i] = zeta[0] + Complex(Real(i)*Real(0.04));
    rho[i]  = rho[0];
    epsilon[i] = epsilon[i-1]*Real(0.1);
  }
  
  gaussian(b);
  b /= sqrt(norm2(b));


  
  MInvSUMR(*g5eps,
	   b,
	   x,
	   zeta,
	   rho,
	   epsilon,
	   Real(1.0e-6),
	   100,
	   n_count);

  // Check back the solutions
  for(int i=0; i < numroot; i++) { 
    LatticeFermion t1;

    (*g5eps)(t1, x[i], PLUS);
    t1 *= rho[i];
    t1 += zeta[i]*x[i];
    
    t1 -= b;
    QDPIO::cout << "|| b - A x || = " << sqrt(norm2(t1)) << endl;
    QDPIO::cout << "|| b - A x || / || b || = " << sqrt(norm2(t1))/sqrt(norm2(b)) << endl;

  }

  QDP_finalize();
    
  exit(0);
}
