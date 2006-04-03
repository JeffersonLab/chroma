// $Id: t_g5eps_bj.cc,v 3.0 2006-04-03 04:59:14 edwards Exp $

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
  Chroma::initialize(&argc, &argv);

  App_input_t input;
  XMLReader xml_in(Chroma::getXMLInputFileName());

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

  XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();
  push(xml_out,"t_g5eps_bj");


  // Measure the plaquette on the gauge
  MesPlq(xml_out, "Observables", u);
  xml_out.flush();

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

  LatticeFermion chi1, chi2, tmp, g5_psi;
  (*D_op)(chi1, psi, PLUS);
  (*g5eps)(chi2, psi, PLUS);

  // chi2 = gamma_5 eps(H) psi
  //
  // now form (1/2)(( 1 + mu ) + (1-mu)gamma_5 eps(H)) psi
  chi2 *= (Real(1) - input.param.FermActHandle->getMass());
  chi2 += (Real(1) + input.param.FermActHandle->getMass())*psi;
  chi2 *= 0.5;

  
  
  LatticeFermion diff = chi2 - chi1;
  Double norm_diff = sqrt(norm2(diff));

  QDPIO::cout << "||PLUS: lovlapms - explicit epsilon construct ||=" << norm_diff << endl;

  (*D_op)(chi1, psi, MINUS);

  g5_psi = Gamma(G5)*psi;

  (*g5eps)(tmp, g5_psi, PLUS);

  // Do 
  // chi2 = 1/2 ( ( 1 + mu) psi + (1 - mu) epsilon gamma_5 psi
  //
  //      = 1/2 ( ( 1 + mu ) psi + (1 - mu ) gamma_5 g5_eps gamma_5 psi)
  //      = 1/2 ( ( 1 + mu ) psi + (1 - mu ) gamma_5 tmp;
  chi2 = Gamma(G5)*tmp;
  chi2 *= (Real(1) - input.param.FermActHandle->getMass());
  chi2 += (Real(1) + input.param.FermActHandle->getMass())*psi;
  chi2 *= 0.5;

  diff = chi1 - chi2;
  norm_diff = sqrt(norm2(diff));

  QDPIO::cout << "||MINUS: lovlapms - explicit epsilon construct ||=" << norm_diff << endl;


  // Now do U^dagger = epsilon^dagger gamma_5 
  (*g5eps)(chi2, psi, MINUS);
  
  // chi2 = 1/2( (1 + mu) + (1 - mu) epsilon^{dagger} gamma_5 psi
  //      = 1/2( (1 + mu) + (1 - mu) gamma_5 epsilon gamma_5 gamma_5
  chi2 *= (Real(1) - input.param.FermActHandle->getMass());
  chi2 += (Real(1) + input.param.FermActHandle->getMass())*psi;
  chi2 *= 0.5;
  
  diff = chi1 - chi2;
  norm_diff = sqrt(norm2(diff));

  QDPIO::cout << "||MINUS2: lovlapms - explicit epsilon construct ||=" << norm_diff << endl;


  // Test unitarity
  (*g5eps)(chi1, psi, PLUS);
  (*g5eps)(chi2, chi1, MINUS);

  diff = psi - chi2;
  norm_diff = sqrt(norm2(diff));

  QDPIO::cout << "||UNITARITY: 1 - U^{dag} U ||=" << norm_diff << endl;

  (*g5eps)(chi1, psi, PLUS);
  (*g5eps)(chi2, psi, MINUS);
  diff = chi1 - chi2;
  norm_diff = sqrt(norm2(diff));

  QDPIO::cout << "||HERMITICITY: U - U^{dag} ||=" << norm_diff << endl;


  chi1 = Gamma(G5)*psi;
  (*g5eps)(tmp, chi1, PLUS);
  chi1 = Gamma(G5)*tmp;

  (*g5eps)(chi2, psi, MINUS);
  diff = chi1 - chi2;
  norm_diff = sqrt(norm2(diff));

  QDPIO::cout << "||gamma_5 HERMITICITY: g5 U g5 - U^+||=" << norm_diff << endl;

  chi1 = Gamma(G5)*psi;
  (*g5eps)(tmp, chi1, MINUS);
  chi1 = Gamma(G5)*tmp;

  (*g5eps)(chi2, psi, PLUS);
  diff = chi1 - chi2;
  norm_diff = sqrt(norm2(diff));

  QDPIO::cout << "||gamma_5 HERMITICITY: g5 U^+ g5  - U ||=" << norm_diff << endl;

  // Get g5 (g5 eps) psi = eps psi
  (*g5eps)(tmp, psi, PLUS);
  chi1 = Gamma(G5)*tmp;

  // Get (g5 eps)^{dag} g5 psi = eps^{dag} psi
  tmp = Gamma(G5)*psi;
  (*g5eps)(chi2, tmp, MINUS);
  
  diff = chi1 - chi2;
  norm_diff = sqrt(norm2(diff));

  QDPIO::cout << "|| sgn HERMITICITY: eps - eps^{dag} || = " << norm_diff <<endl;


  // Chi2 is eps^{dag}
  // get tmp = gamma_5 eps eps^{dag} psi
  (*g5eps)(tmp, chi2, PLUS);

  // chi2 = gamma_5 tmp = eps eps^{dag} psi
  chi2 = Gamma(G5)*tmp;

  diff = psi - chi2;
  norm_diff = sqrt(norm2(diff));

  QDPIO::cout << "||sgn UNITARITY: 1 - eps eps^{dag} || = " << norm_diff << endl;
  pop(xml_out);
  Chroma::finalize();
    
  exit(0);
}
