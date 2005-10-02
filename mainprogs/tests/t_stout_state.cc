// $Id: t_stout_state.cc,v 2.2 2005-10-02 03:08:50 bjoo Exp $

#include <iostream>
#include <cstdio>

#include "chroma.h"
#include "actions/ferm/fermacts/stout_state.h"
#include "actions/ferm/fermacts/stout_fermact_params.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"

using namespace Chroma;

int main(int argc, char *argv[])
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);

  // Setup the layout
  const int foo[] = {4,4,4,8};
  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements
  Layout::setLattSize(nrow);
  Layout::create();

  XMLFileWriter xml("t_mesplq.xml");
  push(xml, "t_mesplq");

  push(xml,"lattis");
  write(xml,"Nd", Nd);
  write(xml,"Nc", Nc);
  write(xml,"nrow", nrow);
  pop(xml);

  //! Example of calling a plaquette routine
  /*! NOTE: the STL is *not* used to hold gauge fields */
  multi1d<LatticeColorMatrix> u(Nd);
  Double w_plaq, s_plaq, t_plaq, link;

  QDPIO::cout << "Start gaussian\n";
  for(int m=0; m < u.size(); ++m)
    gaussian(u[m]);

  // Reunitarize the gauge field
  for(int m=0; m < u.size(); ++m)
    reunit(u[m]);

  // Try out the plaquette routine
  MesPlq(u, w_plaq, s_plaq, t_plaq, link);
  QDPIO::cout << "w_plaq = " << w_plaq << endl;
  QDPIO::cout << "link = " << link << endl;

  // Test polyakov routine
  multi1d<DComplex> pollp(Nd);
  for(int mu = 0; mu < Nd; ++mu)
    polylp(u, pollp[mu], mu);

  // Write out the results
  push(xml,"Observables");
  write(xml,"w_plaq", w_plaq);
  write(xml,"link", link);
  write(xml,"pollp", pollp);
  pop(xml);


  // Test gauge invariance
  rgauge(u);

  MesPlq(u, w_plaq, s_plaq, t_plaq, link);
  QDPIO::cout << "After GT w_plaq = " << w_plaq << endl;
  QDPIO::cout << "After GT link = " << link << endl;


  // Call the old stout smear routine 
  Real rho=0.1;
  int  n_smear=0;

  multi1d<LatticeColorMatrix> u_rg(Nd);
  u_rg = u;

  // New way 
  StoutConnectState s_state(u, rho, n_smear);

  //SimpleConnectState s_state(u);
  
  rgauge(u_rg);

  //SimpleConnectState s_state2(u_rg);
  // Make fat links from random gauge transformed u
  StoutConnectState s_state2(u_rg, rho, n_smear);


  // Try out the plaquette routine
  MesPlq(s_state.getLinks(), w_plaq, s_plaq, t_plaq, link);
  QDPIO::cout << "w_plaq ("<< n_smear << " levels of new stout smearing) = " << w_plaq << endl;
  QDPIO::cout << "link (" << n_smear << " levels of new stout smearing) = " << link << endl;


  // Try out the plaquette routine
  MesPlq(s_state2.getLinks(), w_plaq, s_plaq, t_plaq, link);
  QDPIO::cout << "w_plaq (After GT2  new stout smearing) = " << w_plaq << endl;
  QDPIO::cout << "link (After GT2 new stout smearing) = " << link << endl;


  multi1d<LatticeColorMatrix> fat_force(Nd);
  multi1d<LatticeColorMatrix> fat_force2(Nd);

  fat_force=0;
  fat_force2=0;

  LatticeFermion phi;
  LatticeFermion X;
  LatticeFermion Y;

  gaussian(phi);
  
  Real Mass = 0.2;
  int n_count;
  Real RsdCG=Real(1.0e-5);
  int MaxCG=200;
  // Now create a linop
  X=zero;  
  UnprecWilsonLinOp M1(s_state.getLinks(), Mass);
  InvCG2(M1, phi, X, RsdCG, MaxCG, n_count);
  QDPIO::cout << "n_count is " << n_count << endl;
  M1(Y, X, PLUS);
  M1.deriv(fat_force, X, Y, MINUS);
    
  X=zero;
  UnprecWilsonLinOp M2(s_state2.getLinks(), Mass);
  InvCG2(M2, phi, X, RsdCG, MaxCG, n_count);
  QDPIO::cout << "n_count is " << n_count << endl;
  M2(Y, X, PLUS);
  M2.deriv(fat_force2, X,Y, MINUS);
  
  Double F_norm;

  F_norm = norm2(fat_force);
  QDPIO::cout << "F_norm for fat_force is " << F_norm << endl;
  
  // Get force
  F_norm = norm2(fat_force2);
  QDPIO::cout << "RG Trans F_norm  for fat_force2 is " << F_norm << endl;
  

  // Now derive wrt fat links
  s_state.deriv(fat_force);
  s_state2.deriv(fat_force2);

  
  F_norm = norm2(fat_force);
  QDPIO::cout << "F_norm for fat_force is " << F_norm << endl;
  

  F_norm = norm2(fat_force2);
  QDPIO::cout << "RG Trans F_norm  for fat_force2 is " << F_norm << endl;
  


  pop(xml);


  
  XMLReader stout_xml("stout_xml");
  StoutFermActParams p(stout_xml, "/foo/FermionAction");
  XMLFileWriter fred("./foo.xml");
  push(fred, "foo");
  write(fred, "FermionAction", p);
  pop(fred);
  

  std::string fermion_action_name;
  read(stout_xml, "/foo/FermionAction/FermAct", fermion_action_name);
  QDPIO::cout <<"FermAct = " << fermion_action_name << endl;
  try { 
    volatile bool r =  UnprecStoutWilsonTypeFermActEnv::registered;
    r &= UnprecWilsonFermActEnv::registered;
    
    FermAct4D<LatticeFermion>* F=UnprecStoutWilsonTypeFermActEnv::createFermAct4D(stout_xml, "/foo/FermionAction");
    
    Handle< const ConnectState > state(F->createState(u));
    Handle< const LinearOperator<LatticeFermion> > M(F->linOp(state));
    
    (*M)(X,Y, PLUS);
  }
  catch( const std::string& e) { 
    QDPIO::cout << "caught: " << e << endl;
  }
  

  // Time to bolt
  Chroma::finalize();

  exit(0);
}
