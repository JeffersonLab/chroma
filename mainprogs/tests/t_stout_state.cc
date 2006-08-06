// $Id: t_stout_state.cc,v 3.1 2006-08-06 16:40:24 edwards Exp $

#include <iostream>
#include <cstdio>

#include "chroma.h"
#include "actions/ferm/fermacts/stout_fermstate_w.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermbcs/periodic_fermbc.h"

using namespace Chroma;

int main(int argc, char *argv[])
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);

  // Setup the layout
  multi1d<int> nrow;

  // Smearing Parameters 
  XMLReader xml_in(Chroma::getXMLInputFileName());

  StoutFermStateParams stoutParams;

//  stoutParams.rho = 0.1;
//  stoutParams.n_smear = 2;
//  stoutParams.orthog_dir = 2;

  try {
    read(xml_in, "/t_stout_state", stoutParams);
    read(xml_in, "/t_stout_state/nrow", nrow);
  }
  catch(const std::string& err) { 
    QDPIO::cerr << "Caught Error while reading XML: " << err << endl;
    QDP_abort(1);
  }

  Layout::setLattSize(nrow);
  Layout::create();

  XMLFileWriter xml(Chroma::getXMLOutputFileName());
  push(xml, "t_stout_state");

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


  // -----------------  CHECK SMEARING ----------------------------------------
  QDPIO::cout << endl << "Stout Smearing Checks " << endl;


  write(xml, "SmearingParams", stoutParams);

  // ------------------ REGRESSION TEST THE  STOUT SMEARING IN THE 
  // ------------------ STOUT STATE against the independently cosded routine


  // Random Gauge Transformed field
  multi1d<LatticeColorMatrix> u_rg(Nd);
  LatticeColorMatrix g; // Gauge transformation matrices

  // Do the gauge transformation
  u_rg = u;
  rgauge(u_rg,g);
  
  // Create state - both untransformed and transformed 
  typedef LatticeFermion               T;
  typedef multi1d<LatticeColorMatrix>  P;
  typedef multi1d<LatticeColorMatrix>  Q;

  Handle< FermBC<T,P,Q> > fbc(new PeriodicFermBC<T,P,Q>);
  Handle< CreateFermState<T,P,Q> > cfs(new CreateStoutFermState(fbc, stoutParams));

  Handle< FermState<T,P,Q> > s_state1((*cfs)(u));
  Handle< FermState<T,P,Q> > s_state2((*cfs)(u_rg));

  // Get the plaquette
  MesPlq(s_state1->getLinks(), w_plaq, s_plaq, t_plaq, link);
  QDPIO::cout << "w_plaq ("<< stoutParams.n_smear 
	      << " levels of new stout smearing) = " << w_plaq << endl;
  write(xml, "new_smearing_from_state", w_plaq);
  
  // Try out the plaquette routine
  MesPlq(s_state2->getLinks(), w_plaq, s_plaq, t_plaq, link);
  QDPIO::cout << "w_plaq (After gauge transf and " << stoutParams.n_smear 
	      << " levels new stout smearing) = " << w_plaq << endl << endl;
  write(xml, "new_smearing_from_state_gtrans", w_plaq);
  pop(xml);

  // Now get the forces
  multi1d<LatticeColorMatrix> fat_force(Nd);  //  original 
  multi1d<LatticeColorMatrix> fat_force2(Nd); //  for the RG transform

  fat_force=0;
  fat_force2=0;

  LatticeFermion phi;
  LatticeFermion X;
  LatticeFermion Y;

  // Untransformed
  gaussian(phi);
  
  Real Mass = 0.2;
  int n_count;
  Real RsdCG=Real(1.0e-7);
  int MaxCG=200;

  
  // Get Force for untransformed field
  X=zero;  
  UnprecWilsonLinOp M1(s_state1, Mass);
  InvCG2(M1, phi, X, RsdCG, MaxCG);
  M1(Y, X, PLUS);
  M1.deriv(fat_force, X, Y, MINUS);
  
  

  // Get Force for transformed field 
  LatticeFermion phi2 = g*phi; // Transform source fermion
  X=zero;
  UnprecWilsonLinOp M2(s_state2, Mass);
  InvCG2(M2, phi2, X, RsdCG, MaxCG);
  M2(Y, X, PLUS);
  M2.deriv(fat_force2, X,Y, MINUS);


  QDPIO::cout << "Force norms before derivative with respect to thin links" << endl;
  QDPIO::cout << "========================================================" << endl << endl;

  // Get Force norms
  Double F_norm;

  push(xml, "ForcesCheck");

  F_norm = norm2(fat_force);
  QDPIO::cout << "F_norm for fat force is " << F_norm << endl;
  write(xml, "forceNormPreGaugeDeriv", F_norm);

  F_norm = norm2(fat_force2);
  QDPIO::cout << "F_norm for RG transformed fat_force is " << F_norm << endl;
  write(xml, "forceNormPreGaugeDerivGt", F_norm);

  // Now do the recursive derivative wrt thin links.
  s_state1->deriv(fat_force);
  s_state2->deriv(fat_force2);

  QDPIO::cout << endl << endl;
  QDPIO::cout << "Force norms after derivative with respect to thin links" << endl;
  QDPIO::cout << "========================================================" << endl << endl;

  // Get Force norms
  F_norm = norm2(fat_force);
  QDPIO::cout << "F_norm for fat force is " << F_norm << endl;
  write(xml, "ForceNormPostGaugeDeriv", F_norm);

  F_norm = norm2(fat_force2);
  QDPIO::cout << "F_norm for RG transformed fat_force is " << F_norm << endl;
  write(xml, "ForceNormPostGaugeDerivGt", F_norm);

  multi1d<LatticeColorMatrix> force_diff(Nd);

  for(int mu=0; mu < Nd; mu++) { 
    force_diff[mu]  = fat_force[mu] - adj(g)*fat_force2[mu]*g;
    F_norm = sqrt(norm2(force_diff[mu]));
    QDPIO::cout << "|| force - RG force in dir "<< mu <<" ||=  "<< F_norm <<   endl;
    ostringstream tagname;
    tagname << "force_diff_norm_" << mu;
    write(xml, tagname.str(), F_norm);
  }
  pop(xml);

  F_norm = sqrt(norm2(force_diff));
  QDPIO::cout << "Total difference between original and gauge transformed force: " << F_norm << endl;
  push(xml, "ForceDiffNorm");
  write(xml, "totalForceDiffNorm", F_norm);
  pop(xml);


  pop(xml);

  // Time to bolt
  Chroma::finalize();

  exit(0);
}
