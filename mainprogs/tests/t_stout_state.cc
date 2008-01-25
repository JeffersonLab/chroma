// $Id: t_stout_state.cc,v 3.11 2008-01-25 22:23:24 edwards Exp $

#include <iostream>
#include <cstdio>

#include "chroma.h"

#include "actions/ferm/fermstates/stout_fermstate_params.h"
#include "actions/ferm/fermstates/stout_fermstate_w.h"
#include "actions/ferm/invert/invcg2.h"

using namespace Chroma;

typedef LatticeFermion T;
typedef multi1d<LatticeColorMatrix> P;
typedef multi1d<LatticeColorMatrix> Q;

int main(int argc, char *argv[])
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);

  // Setup the layout
  const int foo[] = {4,4,4,4};
  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements


  Real rho=0.22;
  int  n_smear=1;
  int orthog_dir=15;

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

  XMLReader file_xml, record_xml;
  Cfg_t cfg;
  cfg.cfg_type=CFG_TYPE_WEAK_FIELD;

  gaugeStartup(file_xml, record_xml, u, cfg);

  // Try out the plaquette routine
  MesPlq(u, w_plaq, s_plaq, t_plaq, link);
  QDPIO::cout << "w_plaq = " << w_plaq << endl;
  QDPIO::cout << "link = " << link << endl;

  // Write out the results
  push(xml,"Observables");
  write(xml,"w_plaq", w_plaq);
  write(xml,"link", link);
  pop(xml);
   

  // -----------------  CHECK SMEARING ----------------------------------------
  QDPIO::cout << endl << "Stout Smearing Checks " << endl;


  push(xml, "SmearingParams");
  write(xml, "rho", rho);
  write(xml, "n_smear", n_smear);
  write(xml, "orthog_dir", orthog_dir);
  pop(xml);

  // smeared and unsmeared gauge fields
  multi1d<LatticeColorMatrix> u_smear(Nd);
  multi1d<LatticeColorMatrix> u_tmp(Nd);

  // Setup stout state smearing params.
  StoutFermStateParams s_p;
  s_p.n_smear = n_smear;
  s_p.rho.resize(Nd, Nd);
  s_p.smear_in_this_dirP.resize(Nd);

  for(int mu=0; mu < Nd; mu++) { 
    for(int nu=0; nu < Nd; nu++) {
      if( mu != nu) { 
	s_p.rho(mu,nu) = rho;
      }
      else {
	s_p.rho(mu,nu) = 0;
      }
    }

    s_p.smear_in_this_dirP[mu] = ( mu == orthog_dir ) ? false : true;
  }

  // Get the unsmeared fields into u_tmp
  u_tmp = u;
  for(int i=0; i < n_smear; i++) {
    for(int mu=0; mu < Nd; mu++){
      if( mu != orthog_dir) { 
	Stouting::stout_smear(u_smear[mu], u_tmp, mu, s_p.smear_in_this_dirP, s_p.rho);
      }
      else {
	u_smear[mu] = u_tmp[mu];
      }
    }

    u_tmp = u_smear;
  }
  
  MesPlq(u_smear, w_plaq, s_plaq, t_plaq, link);
  QDPIO::cout << "w_plaq ("<< n_smear << " levels of old stout smearing) = " << w_plaq << endl;

  push(xml, "CheckStoutStateSmear");
  write(xml, "w_plaq_old_smear", w_plaq);
  pop(xml);

  // ------------------ REGRESSION TEST THE  STOUT SMEARING IN THE 
  // ------------------ STOUT STATE against the independently cosded routine


  // Random Gauge Transformed field
  multi1d<LatticeColorMatrix> u_rg(Nd);
  LatticeColorMatrix g; // Gauge transformation matrices

  // Do the gauge transformation
  u_rg = u;
  rgauge(u_rg,g);


  // Create the stout ferm states
  typedef LatticeFermion T;
  typedef multi1d<LatticeColorMatrix> P;
  typedef multi1d<LatticeColorMatrix> Q;

  Handle< FermBC<T,P,Q> > fbc( new PeriodicFermBC<T,P,Q>() );
  // Create  Periodic FermBC
 
  Handle< StoutFermState<T,P,Q> > s_state(  new StoutFermState<T,P,Q>(fbc, s_p, u) );
  Handle< StoutFermState<T,P,Q> > s_state2( new StoutFermState<T,P,Q>(fbc, s_p, u_rg) );

  

  // Get the plaquette
  MesPlq((*s_state).getLinks(), w_plaq, s_plaq, t_plaq, link);
  QDPIO::cout << "w_plaq ("<< n_smear << " levels of new stout smearing) = " << w_plaq << endl;
  write(xml, "new_smearing_from_state", w_plaq);
  
  // Try out the plaquette routine
  MesPlq((*s_state2).getLinks(), w_plaq, s_plaq, t_plaq, link);
  QDPIO::cout << "w_plaq (After gauge transf and " << n_smear << " levels new stout smearing) = " << w_plaq << endl << endl;

  write(xml, "new_smearing_from_state_gtrans", w_plaq);

  for(int mu=0; mu < Nd; mu++) {
    QDPIO::cout << "mu: " << mu << endl;
    LatticeColorMatrix Q1, Q2;
    LatticeColorMatrix QQ1,QQ2;

    LatticeColorMatrix C_tmp;
    LatticeDouble c0,c1,c0_rg,c1_rg;

    Stouting::getQsandCs(u_rg, Q2, QQ2, C_tmp, mu, s_p.smear_in_this_dirP, s_p.rho);
    Stouting::getQsandCs(u,Q1,QQ1, C_tmp, mu, s_p.smear_in_this_dirP, s_p.rho);

    QDPIO::cout << "Gauge invatiance check for Q: " << norm2( Q2-g*Q1*adj(g) ) 
		<< endl;
    QDPIO::cout << "Gauge invariance check for Q^2: " 
		<< norm2( QQ2 - g*QQ1*adj(g)) << endl;

   

    multi1d<LatticeComplex> f(3), f_rg(3);
    multi1d<LatticeComplex> b_1,b_2;
    Stouting::getFsAndBs(Q2,QQ2,f_rg,b_1,b_2,false);
    Stouting::getFsAndBs(Q1,QQ1,f,b_1,b_2,false);

    LatticeColorMatrix expiQ = (f[0]+f[1]*Q1+f[2]*QQ1);
    LatticeColorMatrix rg_expiQ    = (f_rg[0]+f_rg[1]*Q2+f_rg[2]*QQ2);
    QDPIO::cout << "Gauge invariance check for e(iQ): " 
		<< norm2( rg_expiQ - g*expiQ*adj(g)) << endl;


    LatticeColorMatrix stout_link = expiQ*u[mu];
    LatticeColorMatrix rg_stout_link = rg_expiQ*u_rg[mu];
    QDPIO::cout << "Gauge invariance check for stout_link: " 
		<< norm2( rg_stout_link - g*stout_link*shift(adj(g),FORWARD,mu) ) 
		<< endl;

    QDPIO::cout << "Diff getLink(): " << norm2(stout_link-s_state->getLinks()[mu])
		<<endl;

    QDPIO::cout << "Diff getLink() gtrans ["<<mu<<"] = " 
		<< norm2(rg_stout_link- s_state2->getLinks()[mu]) << endl;

    LatticeColorMatrix stout_smeared;
    Stouting::stout_smear(stout_smeared, u, mu, s_p.smear_in_this_dirP, s_p.rho);
    LatticeColorMatrix rg_stout_smeared;
    Stouting::stout_smear(rg_stout_smeared, u_rg, mu, s_p.smear_in_this_dirP, s_p.rho);

    QDPIO::cout << "NonStateStoutSmear - StateStoutSmear: " << norm2(stout_smeared - s_state->getLinks()[mu]) << endl;
    QDPIO::cout << "RG: NOnStateStoutSmeared -StateStoutSmeared: " <<norm2( rg_stout_smeared - s_state2->getLinks()[mu]) << endl;
    

  }
  
  // Test the forces. Assume an action of the form
  // 2 Re Tr U X where X is a fixed random SU3 matrix
  // Then the (fat) force is U (or in the case of no smearing
  // the thin force is U
  //  U Force has to transform as G(x) U adj(G+mu) (G+mu) X adj(G)
  multi1d<LatticeColorMatrix> X(Nd);
  for(int mu=0; mu < Nd; mu++) { 
    gaussian(X[mu]);
    reunit(X[mu]);
  }
  
  // Check that X commutes with U
  multi1d<LatticeColorMatrix> F1(Nd),F2(Nd);
  for(int mu=0; mu < Nd; mu++) { 
    F1[mu] = X[mu];
    F2[mu] = shift(g,FORWARD,mu)*X[mu]*adj(g);
  }

  s_state->deriv(F1);
  s_state2->deriv(F2);

  for(int mu=0; mu < Nd; mu++) { 
    QDPIO::cout << "tr(RG F - F) = "<< norm2(trace(F2[mu]-F1[mu])) << endl;
  }

  for(int mu=0; mu < Nd; mu++) { 
    QDPIO::cout << "RG F - gtrans(F) = "<< norm2(F2[mu]-g*F1[mu]*adj(g)) << endl;
  }

  
#if 0
  // Now get the forces
  multi1d<LatticeColorMatrix> fat_force1(Nd);  //  original 
  multi1d<LatticeColorMatrix> fat_force2(Nd); //  for the RG transform

  fat_force1 = zero;
  fat_force2 = zero;

  LatticeFermion phi;
  LatticeFermion X;
  LatticeFermion Y;

  // Untransformed
  gaussian(phi);
  
  Real Mass = 0.2;
  Real RsdCG=Real(1.0e-7);
  int MaxCG=200;

  // Get Force for untransformed field

  X=zero;  
  UnprecWilsonLinOp M1(s_state, Mass);
  InvCG2<T>(M1, phi, X, RsdCG, MaxCG);
  M1(Y, X, PLUS);
  M1.deriv(fat_force, X, Y, MINUS);
  
  

  // Get Force for transformed field 
  LatticeFermion phi2 = g*phi; // Do not transfrom
  X=zero;
  UnprecWilsonLinOp M2(s_state2, Mass);
  InvCG2<T>(M2, phi2, X, RsdCG, MaxCG);
  M2(Y, X, PLUS);
  M2.deriv(fat_force2, X,Y, MINUS);

  
  for(int mu=0; mu < Nd; mu++) {
    multi1d<LatticeColorMatrix>& l = fat_force;
    multi1d<LatticeColorMatrix>& l2 = fat_force2;
    
    LatticeColorMatrix tmp_m = s_state2->getLinks()[mu]*l2[mu];
    LatticeColorMatrix tmp_m2 = s_state->getLinks()[mu]*l[mu];

    // taproj(tmp_m);
    //taproj(tmp_m2);

    LatticeColorMatrix diff_mat = adj(g)*tmp_m2*shift(g,FORWARD, mu) - tmp_m;

    QDPIO::cout << "Diff ["<<mu<<"] = " << norm2(diff_mat) << endl;
  }

#endif

#if 0
  QDPIO::cout << "Force norms before derivative with respect to thin links" << endl;
  QDPIO::cout << "========================================================" << endl << endl;

  // Get Force norms
  Double F_norm;

  push(xml, "ForcesCheck");

  F_norm = norm2(fat_force1);
  QDPIO::cout << "F_norm for fat force is " << F_norm << endl;
  write(xml, "forceNormPreGaugeDeriv", F_norm);

  F_norm = norm2(fat_force2);
  QDPIO::cout << "F_norm for RG transformed fat_force is " << F_norm << endl;
  write(xml, "forceNormPreGaugeDerivGt", F_norm);


  // Now do the recursive derivative wrt thin links.
  (*s_state).deriv(fat_force);
  (*s_state2).deriv(fat_force2);

  QDPIO::cout << endl << endl;
  QDPIO::cout << "Force norms after derivative with respect to thin links" << endl;
  QDPIO::cout << "========================================================" << endl << endl;

  // Get Force norms
  F_norm = norm2(fat_force1);
  QDPIO::cout << "F_norm for fat force is " << F_norm << endl;
  write(xml, "ForceNormPostGaugeDeriv", F_norm);

  F_norm = norm2(fat_force2);
  QDPIO::cout << "F_norm for RG transformed fat_force is " << F_norm << endl;
  write(xml, "ForceNormPostGaugeDerivGt", F_norm);

  multi1d<LatticeColorMatrix> force_diff(Nd);

  for(int mu=0; mu < Nd; mu++) 
  { 
    force_diff[mu]  = fat_force1[mu] - adj(g)*fat_force2[mu]*g;
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

#endif
  pop(xml);

  // Time to bolt
  Chroma::finalize();

  exit(0);
}
