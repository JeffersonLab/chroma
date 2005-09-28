// $Id: t_stout_state.cc,v 2.1 2005-09-28 03:24:19 bjoo Exp $

#include <iostream>
#include <cstdio>

#include "chroma.h"
#include "actions/ferm/fermacts/stout_state.h"

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

#if 1
  QDPIO::cout << "Start gaussian\n";
  for(int m=0; m < u.size(); ++m)
    gaussian(u[m]);

  // Reunitarize the gauge field
  for(int m=0; m < u.size(); ++m)
    reunit(u[m]);
#else
  {
    XMLReader gauge_xml;
    readSzin(gauge_xml, u, string("CFGIN"));
    xml << gauge_xml;
  }

  unitarityCheck(u);
#endif

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
  for(int mu = 0; mu < Nd; ++mu)
    polylp(u, pollp[mu], mu);

  QDPIO::cout << "w_plaq = " << w_plaq << endl;
  QDPIO::cout << "link = " << link << endl;

  push(xml,"Observables_after_gt");
  write(xml,"w_plaq", w_plaq);
  write(xml,"link", link);
  write(xml,"pollp", pollp);
  pop(xml);


  // Call the old stout smear routine 
  Real rho=0.45;
  int  n_smear=2;
  int  j_decay=Nd-1;

  multi1d<LatticeColorMatrix> u_stout_old(Nd);
  u_stout_old=u;

  for(int i=0; i < n_smear; i++) {
    multi1d<LatticeColorMatrix> u_tmp(Nd);

    for(int mu=0; mu < Nd; mu++) { 

      if( mu != j_decay ) {

	stout_smear(u_tmp[mu], u_stout_old, mu, 
		    rho, j_decay);
      }
      else { 
	u_tmp[mu] = u_stout_old[mu];
      }
    }

    u_stout_old = u_tmp;
  }
  
  // Try out the plaquette routine
  MesPlq(u_stout_old, w_plaq, s_plaq, t_plaq, link);
  QDPIO::cout << "w_plaq (old stout smearing) = " << w_plaq << endl;
  QDPIO::cout << "link (old stout smearing) = " << link << endl;
  
  // Test polyakov routine
  for(int mu = 0; mu < Nd; ++mu)
    polylp(u_stout_old, pollp[mu], mu);
  
  // Write out the results
  push(xml,"OldStoutedObservables");
  write(xml,"w_plaq", w_plaq);
  write(xml,"link", link);
  write(xml,"pollp", pollp);
  pop(xml);
  



  // New way 
  StoutConnectState s_state(u, rho, n_smear, j_decay);


  const multi1d<LatticeColorMatrix>& u_stout_new=s_state.getLinks();

  // Try out the plaquette routine
  MesPlq(u_stout_new, w_plaq, s_plaq, t_plaq, link);
  QDPIO::cout << "w_plaq (new stout smearing) = " << w_plaq << endl;
  QDPIO::cout << "link (new stout smearing) = " << link << endl;
 
  // Test polyakov routine
  for(int mu = 0; mu < Nd; ++mu)
    polylp(u_stout_new, pollp[mu], mu);
  
  // Write out the results
  push(xml,"NewStoutedObservables");
  write(xml,"w_plaq", w_plaq);
  write(xml,"link", link);
  write(xml,"pollp", pollp);
  pop(xml);



  pop(xml);

  // Time to bolt
  Chroma::finalize();

  exit(0);
}
