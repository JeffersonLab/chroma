#include "chroma.h"
#include "update/field_state.h"
#include "update/hamiltonian.h"
#include "update/abs_symp_updates.h"
#include "update/symp_updates.h"
#include "update/hyb_int.h"
#include <iostream>

using namespace QDP;
using namespace std;


int main(int argc, char *argv[])
{
  // Initialise QDP
  QDP_initialize(&argc, &argv);

  // Setup a small lattice
  const int nrow_arr[] = {4, 4, 4, 8};
  multi1d<int> nrow(Nd);
  nrow=nrow_arr;
  Layout::setLattSize(nrow);
  Layout::create();
  
  multi1d<LatticeColorMatrix> u(Nd);
  {
    XMLReader file_xml;
    XMLReader config_xml;
    Cfg_t foo; foo.cfg_type=CFG_TYPE_DISORDERED;
    gaugeStartup(file_xml, config_xml, u, foo);
  }
   
  multi1d<LatticeColorMatrix> p(Nd);
  
  Double w_plaq, s_plaq, t_plaq, link;
  MesPlq(u, w_plaq, s_plaq, t_plaq, link);
  QDPIO::cout << "w_plaq before mc_state = " << w_plaq << endl;

  PureGaugeFieldState mc_state(p, u);
  MesPlq(mc_state.getQ(), w_plaq, s_plaq, t_plaq, link);
  QDPIO::cout << "w_plaq after mc_state = " << w_plaq << endl;


  // Get Periodic Gauge Boundaries
  Handle<GaugeBC> gbc(new PeriodicGaugeBC);

  Real betaMC = 5.9;
  Real betaMD = 6.0;

  // Get a Wilson Gauge Action
  WilsonGaugeAct S_pg_MC(gbc, betaMC);
  WilsonGaugeAct S_pg_MD(gbc, betaMD);

  ExactPureGaugeHamiltonian<WilsonGaugeAct> H_MC(S_pg_MC);
  ExactPureGaugeHamiltonian<WilsonGaugeAct> H_MD(S_pg_MD);

  Double KE, PE;
  H_MC.mesE(mc_state, KE, PE);
  QDPIO::cout << "MesE_MC: KE= " << KE << " PE=" << PE << endl;

  // Instantiate a SympUpdates class
  // Same Hamiltonian as HMC for now.
  PureGaugeSympUpdates leaps(H_MC);

  // Step Sizes
  Real dt = 0.1;
  Real dtby2 = dt/Real(2);

  
  // Prototyp HMC Step

  // Refresh Momenta
  H_MC.refreshP(mc_state);

  // Save State
  PureGaugeFieldState old_state(mc_state);
  PureGaugePQPLeapFrog<PureGaugeSympUpdates> leapfrog(leaps, Real(0.1), Real(1));

  XMLBufferWriter foo;
  leapfrog(mc_state, foo);


  Double KE_old, PE_old;
  H_MC.mesE(old_state, KE_old, PE_old);
  H_MC.mesE(mc_state,  KE, PE);

  QDPIO::cout << "OLD   KE: " << KE_old << " PE: " << PE_old << endl;
  QDPIO::cout << "NEW   KE: " << KE     << " PE: " << PE << endl;
  QDPIO::cout << "Delta KE: " << KE - KE_old << " PE: " << PE - PE_old << endl;

  XMLFileWriter jim("OUT");
  jim << foo;

  // Finish
  QDP_finalize();
  exit(0);
}

