#include "chroma.h"
#include "update/field_state.h"
#include "update/hamiltonian.h"
#include "update/abs_symp_updates.h"
#include "update/symp_updates.h"
#include "update/hyb_int.h"
#include "update/hmc_classes.h"
#include "update/pg_hmc.h"

#include <iostream>

using namespace QDP;
using namespace std;


int main(int argc, char *argv[])
{
  // Initialise QDP
  QDP_initialize(&argc, &argv);

  // Setup a small lattice
  const int nrow_arr[] = {4, 4, 4, 4};
  multi1d<int> nrow(Nd);
  nrow=nrow_arr;
  Layout::setLattSize(nrow);
  Layout::create();
  
  multi1d<LatticeColorMatrix> u(Nd);
  {
    XMLReader file_xml;
    XMLReader config_xml;
    Cfg_t foo; foo.cfg_type=CFG_TYPE_UNIT;
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

  Real betaMC = Real(5.7);
  Real betaMD = Real(5.7);

  // Get a Wilson Gauge Action
  // For the Monte Carlo
  WilsonGaugeAct S_pg_MC(gbc, betaMC);

  // For the MD
  WilsonGaugeAct S_pg_MD(gbc, betaMD);

  // Generate Hamiltonians 
  ExactPureGaugeHamiltonian<WilsonGaugeAct> H_MC(S_pg_MC);
  ExactPureGaugeHamiltonian<WilsonGaugeAct> H_MD(S_pg_MD);

  // Generate the symplectic updates with respect to H_MD
  PureGaugeSympUpdates leaps(H_MD);

  // Step Sizes
  Real dt = Real(0.01);
  Real tau = Real(0.1);

  PureGaugePQPLeapFrog leapfrog(leaps, dt, tau);

  // Create the HMC 
  PureGaugeHMCTraj HMC(H_MC, leapfrog);

  XMLFileWriter monitorHMC("FOO");
  push(monitorHMC, "HMCTest");

  // Thermalise the HMC always accepting
  for(int i=0; i < 100; i++) { 
    HMC(mc_state, true, monitorHMC);
    MesPlq(mc_state.getQ(), w_plaq, s_plaq, t_plaq, link);
    QDPIO::cout << "Traj: " << HMC.getTrajNum() << " w_plaq=" << w_plaq << endl;
  }



  pop(monitorHMC);
  // Finish
  QDP_finalize();
  exit(0);
}

