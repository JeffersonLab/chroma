#include "chroma.h"

#include <iostream>

using namespace Chroma;
using namespace std;


int main(int argc, char *argv[])
{
  // Initialise QDP
  Chroma::initialize(&argc, &argv);

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
    Cfg_t foo; foo.cfg_type=CFG_TYPE_DISORDERED;
    // Cfg_t foo; foo.cfg_type=CFG_TYPE_SZIN; foo.cfg_file="./CFGOUT";
    gaugeStartup(file_xml, config_xml, u, foo);
  }

    
  multi1d<LatticeColorMatrix> p(Nd);

#if 0 
  for(int mu = 0; mu < Nd; mu++) { 
    
    gaussian(p[mu]);
    taproj(p[mu]);
  
  }
#endif 

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
  ExactPureGaugeHamiltonian H_MC(S_pg_MC);
  ExactPureGaugeHamiltonian H_MD(S_pg_MD);

  // Generate the symplectic updates with respect to H_MD
  PureGaugeSympUpdates leaps(H_MD);

  // Test the integrator -- for energy conservation, and reversibility 
#if 0
  XMLFileWriter lf_xml("./LEAPFROG_TESTS");
  push(lf_xml, "LeapFrogTests");

  // Energy conservation
  // Do trajectories of length 1, changing dtau from 0.01 to 0.2
  for(Real dtau=0.01; toBool(dtau < 0.2); dtau +=Real(0.01)) {
    Real tau = Real(1);
    PureGaugePQPLeapFrog lf(leaps, dtau, tau);
    PureGaugeFieldState old_state(p, u);
    PureGaugeFieldState working_state(p, u);

    Double KE_old, PE_old;
    H_MD.mesE(working_state, KE_old, PE_old);
 
    // DO a traj
    lf(working_state, lf_xml);

    Double KE_new, PE_new;
    H_MD.mesE(working_state, KE_new, PE_new);

    // Flip Momenta
    for(int mu = 0; mu < Nd; mu++) { 
      working_state.getP()[mu] *= Real(-1);
    }

    // Do reverse traj
    // DO a traj
    lf(working_state, lf_xml);

    // Flip Momenta
    for(int mu = 0; mu < Nd; mu++) { 
      working_state.getP()[mu] *= Real(-1);
    }
    
    Double KE_new2, PE_new2;
    H_MD.mesE(working_state, KE_new2, PE_new2);

    Double deltaKE = KE_new - KE_old;
    Double deltaPE = PE_new - PE_old;
    Double deltaH  = deltaKE + deltaPE;
    Double ddKE = KE_new2 - KE_old;
    Double ddPE = PE_new2 - PE_old;
    Double ddH  = ddKE - ddPE;
    
    push(lf_xml, "elem");
    write(lf_xml, "tau", tau);
    write(lf_xml, "dt", dtau);
    write(lf_xml, "delH", deltaH);
    write(lf_xml, "delKE", deltaKE);
    write(lf_xml, "delPE", deltaPE);
    write(lf_xml, "delDelH", ddH);
    write(lf_xml, "ddPE", ddPE);
    write(lf_xml, "ddKE", ddKE);
    pop(lf_xml);

    QDPIO::cout << " dt = " << dtau << " deltaH = " << deltaH <<  endl;
    QDPIO::cout << "       delta KE = " << deltaKE 
		<< "       delta PE = " << deltaPE 
                << "       dPE/dKE = " << deltaPE/deltaKE << endl;

    QDPIO::cout << "       delta delta H  = " << ddH 
		<< "       delta delta KE = " << ddKE
		<< "       delta delta PE = " << ddPE<< endl << endl;
  }    
  pop(lf_xml);
  lf_xml.close();
#endif 
  
  // Step Sizes
  Real tau = Real(1);
  Real dt = Real(0.1);
  PureGaugePQPLeapFrog leapfrog(leaps, dt, tau);
  

  // Create the HMC 
  PureGaugeHMCTraj HMC(H_MC, leapfrog);
  
  XMLFileWriter monitorHMC("HMC");
  push(monitorHMC, "HMCTest");
  PureGaugeFieldState mc_state(p, u);
 


  for(int i=0; i < 100000; i++) { 

    // Do trajectory
    HMC(mc_state, true, monitorHMC);

    // Do measurements
    push(monitorHMC, "Measurements");

    // -1 because it is the last trajectory
    write(monitorHMC, "HMCtraj", HMC.getTrajNum()-1);

    Double w_plaq, s_plaq,t_plaq, link;
    MesPlq(mc_state.getQ(), w_plaq, s_plaq, t_plaq, link);
    write(monitorHMC, "w_plaq", w_plaq);

    pop(monitorHMC);

    QDPIO::cout << "Traj: " << HMC.getTrajNum()-1 << " w_plaq = " << w_plaq << endl;
    
  }

  pop(monitorHMC);
  monitorHMC.close();
  // Finish

  Chroma::finalize();
  exit(0);
}

