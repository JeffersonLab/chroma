#include "chroma.h"


#include <iostream>

using namespace QDP;
using namespace std;


int main(int argc, char *argv[])
{
  // Initialise QDP
  QDP_initialize(&argc, &argv);

  // Setup a small lattice
  const int nrow_arr[] = {2, 2, 2, 2};
  multi1d<int> nrow(Nd);
  nrow=nrow_arr;
  Layout::setLattSize(nrow);
  Layout::create();
  
  multi1d<LatticeColorMatrix> u(Nd);
  {
    XMLReader file_xml;
    XMLReader config_xml;
    // Cfg_t foo; foo.cfg_type=CFG_TYPE_UNIT;
    Cfg_t foo; foo.cfg_type=CFG_TYPE_SZIN; foo.cfg_file="./CFGIN";
    gaugeStartup(file_xml, config_xml, u, foo);
  }
    
  multi1d<LatticeColorMatrix> p(Nd);
  multi1d<LatticeFermion> phi(1);


  // Get Periodic Gauge Boundaries
  Handle<GaugeBC> gbc(new PeriodicGaugeBC);
  Real betaMC = Real(5.7);
  WilsonGaugeAct S_pg_MC(gbc, betaMC);

  multi1d<int> boundary(Nd);
  boundary[0] = 1;
  boundary[1] = 1;
  boundary[2] = 1;
  boundary[3] = -1;

  // Now set up a ferm act
  Handle<FermBC<LatticeFermion> > fbc(new SimpleFermBC<LatticeFermion>(boundary));
  
  // Now make up an array of handles
  Real m = Real(0.1);
  UnprecWilsonFermAct S_f_MC(fbc, m);
  
  Real Nd_pm = Real(Nd) + m;

  // Now the inverter params -- this will need cleaning up
  InvertParam_t  inv_params_MC;
  inv_params_MC.invType = CG_INVERTER;
  inv_params_MC.RsdCG = Real(1.0e-7);
  inv_params_MC.MaxCG = 500;

  // Now the Hamiltonian
  TwoFlavorDegenFermHamiltonian<WilsonGaugeAct,UnprecWilsonFermAct> H_MC( S_pg_MC, S_f_MC, inv_params_MC);


  LatticeFermion psi;
  readSzinFerm(psi, "./PSI");
  Double psi_norm = norm2(psi);
  QDPIO::cout << "|| psi ||^2 = " << norm2(psi) << endl;
  multi1d<LatticeColorMatrix> ds_u_f(Nd);
  for(int mu=0; mu < Nd; mu++) { 
    ds_u_f[mu] = zero;
  }

  Handle<const ConnectState> g_state(S_f_MC.createState(u));
  S_f_MC.dsdu(ds_u_f, g_state, psi);

  XMLFileWriter force_xml("./FERM_FORCE");
  push(force_xml, "fforce");
  write(force_xml, "ds_u_f", ds_u_f);

 
  
  GaugeFermFieldState mc_state(p,u,phi);
  MomRefreshGaussian(mc_state, H_MC);
  PseudoFermionHeatbath(mc_state, H_MC);

  XMLFileWriter monitor("FORCE MONITOR");
  multi1d<LatticeColorMatrix> force(Nd);

  GaugeFermSympUpdates leaps(H_MC); 
  XMLFileWriter lf_xml("./LEAPFROG_TESTS");
  push(lf_xml, "LeapFrogTests");

  // Energy conservation
  // Do trajectories of length 1, changing dtau from 0.01 to 0.2
  for(Real dtau=0.005; toBool(dtau < 0.4); dtau +=Real(0.005)) {
    Real tau = Real(1);
    GaugeFermPQPLeapFrog lf(leaps, dtau, tau);
    GaugeFermFieldState old_state(mc_state);
    GaugeFermFieldState working_state(mc_state);

    Double KE_old, GE_old, FE_old, PE_old;
    H_MC.mesE(working_state, KE_old, GE_old, FE_old);
    PE_old = GE_old + FE_old;

    // DO a traj
    lf(working_state, lf_xml);

    Double KE_new, GE_new, FE_new, PE_new;
    H_MC.mesE(working_state, KE_new, GE_new, FE_new);
    PE_new = GE_new + FE_new;

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
    
    Double KE_new2, GE_new2, FE_new2, PE_new2;
    H_MC.mesE(working_state, KE_new2, GE_new2, FE_new2);
    PE_new2 = GE_new2 + FE_new2;


    Double deltaKE = KE_new - KE_old;
    Double deltaGE = GE_new - GE_old;
    Double deltaFE = FE_new - FE_old;
    Double deltaPE = PE_new - PE_old;
    Double deltaH  = fabs(deltaKE + deltaPE);
    Double ddKE = KE_new2 - KE_old;
    Double ddPE = PE_new2 - PE_old;
    Double ddGE = GE_new2 - GE_old;
    Double ddFE = FE_new2 - FE_old;
    Double ddH  = ddKE + ddGE +ddFE;
    
    push(lf_xml, "elem");
    write(lf_xml, "tau", tau);
    write(lf_xml, "dt", dtau);
    write(lf_xml, "delH", deltaH);
    write(lf_xml, "delKE", deltaKE);
    write(lf_xml, "delPE", deltaPE);
    write(lf_xml, "delGE", deltaGE);
    write(lf_xml, "delFE", deltaFE);

    write(lf_xml, "ddH", ddH);
    write(lf_xml, "ddKE", ddKE);
    write(lf_xml, "ddPE", ddPE);
    write(lf_xml, "ddGE", ddGE);
    write(lf_xml, "ddFE", ddFE);
    pop(lf_xml);

    QDPIO::cout << " dt = " << dtau << " deltaH = " << deltaH <<  endl;
    QDPIO::cout << " delta KE = " << deltaKE 
                << " delta PE = " << deltaPE
		<< " delta GE = " << deltaGE 
                << " delta FE = " << deltaFE
                << " dPE/dKE = " << deltaPE/deltaKE << endl;

    QDPIO::cout << " delta delta H  = " << ddH 
		<< " delta delta KE = " << ddKE
		<< " delta delta PE = " << ddPE<< endl << endl;
  }    
  pop(lf_xml);
  lf_xml.close();

  Real tau = Real(1);
  Real dt = Real(0.1);
  GaugeFermPQPLeapFrog leapfrog(leaps, dt, tau);
  GaugeFermHMCTraj HMC(H_MC, leapfrog);
  XMLFileWriter monitorHMC("HMC");
  push(monitorHMC, "HMCTest");

  mc_state.getQ() = u;
  
  for(int i=0; i < 500; i++) {
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

 
  
   // Finish
  QDP_finalize();

  exit(0);
}

