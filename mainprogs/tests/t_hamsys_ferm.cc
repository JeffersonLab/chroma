#include "chroma.h"

#include <iostream>

using namespace Chroma;
using namespace std;

void wilson_dsdu(const UnprecWilsonFermAct& S,
		 multi1d<LatticeColorMatrix> & ds_u,
		 Handle<const ConnectState> state,
		 const LatticeFermion& psi) 
{
  START_CODE();
  
   // Get at the U matrices
  const multi1d<LatticeColorMatrix>& u = state->getLinks();

  // Get a linear operator
  Handle<const LinearOperator<LatticeFermion> > M(S.linOp(state));

  // Compute MY
  LatticeFermion Y;
  (*M)(Y, psi, PLUS);

    // Usually this is Kappa. In our normalisation it is 0.5 
    // I am adding in a factor of -1 to be consistent with the sign
    // convention for the preconditioned one. (We can always take this out
    // later
    Real prefactor=Real(0.5);

    // Two temporaries
    LatticeFermion f_tmp;
    LatticeColorMatrix u_tmp;
    for(int mu = 0; mu < Nd; mu++)
    { 
      // f_tmp = (1 + gamma_mu) Y 
      f_tmp = Gamma(1<<mu)*Y;
      f_tmp += Y;

      //   trace_spin ( ( 1 + gamma_mu ) Y_x+mu X^{dag}_x )
//      u_tmp = traceSpin(outerProduct(shift(f_tmp, FORWARD, mu),psi));
      LatticeFermion foo = shift(f_tmp, FORWARD, mu);
      u_tmp = traceSpin(outerProduct(foo,psi));

      // f_tmp = -(1 -gamma_mu) X
      f_tmp = Gamma(1<<mu)*psi;
      f_tmp -= psi;

      //  +trace_spin( ( 1 - gamma_mu) X_x+mu Y^{dag}_x)
//      u_tmp -= traceSpin(outerProduct(shift(f_tmp, FORWARD, mu),Y));
      foo = shift(f_tmp, FORWARD, mu);
      u_tmp -= traceSpin(outerProduct(foo,Y));
    
      // accumulate with prefactor
      ds_u[mu] = prefactor*( u[mu]*u_tmp );
    }
    
    END_CODE();
}

void funky_new_dsdu(const UnprecLinearOperator<LatticeFermion, multi1d<LatticeColorMatrix> >& M,
		   multi1d<LatticeColorMatrix> & ds_u,
		   const LatticeFermion& X) {
  LatticeFermion Y;
  M(Y,X,PLUS);

  // The 2 flavour derivative is:
  // -Xdag ( Mdag^{dot} M + M^dag M^{dot} ) X
  // = -Xdag Mdag^dot M X  - Xdag Mdat M^{dot}X
  // = -Xdag Mdag^dot Y  - Ydag M^dot X
  // = -( Xdag Mdag^dot Y  + Ydat M^dot X )

  // Do the X^dag Mdag^dot Y term
  M.deriv(ds_u, X, Y, MINUS);

  // Do the Ydag M^dot X term 
  multi1d<LatticeColorMatrix> F_tmp(Nd);
  M.deriv(F_tmp, Y, X, PLUS);

  // Add them and put in the minus sign
  for(int mu=0; mu < Nd; mu++) { 
    
    ds_u[mu] += F_tmp[mu];
    ds_u[mu] *= Real(-1);

  }


}


int main(int argc, char *argv[])
{
  // Initialise QDP
  Chroma::initialize(&argc, &argv);

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
    Cfg_t foo; foo.cfg_type=CFG_TYPE_DISORDERED;
    // Cfg_t foo; foo.cfg_type=CFG_TYPE_SZIN; foo.cfg_file="./CFGIN";
    gaugeStartup(file_xml, config_xml, u, foo);
  }

  XMLFileWriter xml_out("./XMLDAT");
  push(xml_out, "t_hamsys_ferm");

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

  Handle<const ConnectState> state(S_f_MC.createState(u));

  multi1d<LatticeColorMatrix> dsdu_1(Nd);
  multi1d<LatticeColorMatrix> dsdu_2(Nd);

  LatticeFermion psi;
  gaussian(psi);
  //  for(int mu=0; mu < Nd; mu++) { 
  //  dsdu_1[mu] = zero;
  //  dsdu_2[mu] = zero;
  // }

  //  S_f_MC.dsdu(dsdu_1, state, psi);
  //S_f_MC.dsdu2(dsdu_2, state, psi);

  wilson_dsdu(S_f_MC, dsdu_1, state, psi);

  Handle<const UnprecLinearOperator<LatticeFermion, multi1d<LatticeColorMatrix> > > M_w(S_f_MC.linOp(state));


  funky_new_dsdu(*M_w, dsdu_2, psi);

  for(int mu=0; mu < Nd; mu++) { 
    taproj(dsdu_1[mu]);
    taproj(dsdu_2[mu]);

    push(xml_out, "dsdu");
    write(xml_out, "dsdu_1", dsdu_1[mu]);
    write(xml_out, "dsdu_2", dsdu_2[mu]);
    pop(xml_out);
  
    LatticeColorMatrix dsdu_diff=dsdu_1[mu] - dsdu_2[mu];

    Double sum_diff=norm2(dsdu_diff);
    QDPIO::cout << "Mu = " << mu << " Sum Diff=" << sum_diff << endl;

    push(xml_out, "ForceDiff");
    write(xml_out, "mu", mu);
    write(xml_out, "dsdu_diff", dsdu_diff);
    pop(xml_out);
  }

  pop(xml_out);
  xml_out.close();


  /*
  Real Nd_pm = Real(Nd) + m;

  // Now the inverter params -- this will need cleaning up
  InvertParam_t  inv_params_MC;
  inv_params_MC.invType = CG_INVERTER;
  inv_params_MC.RsdCG = Real(1.0e-7);
  inv_params_MC.MaxCG = 500;

  // Now the Hamiltonian
  TwoFlavorDegenFermHamiltonian<WilsonGaugeAct,UnprecWilsonFermAct> H_MC( S_pg_MC, S_f_MC, inv_params_MC);


  GaugeFermFieldState mc_state(p,u,phi);
  MomRefreshGaussian(mc_state, H_MC);
  PseudoFermionHeatbath(mc_state, H_MC);

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

    // -1 because it belongs to the previous trajectory
    write(monitorHMC, "HMCtraj", HMC.getTrajNum()-1);

    Double w_plaq, s_plaq,t_plaq, link;
    MesPlq(mc_state.getQ(), w_plaq, s_plaq, t_plaq, link);
    write(monitorHMC, "w_plaq", w_plaq);

    pop(monitorHMC);

    QDPIO::cout << "Traj: " << HMC.getTrajNum()-1 << " w_plaq = " << w_plaq << endl;
    
  }

  pop(monitorHMC);
  monitorHMC.close();
  // Finish.
  */

    // Finish
  Chroma::finalize();

  exit(0);
}

