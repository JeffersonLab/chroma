#include "chroma.h"

#include <iostream>

using namespace Chroma;

void prec_wilson_dsdu(const EvenOddPrecWilsonFermAct& S,
		      const Real& Mass,
		      multi1d<LatticeColorMatrix> & ds_u,
		      Handle<const ConnectState> state,
		      const LatticeFermion& psi) 
{
  START_CODE();
  
  ds_u.resize(Nd);

  Real prefactor = Real(1)/(4*(Real(Nd) + Mass));
  
  LatticeColorMatrix utmp_1=zero;
  LatticeFermion phi=zero;
  LatticeFermion rho=zero;
  LatticeFermion sigma=zero;;
  
  LatticeFermion ftmp_2;
  
  // Do the usual Wilson fermion dS_f/dU
  // const LinearOperatorProxy<LatticeFermion> A(linOp(u));
  const Handle< const LinearOperator<LatticeFermion> >&  M(S.linOp(state));
  
  // Need the wilson dslash
  // Use u from state with BC's on 
  const multi1d<LatticeColorMatrix>& u = state->getLinks();
  WilsonDslash  D(u);
  
  //  phi = M(u)*psi
  
  (*M)(phi, psi, PLUS);
  
  /* rho = Dslash(0<-1) * psi */
  D.apply(rho, psi, PLUS, 0);
  
  /* sigma = Dslash_dag(0 <- 1) * phi */
  D.apply(sigma, phi, MINUS, 0);
  
  for(int mu = 0; mu < Nd; ++mu) {
      
    // ftmp_2(x) = -(psi(x) - ftmp_2(x)) = -(1 - gamma(mu))*psi( x )
    ftmp_2[rb[1]] = Gamma(1<<mu) * psi;
    ftmp_2[rb[1]]  -= psi;
    
    
    // utmp_1 = - Trace_spin [ ( 1 - gamma(mu) )*psi_{x+mu)*sigma^{dagger} ]
    //        = - Trace_spin [ sigma^{dagger} ( 1 - gamma_mu ) psi_{x+mu} ]
    utmp_1[rb[0]] = -traceSpin( outerProduct( shift(ftmp_2, FORWARD, mu), sigma) );
    
    
    // ftmp_2 = phi + ftmp_2 = (1 + gamma(mu))*phi( x) 
    ftmp_2[rb[1]] = Gamma(1<<mu) * phi;
    ftmp_2[rb[1]] += phi;
    
    // utmp_1 += ( 1 + gamma(mu) )*phi_{x+mu)*rho^{dagger}_x 
    utmp_1[rb[0]] += traceSpin( outerProduct( shift(ftmp_2, FORWARD, mu), rho) );
    
    // ds_u[mu][0] += u[mu][0] * utmp_1 
    //              = u[mu][0] [   ( 1 - gamma(mu) )*psi_{x+mu)*sigma^{dagger}_x
    //                           + ( 1 + gamma(mu) )*phi_{x+mu)*rho^{dagger}_x   ]
    ds_u[mu][rb[0]] = prefactor * u[mu] * utmp_1;
    
    // Checkerboard 1
    
    // ftmp_2 = -(rho - ftmp_2) = -(1 - gamma(mu))*rho( x ) 
    ftmp_2[rb[0]] = Gamma(1<<mu)*rho;
    ftmp_2[rb[0]] -= rho;
    
    // utmp_1 = ( 1 - gamma(mu) )*rho_{x+mu)*phi^{dagger}_x
    utmp_1[rb[1]] = -traceSpin( outerProduct( shift(ftmp_2, FORWARD, mu), phi) );
    
    // ftmp_2 = (gamma(mu))*sigma 
    ftmp_2[rb[0]] = Gamma(1<<mu)*sigma;
    ftmp_2[rb[0]] += sigma;
    
    
    utmp_1[rb[1]] += traceSpin( outerProduct( shift(ftmp_2, FORWARD, mu), psi) );
    ds_u[mu][rb[1]] = prefactor * u[mu] * utmp_1;
    
  }
  
  END_CODE();
}

void funky_new_dsdu(const EvenOddPrecLinearOperator<LatticeFermion, multi1d<LatticeColorMatrix> >& M,
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
  ds_u = zero;
  M.deriv(ds_u, X, Y, MINUS);

  // Do the Ydag M^dot X term 
  multi1d<LatticeColorMatrix> F_tmp(Nd);
  F_tmp = zero;
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
    Cfg_t foo; foo.cfg_type=CFG_TYPE_UNIT;
    // Cfg_t foo; foo.cfg_type=CFG_TYPE_SZIN; foo.cfg_file="./CFGIN";
    gaugeStartup(file_xml, config_xml, u, foo);
  }

  XMLFileWriter xml_out("./XMLDAT");
  push(xml_out, "t_prec_wilson_force");

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
  boundary[3] =-1;

  // Now set up a ferm act
  Handle<FermBC<LatticeFermion> > fbc(new SimpleFermBC<LatticeFermion>(boundary));
  
  // Now make up an array of handles
  Real m = Real(0.1);
  EvenOddPrecWilsonFermAct S_f_MC(fbc, m);

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

  prec_wilson_dsdu(S_f_MC, m, dsdu_1, state, psi);

  Handle<const EvenOddPrecLinearOperator<LatticeFermion, multi1d<LatticeColorMatrix> > > M_w(S_f_MC.linOp(state));


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

    // Finish
  Chroma::finalize();

  exit(0);
}

