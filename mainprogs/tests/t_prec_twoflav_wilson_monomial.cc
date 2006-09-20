#include "chroma.h"

using namespace Chroma;

//! To insure linking of code, place the registered code flags here
/*! This is the bit of code that dictates what fermacts are in use */
bool linkage_hack()
{
  bool foo = true;

  // 4D actions
  foo &= EvenOddPrecWilsonFermActEnv::registerAll();
 
  // 4D Monomials
  foo &= EvenOddPrecTwoFlavorWilsonFermMonomialEnv::registerAll();
  return foo;
}

//! Old dsdu routine
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

  // Dump output
  XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();
  push(xml_out, "t_monomial");

  // Read Parameters
  std::string monomial_name;           // String for Factory
  XMLReader param_in(Chroma::getXMLInputFileName());
  Handle<ExactMonomial<multi1d<LatticeColorMatrix>,
                       multi1d<LatticeColorMatrix> > >
    S_w;

  try { 
    
    // Snarf it all
    XMLReader paramtop(param_in, "/MonomialTest");
    
    // Get the string for the factory
    read(paramtop, "Monomial", S_w);
  }
  catch(const string& e) { 
    QDPIO::cerr << "Error reading XML" << endl;
    QDP_abort(1);
  }


  // Fictitious momenta for now
  multi1d<LatticeColorMatrix> p(Nd);
  for(int mu=0; mu<Nd; mu++) { 
    gaussian(p[mu]);
  }

  // Create a field state
  GaugeFieldState gauge_state(p,u);

  // Refresh Pseudofermions
  S_w->refreshInternalFields(gauge_state);

  // Compute Force from Monomial
  multi1d<LatticeColorMatrix> dsdq(Nd);
  S_w->dsdq(dsdq, gauge_state);

  // Compute action from monomial
  Double monomial_S = S_w->S(gauge_state);
  QDPIO::cout << "monomial_S = " << monomial_S << endl;

  // Check against normal version
  // Downcast for direct handling
  EvenOddPrecTwoFlavorWilsonFermMonomial& S_down = dynamic_cast<EvenOddPrecTwoFlavorWilsonFermMonomial&>(*S_w);

  // Use the debug accessors
  const LatticeFermion& phi = S_down.debugGetPhi();

  LatticeFermion X=zero;
  // This will repeat the solve... need to put caching in at some stage
  S_down.debugGetX(X, gauge_state);


  // Compute force the old fashioned way
  // Get at the FermAct
  const EvenOddPrecWilsonFermAct& S_f = S_down.debugGetFermAct();
 
  // Need this connec state to make a linOp and for the old
  // Force routine
  Handle<const ConnectState> c_state(S_f.createState(gauge_state.getQ()));

  // Get a linOp for the subset
  Handle<const EvenOddPrecLinearOperator<LatticeFermion, multi1d<LatticeColorMatrix> > >  M(S_f.linOp(c_state));

  // My energy measurement
  Double my_S = innerProductReal(phi,X, M->subset());
  QDPIO::cout << "My S = " << my_S << endl;
  
  // Call the old force routine
  multi1d<LatticeColorMatrix> dsdq2(Nd);

  // Call old force term
  prec_wilson_dsdu(S_f, S_f.getQuarkMass(), dsdq2, c_state, X);

  // Compare to new force 
  for(int mu=0; mu < Nd; mu++) { 
    taproj(dsdq[mu]);
    taproj(dsdq2[mu]);
    
    push(xml_out, "dsdu");
    write(xml_out, "dsdu_1", dsdq[mu]);
    write(xml_out, "dsdu_2", dsdq2[mu]);
    pop(xml_out);
    
    LatticeColorMatrix dsdu_diff=dsdq[mu] - dsdq2[mu];
    
    Double sum_diff=norm2(dsdu_diff);
    QDPIO::cout << "Mu = " << mu << " Sum Diff=" << sum_diff << endl;
    
    push(xml_out, "ForceDiff");
    write(xml_out, "mu", mu);
    write(xml_out, "dsdu_diff", dsdu_diff);
    pop(xml_out);
  }
  
  // End 
  pop(xml_out);
  xml_out.close();

    // Finish
  Chroma::finalize();

  exit(0);
}
