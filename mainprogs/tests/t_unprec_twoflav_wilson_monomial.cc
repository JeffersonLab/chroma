#include "chroma.h"

using namespace Chroma;

//! To insure linking of code, place the registered code flags here
/*! This is the bit of code that dictates what fermacts are in use */
bool linkage_hack()
{
  bool foo = true;

  // 4D actions
  foo &= UnprecWilsonFermActEnv::registerAll();
 
  // 4D Monomials
  foo &= UnprecTwoFlavorWilsonFermMonomialEnv::registerAll();
  return foo;
}

//! Old dsdu routine
void wilson_dsdu(const UnprecWilsonFermAct& S,
		 multi1d<LatticeColorMatrix>& ds_u,
		 multi1d<LatticeColorMatrix>& u_,
		 const LatticeFermion& psi) 
{
  START_CODE();
  
  // Apply BC's
  Handle<const ConnectState> state(S.createState(u_));

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

  Handle< ExactMonomial<multi1d<LatticeColorMatrix>, 
    multi1d<LatticeColorMatrix> > > S_w;
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
  UnprecTwoFlavorWilsonFermMonomial& S_down = dynamic_cast<UnprecTwoFlavorWilsonFermMonomial&>(*S_w);

  // Use the debug accessors
  const LatticeFermion& phi = S_down.debugGetPhi();

  LatticeFermion X=zero;
  // This will repeat the solve... need to put caching in at some stage
  S_down.debugGetX(X, gauge_state);

  // My energy measurement
  Double my_S = innerProductReal(phi,X);
  QDPIO::cout << "My S = " << my_S << endl;

  // Compute force the old fashioned way
  // Get at the FermAct
  const UnprecWilsonFermAct& S_f = S_down.debugGetFermAct();
  multi1d<LatticeColorMatrix> dsdq2(Nd);
  
  // Call the old force routine
  wilson_dsdu(S_f, dsdq2, gauge_state.getQ(), X);

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
