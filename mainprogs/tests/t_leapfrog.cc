#include "chroma.h"
#include <string>

using namespace QDP;
using namespace Chroma;
using namespace std;

//! To insure linking of code, place the registered code flags here
/*! This is the bit of code that dictates what fermacts are in use */
bool linkage_hack()
{
  bool foo = true;


  // GaugeBC's
  foo &= SimpleGaugeBCEnv::registered;
  foo &= PeriodicGaugeBCEnv::registered;

  // GaugeActs
  foo &= WilsonGaugeActEnv::registered;

  // Gauge Monomials
  foo &= WilsonGaugeMonomialEnv::registered;

  // 4D Ferm actions
  foo &= EvenOddPrecWilsonFermActEnv::registered;
 
  // 4D Ferm Monomials
  foo &= EvenOddPrecTwoFlavorWilsonFermMonomialEnv::registered;

  // MD Integrators
  foo &= LatColMatPQPLeapfrogIntegratorEnv::registered;
  return foo;
}


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
    Cfg_t foo; foo.cfg_type=CFG_TYPE_DISORDERED;
    // Cfg_t foo; foo.cfg_type=CFG_TYPE_SZIN; foo.cfg_file="./CFGIN";
    gaugeStartup(file_xml, config_xml, u, foo);
  }

  // Dump output
  XMLFileWriter xml_out("./XMLDAT");
  push(xml_out, "t_gauge_ferm_monomials");

  // Read Parameters
  multi1d<int> boundary(Nd);           // Ferm BC's
  std::string monomial_name;           // String for Factory
  XMLReader param_in("DATA");
  // Snarf it all
  XMLReader paramtop(param_in, "/HamiltonianTest");
  Handle<AbsHamiltonian<multi1d<LatticeColorMatrix>,multi1d<LatticeColorMatrix> > > H(new ExactLatColMatHamiltonian(paramtop, "./Hamiltonian"));

  std::string integrator_name;
  try { 
    read(paramtop, "./MDIntegrator/Name", integrator_name);
  }
  catch(const std::string& e) { 
    QDPIO::cout << "Error reading XML: " << integrator_name << endl;
    QDP_abort(1);
  }


  // Get the Leapfrog 
  Handle< AbsMDIntegrator<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > MD( TheMDIntegratorFactory::Instance().createObject(integrator_name, paramtop, "./MDIntegrator", H) );

  // Fictitious momenta for now
  multi1d<LatticeColorMatrix> p(Nd);

  // Get some noise into the momenta...
  for(int mu=0; mu<Nd; mu++) { 
    gaussian(p[mu]);
    p[mu] *= sqrt(0.5);
    taproj(p[mu]);
  }

  // Create a field state
  GaugeFieldState gauge_state(p,u);
 
  
  ExactLatColMatHamiltonian& H_exact = dynamic_cast<ExactLatColMatHamiltonian&    >(*H);
  Double KE_old, PE_old;

  // Put some noise into the pseudofermions...
  H_exact.refresh(gauge_state);
  H_exact.mesE(gauge_state, KE_old, PE_old);
  QDPIO::cout << "Initial energies: KE =" << KE_old << " PE = " << PE_old <<endl;

  // Do a trajectory
  (*MD)(gauge_state);

  Double KE_new, PE_new;
  H_exact.mesE(gauge_state, KE_new, PE_new);
  QDPIO::cout << "Final energies: KE =" << KE_new << " PE = " << PE_new <<endl;
  

  Double deltaPE = PE_new - PE_old;
  Double deltaKE = KE_new - KE_old;

  QDPIO::cout << "DeltaPE = " << deltaPE << endl;
  QDPIO::cout << "DeltaKE = " << deltaKE << endl;

  QDPIO::cout << "DeltaE = " << deltaKE + deltaPE <<endl;
  pop(xml_out);
  xml_out.close();

    // Finish
  QDP_finalize();

  exit(0);
}
