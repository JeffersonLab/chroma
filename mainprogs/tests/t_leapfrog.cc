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
  foo &= UnprecWilsonFermActEnv::registered;

  // 4D Ferm Monomials
  foo &= UnprecTwoFlavorWilsonTypeFermMonomialEnv::registered;
  foo &= EvenOddPrecTwoFlavorWilsonTypeFermMonomialEnv::registered;

  // 5D Ferm Monomials
  foo &= UnprecTwoFlavorWilsonTypeFermMonomial5DEnv::registered;
  foo &= EvenOddPrecTwoFlavorWilsonTypeFermMonomial5DEnv::registered;

  // MD Integrators
  foo &= LatColMatPQPLeapfrogIntegratorEnv::registered;
  return foo;
}


int main(int argc, char *argv[]) 
{
  // Initialise QDP
  QDP_initialize(&argc, &argv);

  // Snarf it all
  XMLReader param_in("DATA");
  XMLReader paramtop(param_in, "/LeapfrogTest");

  multi1d<int> nrow(Nd);

  try { 
    read(paramtop, "nrow", nrow);
  }
  catch(const std::string& e) { 
    QDPIO::cerr << "Unable to read nrow from XML: " << e << endl;
    QDP_abort(1);
  }
    
  Layout::setLattSize(nrow);
  Layout::create();


  // Dump output
  XMLFileWriter xml_out("./XMLDAT");
  push(xml_out, "t_leapfrog");

  // Read Parameters
  multi1d<int> boundary(Nd);           // Ferm BC's
  std::string monomial_name;           // String for Factory


  Cfg_t cfg;
  try { 
    read(paramtop, "./GaugeStartup", cfg);
  }
  catch( const std::string& e ) { 
    QDPIO::cerr << " Error reading XML " << e << endl;
    QDP_abort(1);
  }

  multi1d<LatticeColorMatrix> u(Nd);
  {
    XMLReader file_xml;
    XMLReader config_xml;

    gaugeStartup(file_xml, config_xml, u, cfg);
  }


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
  H_exact.refreshInternalFields(gauge_state);
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
