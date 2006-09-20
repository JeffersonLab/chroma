#include "chroma.h"

using namespace Chroma;

//! To insure linking of code, place the registered code flags here
/*! This is the bit of code that dictates what fermacts are in use */
bool linkage_hack()
{
  bool foo = true;


  // GaugeBC's
  foo &= SimpleGaugeBCEnv::registerAll();
  foo &= PeriodicGaugeBCEnv::registerAll();

  // GaugeActs
  foo &= WilsonGaugeActEnv::registerAll();

  // Gauge Monomials
  foo &= WilsonGaugeMonomialEnv::registerAll();

  // 4D Ferm actions
  foo &= EvenOddPrecWilsonFermActEnv::registerAll();
 
  // 4D Ferm Monomials
  foo &= EvenOddPrecTwoFlavorWilsonFermMonomialEnv::registerAll();
  return foo;
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
    Cfg_t foo; foo.cfg_type=CFG_TYPE_UNIT;;
    // Cfg_t foo; foo.cfg_type=CFG_TYPE_SZIN; foo.cfg_file="./CFGIN";
    gaugeStartup(file_xml, config_xml, u, foo);
  }

  // Dump output
  XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();
  XMLFileWriter& xml_log = Chroma::getXMLLogInstance();
  push(xml_out, "t_gauge_ferm_monomials");
  push(xml_log, "t_gauge_ferm_monomials");

  // Read Parameters
  multi1d<int> boundary(Nd);           // Ferm BC's
  std::string monomial_name;           // String for Factory
  XMLReader param_in(Chroma::getXMLInputFileName());
  // Snarf it all
  XMLReader paramtop(param_in, "/HamiltonianTest");
    

  Handle< ExactLatColMatHamiltonian > H_handle;

  try { 

    read(paramtop, "./Hamiltonian", H_handle);

  }
  catch( const std::string& e) { 
    QDPIO::cerr << "Error Reading Hamiltonian " << e <<  endl;
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
  H_handle->refreshInternalFields(gauge_state);

  // Compute Force from Monomial
  multi1d<LatticeColorMatrix> dsdq(Nd);
  H_handle->dsdq(dsdq, gauge_state);

  // Compute action from monomial
  Double KE, PE;
  H_handle->mesE(gauge_state, KE, PE);
  QDPIO::cout << "KE = " << KE << "  PE = " << PE << endl;

  pop(xml_log);
  pop(xml_out);
  xml_out.close();

    // Finish
  Chroma::finalize();

  exit(0);
}
