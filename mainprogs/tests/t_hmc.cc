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

namespace Chroma { 
  
  struct HMCParams { 

    multi1d<int> nrow;

    Cfg_t start_cfg;

    int n_warm_up;
    int n_prod;
    
    // Polymorphic
    std::string H_MC_xml;
    std::string H_MD_xml;
    std::string Integrator_xml;
  };

  void read(XMLReader& xml, const std::string& path, HMCParams& p) 
  {
    try {
      XMLReader paramtop(xml, path);
      
      read(paramtop, "./nrow", p.nrow);
      read(paramtop, "./GaugeStartup", p.start_cfg);
      read(paramtop, "./NumWarmUp", p.n_warm_up);
      read(paramtop, "./NumHMCTraj", p.n_prod);

      // Now the XML for the Hamiltonians
      XMLReader H_MC_xml(paramtop, "./MC_Hamiltonian");
      std::ostringstream os_H_MC;
      H_MC_xml.print(os_H_MC);
      p.H_MC_xml = os_H_MC.str();
      
      QDPIO::cout << "HMC_xml is: " << endl;
      QDPIO::cout << p.H_MC_xml;

      
      // Separate MD Hamiltonian specified
      XMLReader H_MD_xml(paramtop, "./MD_Hamiltonian");
      std::ostringstream os_H_MD;
      H_MD_xml.print(os_H_MD);
      p.H_MD_xml = os_H_MD.str();

      QDPIO::cout << "HMD_xml is: " << endl;
      QDPIO::cout << p.H_MD_xml;
      
      // Read the Integrator parameters
      XMLReader MD_integrator_xml(paramtop, "./MDIntegrator");
      std::ostringstream os_integrator;
      MD_integrator_xml.print(os_integrator);
      p.Integrator_xml = os_integrator.str();

      QDPIO::cout << "Integrator XML is: " << endl;
      QDPIO::cout << p.Integrator_xml << endl;
    }
    catch( const std::string& e ) { 
      QDPIO::cerr << "Error reading XML : " << e << endl;
      QDP_abort(1);
    }
  }
}; // End namespace 


using namespace Chroma;

int main(int argc, char *argv[]) 
{
  // Initialise QDP
  QDP_initialize(&argc, &argv);

  try {
  // Snarf it all
    XMLReader param_in("DATA");
    XMLReader paramtop(param_in, "/HMCTest");
    
    HMCParams params;
    read( paramtop, "./HMC", params);
      
    Layout::setLattSize(params.nrow);
    Layout::create();


    // Dump output
    XMLFileWriter xml_out("./XMLDAT");
    push(xml_out, "t_hmc");

    
    multi1d<LatticeColorMatrix> u(Nd);
    {
      XMLReader file_xml;
      XMLReader config_xml;
      
      gaugeStartup(file_xml, config_xml, u, params.start_cfg);
    }
    
    
    // Get the MC_Hamiltonian
    std::istringstream H_MC_is(params.H_MC_xml);
    XMLReader H_MC_xml(H_MC_is);
    Handle< ExactAbsHamiltonian< multi1d<LatticeColorMatrix>, 
      multi1d<LatticeColorMatrix> > > H_MC(new ExactLatColMatHamiltonian(H_MC_xml, "/MC_Hamiltonian"));
    
    // Get the MD_Hamiltonian
    std::istringstream H_MD_is(params.H_MD_xml);
    XMLReader H_MD_xml(H_MD_is);
    
    Handle< AbsHamiltonian< multi1d<LatticeColorMatrix>, 
      multi1d<LatticeColorMatrix> > > H_MD(new ExactLatColMatHamiltonian(H_MD_xml, "/MD_Hamiltonian"));
    
    
    std::istringstream Integrator_is(params.Integrator_xml);
    XMLReader MD_xml(Integrator_is);
    
    // Get the Integrator
    std::string integrator_name;
    read(MD_xml, "/MDIntegrator/Name", integrator_name);
      
    // Get the Leapfrog 
    Handle< AbsMDIntegrator<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > MD( TheMDIntegratorFactory::Instance().createObject(integrator_name, MD_xml, "/MDIntegrator", H_MD) );
    
    
    // Get the HMC
    LatColMatHMCTrj theHMCTrj( H_MC, MD );
    
    // Fictitious momenta for now
    multi1d<LatticeColorMatrix> p(Nd);
    
    // Create a field state
    GaugeFieldState gauge_state(p,u);
    
    int trj;
    for(trj=0; trj < params.n_warm_up; trj++) {
      
      // Do the trajectory without accepting 
      theHMCTrj( gauge_state, false );
      
      Double w_plaq, s_plaq, t_plaq, link;
      MesPlq(gauge_state.getQ(), w_plaq, s_plaq, t_plaq, link);
      push(xml_out, "HMCTrjObservables");
      write(xml_out, "traj_no", trj);
      write(xml_out, "w_plaq", w_plaq);
      pop(xml_out);
      QDPIO::cout << "w_plaq = " << w_plaq << endl;
      
    }
    
    for(; trj < params.n_prod + params.n_warm_up; trj++) { 
      
      theHMCTrj( gauge_state, true ) ;

      Double w_plaq, s_plaq, t_plaq, link;
      MesPlq(gauge_state.getQ(), w_plaq, s_plaq, t_plaq, link);
      push(xml_out, "HMCTrjObservables");
      write(xml_out, "traj_no", trj);
      write(xml_out, "w_plaq", w_plaq);
      pop(xml_out);
      QDPIO::cout << "w_plaq = " << w_plaq << endl;
      
    }
    
    pop(xml_out);
  }
  catch( const std::string& e) { 
    QDPIO::cerr << "Caught Exception: " << e << endl;
    QDP_abort(1);
  }

    // Finish
  QDP_finalize();

  exit(0);
}
