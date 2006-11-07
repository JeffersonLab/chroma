#include "chroma.h"
#include <string>
#include "actions/gauge/gaugeacts/rect_gaugeact.h"
#include "actions/gauge/gaugeacts/plaq_plus_spatial_two_plaq_gaugeact.h"
using namespace Chroma;

//! To insure linking of code, place the registered code flags here
/*! This is the bit of code that dictates what fermacts are in use */
bool linkageHack(void)
{
  bool foo = true;
    
  // Gauge Monomials
  foo &= GaugeMonomialEnv::registerAll();
    
  // Ferm Monomials
  foo &= WilsonTypeFermMonomialAggregrateEnv::registerAll();
    
  // MD Integrators
  foo &= LCMMDIntegratorAggregateEnv::registerAll();
    
  // Chrono predictor
  foo &= ChronoPredictorAggregrateEnv::registerAll();

  // Inline Measurements
  foo &= InlineAggregateEnv::registerAll();

  return foo;
}

int main(int argc, char *argv[]) 
{
  Chroma::initialize(&argc, &argv);
  
  START_CODE();

  // Chroma Init stuff -- Open DATA and XMLDAT
  QDPIO::cout << "Linkage = " << linkageHack() << endl;

  // Snarf it all
  XMLReader param_in(Chroma::getXMLInputFileName());
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
  XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();
  XMLFileWriter& xml_log = Chroma::getXMLLogInstance();
  push(xml_out, "t_leapfrog");
  push(xml_log, "t_leapfrog");

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

  QDPIO::cout << "create state" << endl;

  // Create a field state
  GaugeFieldState gauge_state(p,u);
 
  
  QDPIO::cout << "exact ham" << endl;

  ExactLatColMatHamiltonian& H_exact = dynamic_cast<ExactLatColMatHamiltonian&    >(*H);
  Double KE_old, PE_old;

  QDPIO::cout << "fields" << endl;

  // Put some noise into the pseudofermions...
  H_exact.refreshInternalFields(gauge_state);

  QDPIO::cout << "mesE" << endl;

  H_exact.mesE(gauge_state, KE_old, PE_old);
  QDPIO::cout << "Initial energies: KE =" << KE_old << " PE = " << PE_old <<endl;

  QDP::StopWatch swatch;
  swatch.reset();
  swatch.start();
  // Do a trajectory
  (*MD)(gauge_state, (*MD).getTrajLength());
  swatch.stop();
  double total_time = swatch.getTimeInSeconds();
  QDPIO::cout << "Trajectory took: " << total_time << " sec" <<endl;

  QDPIO::cout << "Rect force took: " << RectGaugeActEnv::getTime() 
	      << " sec.  " << Real(100)*RectGaugeActEnv::getTime()/total_time 
	      << " % of total" << endl;

  QDPIO::cout << "Other force took: "
	      << PlaqPlusSpatialTwoPlaqGaugeActEnv::getTime() << " sec.   " 
	      << Real(100)*PlaqPlusSpatialTwoPlaqGaugeActEnv::getTime()/total_time 
	      << " % of total" << endl;

  QDPIO::cout << "Expmat took: "<< ExpMatEnv::getTime()  << " sec.   " 
	      << Real(100)*ExpMatEnv::getTime()/total_time
	      << " % of total" << endl;

  QDPIO::cout << "Reunit took: "<< ReunitEnv::getTime()  << " sec.   " 
	      << Real(100)*ReunitEnv::getTime()/total_time
	      << " % of total" << endl;

  QDPIO::cout << "Taproj took: "<< TaprojEnv::getTime()  << " sec.   " 
	      << Real(100)*TaprojEnv::getTime()/total_time
	      << " % of total" << endl;

  QDPIO::cout << "Stout Smearing Took:"<< StoutLinkTimings::getSmearingTime() << " secs.  "
	      << Real(100)*StoutLinkTimings::getSmearingTime()/total_time
	      <<" % of total " << endl;

  QDPIO::cout << "Stout Force Recursion Took:"<< StoutLinkTimings::getForceTime() << " secs.  "
	      << Real(100)*StoutLinkTimings::getForceTime()/total_time
	      <<" % of total " << endl;

  QDPIO::cout << "Stout Force Functions Took:" << StoutLinkTimings::getFunctionsTime() << " secs   "
	      << Real(100)*StoutLinkTimings::getFunctionsTime()/total_time
	      <<" % of total " << endl;

  double deficit = swatch.getTimeInSeconds() 
    - RectGaugeActEnv::getTime() 
    - PlaqPlusSpatialTwoPlaqGaugeActEnv::getTime() 
    - ExpMatEnv::getTime() 
    - ReunitEnv::getTime() 
    - TaprojEnv::getTime()
    - StoutLinkTimings::getSmearingTime()
    - StoutLinkTimings::getForceTime();

  QDPIO::cout << "Time deficit : " << deficit << " sec.  "<< Real(100)*deficit/ total_time << " % of total" << endl;

  
  Double KE_new, PE_new;
  H_exact.mesE(gauge_state, KE_new, PE_new);
  QDPIO::cout << "Final energies: KE =" << KE_new << " PE = " << PE_new <<endl;
  

  Double deltaPE = PE_new - PE_old;
  Double deltaKE = KE_new - KE_old;

  QDPIO::cout << "DeltaPE = " << deltaPE << endl;
  QDPIO::cout << "DeltaKE = " << deltaKE << endl;

  QDPIO::cout << "DeltaE = " << deltaKE + deltaPE <<endl;
  pop(xml_log);   // t_leapfrog
  pop(xml_out);   // t_leapfrog
// xml_out.close();

  END_CODE();

  // Finish
  Chroma::finalize();
  exit(0);
}
