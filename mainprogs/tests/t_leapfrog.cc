
#include "chroma.h"
#include <string>
#include "actions/gauge/gaugeacts/rect_gaugeact.h"
#include "actions/gauge/gaugeacts/plaq_plus_spatial_two_plaq_gaugeact.h"
#include "update/molecdyn/integrator/lcm_toplevel_integrator.h"

// Specials
#include "update/molecdyn/hamiltonian/exact_hamiltonian.h"

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
  foo &= LCMMDComponentIntegratorAggregateEnv::registerAll();
    
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

  {
    bool monitorForcesP = true;

    if( paramtop.count("./MonitorForces") == 1 ) {
      read(paramtop, "./MonitorForces", monitorForcesP );
    }

    QDPIO::cout << "MonitorForces is " << monitorForcesP << endl;
    setForceMonitoring( monitorForcesP );
  }

  // Try and create an array of monomials:
  try { 
    readNamedMonomialArray(paramtop, "./Monomials");
  }
  catch(const std::string& e) { 
    QDPIO::cout << "Failed to read monomials " << endl;
    QDP_abort(1);
  }

  ExactHamiltonianParams ham_params(paramtop, "./Hamiltonian");
  ExactHamiltonian H(ham_params);

  // create toplevel integrator
  LCMToplevelIntegratorParams int_par(paramtop, "./MDIntegrator");
  LCMToplevelIntegrator the_integrator(int_par);

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


  Double KE_old, PE_old;

  QDPIO::cout << "fields" << endl;


  H.refreshInternalFields(gauge_state);

  QDPIO::cout << "mesE" << endl;

  H.mesE(gauge_state, KE_old, PE_old);

  QDPIO::cout << "Initial energies: KE =" << KE_old << " PE = " << PE_old <<endl;

  QDPIO::cout << "Copying copy list" << endl;
  // Setup fields 
  the_integrator.copyFields();

  QDPIO::cout << "Performing tajectory" << endl;
  QDP::StopWatch swatch;
  swatch.reset();
  swatch.start();
  // Do a trajectory
  the_integrator(gauge_state, the_integrator.getTrajLength());
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

  QDPIO::cout << "Measuring final energies" << endl;
  
  Double KE_new, PE_new;
  H.mesE(gauge_state, KE_new, PE_new);
  QDPIO::cout << "Final energies: KE =" << KE_new << " PE = " << PE_new <<endl;
  

  Double deltaPE = PE_new - PE_old;
  Double deltaKE = KE_new - KE_old;

  QDPIO::cout << "DeltaPE = " << deltaPE << endl;
  QDPIO::cout << "DeltaKE = " << deltaKE << endl;

  QDPIO::cout << "DeltaE = " << deltaKE + deltaPE <<endl;
  push(xml_log, "DeltaE");
   write(xml_log, "DeltaKE", deltaKE);
   write(xml_log, "DeltaPE", deltaPE);
   write(xml_log, "DeltaH", (deltaKE+deltaPE));
  pop(xml_log);
  pop(xml_log);   // t_leapfrog
  pop(xml_out);   // t_leapfrog
// xml_out.close();

  END_CODE();

  // Finish
  Chroma::finalize();
  exit(0);
}
