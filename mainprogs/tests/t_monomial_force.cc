
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
  XMLReader paramtop(param_in, "/MonomialTests");

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
  push(xml_out, "MonomialTimings");
  push(xml_log, "MonomialTimings");

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

  // Try and create an array of monomials:
  try { 
    readNamedMonomialArray(paramtop, "./Monomials");
  }
  catch(const std::string& e) { 
    QDPIO::cout << "Failed to read monomials " << endl;
    QDP_abort(1);
  }

  // Read the list of monomials to test.
  multi1d<std::string> monomial_test_ids;
  read(paramtop, "./MonomialTestIds", monomial_test_ids);

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

  QDPIO::cout << "There are " << monomial_test_ids.size() << " monomials to test" << endl;
  typedef multi1d<LatticeColorMatrix> P;
  typedef multi1d<LatticeColorMatrix> Q;

  multi1d< Handle< Monomial<P, Q> > > mon_handles(monomial_test_ids.size());
  Chroma::IntegratorShared::bindMonomials( monomial_test_ids, mon_handles );

  // Test the forces of each monomial:
  StopWatch swatch;
  double seconds;
  int hits;

  for( int mon = 0; mon < monomial_test_ids.size(); mon++) { 
    // Get a monomial
    Monomial<P,Q>& the_mon = *( mon_handles[mon] );
    QDPIO::cout << "Timing monomial force: " << monomial_test_ids[mon] << endl;
    QDPIO::cout << "Refreshing internal fields: " << endl;

    the_mon.refreshInternalFields(gauge_state);
    P F;

    seconds = 0;
    hits = 1;

    QDPIO::cout << "Calibrating";
    while( seconds < 30 ) { 
      hits *= 2;
      swatch.reset();
      swatch.start();
      for(int i=0; i < hits; i++) { 
	the_mon.dsdq(F,gauge_state);
      }
      swatch.stop();
      seconds = swatch.getTimeInSeconds();
      QDPInternal::globalSum(seconds);
      seconds /= (double)Layout::numNodes();
      QDPIO::cout << "." << flush;
    }

    QDPIO::cout << " " << seconds << "(s)"<< endl;

    QDPIO::cout << "Timing monomial force: " << monomial_test_ids[mon] << " with " << hits << " hits " << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < hits; i++) { 
      the_mon.dsdq(F,gauge_state);
    }
    swatch.stop();
    seconds = swatch.getTimeInSeconds(); 
    QDPInternal::globalSum(seconds);
    seconds /= (double)Layout::numNodes();
   
    QDPIO::cout << "monomial_id = " << monomial_test_ids[mon]
		<< ", Nhits = " << hits 
                << ", Time = " << seconds << "(s)" 
                << ", Time per hit = " << seconds/(double)hits << "(s)"<< endl;

    push(xml_out, "ForceTiming");
    write(xml_out, "monomial_id", monomial_test_ids[mon]);
    write(xml_out, "nhits", hits);
    write(xml_out, "time", seconds);
    write(xml_out, "time_per_hit", (seconds/(double)hits));
    pop(xml_out);
  }

#if 0  
  // Time taproj
  seconds = 0;
  hits = 1;
  QDPIO::cout << "Timing taproj()" << endl;
  QDPIO::cout << "Calibrating";
  while( seconds < 10 ) { 
    hits *= 2;
    swatch.reset();
    swatch.start();
    for(int i=0; i < hits; i++) { 
      for(int mu=0; mu < Nd; mu++) { 
	taproj(gauge_state.getP()[mu]);
      }
    }
    swatch.stop();
    seconds = swatch.getTimeInSeconds();
    QDPInternal::globalSum(seconds);
    seconds /= (double)Layout::numNodes();
    QDPIO::cout << "." << flush;
  }

  QDPIO::cout << " " << seconds << "(s)"<< endl;

  
  QDPIO::cout << "Timing taproj with " << hits << " hits " << endl;
  swatch.reset();
  swatch.start();
  for(int i=0; i < hits; i++) { 
   for(int mu=0; mu < Nd; mu++) { 
	taproj(gauge_state.getP()[mu]);
      }
  }
  swatch.stop();
  seconds = swatch.getTimeInSeconds();
  QDPInternal::globalSum(seconds);
  seconds /= (double)Layout::numNodes();

  QDPIO::cout << "taproj"
	      << ", Nhits = " << hits 
	      << ", Time = " << seconds << "(s)" 
	      << ", Time per hit = " << seconds/(double)hits << "(s)" << endl;

  push(xml_out, "TaprojTiming");
  write(xml_out, "nhits", hits);
  write(xml_out, "time", seconds);
  write(xml_out, "time_per_hit", (seconds/(double)hits));
  pop(xml_out);

  // Time expmat
  seconds = 0;
  hits = 1;
  QDPIO::cout << "Timing expmat()" << endl;
  QDPIO::cout << "Calibrating";
  
  while( seconds < 10 ) { 
    hits *= 2;
    swatch.reset();
    swatch.start();
    for(int i=0; i < hits; i++) { 
      for(int mu=0; mu < Nd; mu++) { 
	expmat(gauge_state.getP()[mu], EXP_EXACT);
      }
    }
    swatch.stop();
    seconds = swatch.getTimeInSeconds();
    QDPInternal::globalSum(seconds);
    seconds /= (double)Layout::numNodes();

    QDPIO::cout << "." << flush;
  }
  QDPIO::cout << " " << seconds << "(s)"<< endl;

  QDPIO::cout << "Timing expmat with " << hits << " hits " << endl;
  swatch.reset();
  swatch.start();
  for(int i=0; i < hits; i++) { 
   for(int mu=0; mu < Nd; mu++) { 
	expmat(gauge_state.getP()[mu], EXP_EXACT);
      }
  }
  swatch.stop();
  seconds = swatch.getTimeInSeconds();
  QDPInternal::globalSum(seconds);
  seconds /= (double)Layout::numNodes();
    
  QDPIO::cout << "expmat"
	      << ", Nhits = " << hits 
	      << ", Time = " << seconds << "(s)" 
	      << ", Time per hit = " << seconds/(double)hits << "(s)" << endl;

  push(xml_out, "ExpmatTiming");
  write(xml_out, "nhits", hits);
  write(xml_out, "time", seconds);
  write(xml_out, "time_per_hit", (seconds/(double)hits));
  pop(xml_out);

  // Time reunit
  seconds = 0;
  hits = 1;
  int numbad;
  QDPIO::cout << "Timing reunit()" << endl;
  QDPIO::cout << "Calibrating";
  
  while( seconds < 10 ) { 
    hits *= 2;
    swatch.reset();
    swatch.start();
    for(int i=0; i < hits; i++) { 
      for(int mu=0; mu < Nd; mu++) { 
	reunit((gauge_state.getQ())[mu], numbad, REUNITARIZE_ERROR);
      }
    }
    swatch.stop();
    seconds = swatch.getTimeInSeconds();
    QDPInternal::globalSum(seconds);
    seconds /= (double)Layout::numNodes();

    QDPIO::cout << "." << flush;
  }
  QDPIO::cout << " " << seconds << "(s)"<< endl;

  QDPIO::cout << "Timing reunit with " << hits << " hits " << endl;
  swatch.reset();
  swatch.start(); 
  for(int i=0; i < hits; i++) { 
    for(int mu=0; mu < Nd; mu++) { 
     reunit((gauge_state.getQ())[mu], numbad, REUNITARIZE_ERROR);
    }
  }
  swatch.stop();
  seconds = swatch.getTimeInSeconds();
  QDPInternal::globalSum(seconds);
  seconds /= (double)Layout::numNodes();
    
  QDPIO::cout << "reunit" 
	      << ", Nhits = " << hits 
	      << ", Time = " << seconds << "(s)" 
	      << ", Time per hit = " << seconds/(double)hits << "(s)" << endl;

  push(xml_out, "ReunitTiming");
  write(xml_out, "nhits", hits);
  write(xml_out, "time", seconds);
  write(xml_out, "time_per_hit", (seconds/(double)hits));
  pop(xml_out);

#endif

  pop(xml_log);   // t_leapfrog
  pop(xml_out);   // t_leapfrog

  END_CODE();

  // Finish
  Chroma::finalize();
  exit(0);
}
