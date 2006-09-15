#include "chroma.h"

#include "io/pg_hmc_io.h"
#include "util/gauge/gaugebc_utils.h"
#include "util/gauge/gaugeact_utils.h"

#include <iostream>

using namespace Chroma;

int main(int argc, char *argv[])
{
  // Initialise QDP
  Chroma::initialize(&argc, &argv);

  XMLReader xml_in(Chroma::getXMLInputFileName());

  struct PureGaugeHMCParams input;

  read(xml_in, "/PureGaugeHMC/params", input);

  XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();
  XMLFileWriter& xml_log = Chroma::getXMLLogInstance();

  push(xml_out, "PureGaugeHMC");
  push(xml_log, "PureGaugeHMC");

  write(xml_out, "params", input);
  write(xml_log, "params", input);

  // Set up the lattice
  Layout::setLattSize(input.MC_iters.nrow);
  Layout::create();

  // Now create the required classes
  
  // Gauge Boundary
  Handle< GaugeBC > gbc( getGaugeBCFromParams( *(input.gauge_bc_handle) ) );


  // Gauge Actions

  // MC Action
  Handle< GaugeAction > S_g_MC_handle( getGaugeActFromParams(gbc,
					     *(input.gauge_act_MC_handle) ) );

  // MD Action
  Handle< GaugeAction > S_g_MD_handle( getGaugeActFromParams(gbc,
					     *(input.gauge_act_MD_handle) ) );


  // Create the Hamiltonians
  ExactPureGaugeHamiltonian H_MC(*(S_g_MC_handle));
  ExactPureGaugeHamiltonian H_MD(*(S_g_MD_handle));

  // Create the MD Integrator
  // Nasty dynamic casting...
  Handle< AbsPureGaugeSympUpdates> leaps_handle;

  Handle< AbsLatColMatHybInt<
      AbsFieldState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >,
      AbsHamiltonian<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
    >
  > integrator_handle;
                              
  // Get the symplectic integrator
  switch( input.MD_params_handle->getType() ) {
  case MD_PQP_LEAPFROG:
    {
      leaps_handle = new PureGaugeSympUpdates(H_MD);
      
      const LeapfrogParams& lf_par =
	dynamic_cast<LeapfrogParams&>(*(input.MD_params_handle));
      
      integrator_handle = new PureGaugePQPLeapFrog(*leaps_handle, 
						   lf_par.getStepSize(),
						   lf_par.getTrajLength());
    
    }
    break;
  default: 
    QDPIO::cerr << "Integrator not yet implemented " << endl;
    QDP_abort(1);
  }

  PureGaugeHMCTraj HMC(H_MC, *integrator_handle);
    
  // Startup the gauge field and momenta
  multi1d<LatticeColorMatrix> u(Nd);
  multi1d<LatticeColorMatrix> p(Nd);

  if( input.MC_startup.restartP == false ) { 
    // Start up from cfg
    XMLReader file_xml;
    XMLReader config_xml;
    
    gaugeStartup(file_xml, config_xml, u, input.MC_startup.start_cfg);
    // setrn(input.MC_startup.seed);
  }
  else { 
    QDPIO::cerr << "Restart from saved state not yet implemented" << endl;
    QDP_abort(1);
  }

  pop(xml_log);
  pop(xml_out);

  Chroma::finalize();
  exit(0);
}

