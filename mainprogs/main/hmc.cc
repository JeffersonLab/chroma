/*! \file
 *  \brief Main code for HMC with dynamical fermion generation
 */

#include "chroma.h"
#include <string>

using namespace Chroma;

namespace Chroma 
{ 
  
  struct MCControl 
  {
    GroupXML_t    cfg;
    QDP::Seed     rng_seed;
    unsigned long start_update_num;
    unsigned long n_warm_up_updates;
    unsigned long n_production_updates;
    unsigned int  n_updates_this_run;
    unsigned int  save_interval;
    std::string   save_prefix;
    QDP_volfmt_t  save_volfmt;
    QDP_serialparallel_t save_pario;
    std::string   inline_measurement_xml;
    bool          repro_checkP;
    int           repro_check_frequency;
    bool          rev_checkP;
    int           rev_check_frequency;
    bool          monitorForcesP;

  };
  
  void read(XMLReader& xml, const std::string& path, MCControl& p) 
  {
    START_CODE();

    try { 
      XMLReader paramtop(xml, path);
      p.cfg = readXMLGroup(paramtop, "Cfg", "cfg_type");
      read(paramtop, "./RNG", p.rng_seed);
      read(paramtop, "./StartUpdateNum", p.start_update_num);
      read(paramtop, "./NWarmUpUpdates", p.n_warm_up_updates);
      read(paramtop, "./NProductionUpdates", p.n_production_updates);
      read(paramtop, "./NUpdatesThisRun", p.n_updates_this_run);
      read(paramtop, "./SaveInterval", p.save_interval);
      read(paramtop, "./SavePrefix", p.save_prefix);
      read(paramtop, "./SaveVolfmt", p.save_volfmt);

      // -- Deal with parallel IO
      p.save_pario = QDPIO_SERIAL; // Default

      // If there is a ParalelIO tag 
      if ( paramtop.count("./ParallelIO") > 0 ) {

	bool parioP=false;
	read(paramtop, "./ParallelIO", parioP);
	if ( parioP ) {
	  QDPIO::cout << "Setting parallel write mode" << std::endl;
	  p.save_pario = QDPIO_PARALLEL;
	}
      }

      // Default values: repro check is on, frequency is 10%
      p.repro_checkP = true;
      p.repro_check_frequency = 10;


      // Now overwrite with user values
      if( paramtop.count("./ReproCheckP") == 1 ) {
	// Read user value if given
	read(paramtop, "./ReproCheckP", p.repro_checkP);
      }

      // If we do repro checking, read the frequency
      if( p.repro_checkP ) { 
	if( paramtop.count("./ReproCheckFrequency") == 1 ) {
	  // Read user value if given
	  read(paramtop, "./ReproCheckFrequency", p.repro_check_frequency);
	}
      }

      // Reversibility checking enabled by default.
      p.rev_checkP = true;
      p.rev_check_frequency = 10;

      // Now overwrite with user values
      if( paramtop.count("./ReverseCheckP") == 1 ) {
	// Read user value if given
	read(paramtop, "./ReverseCheckP", p.rev_checkP);
      }

      // If we do repro checking, read the frequency
      if( p.rev_checkP ) { 
	if( paramtop.count("./ReverseCheckFrequency") == 1 ) {
	  // Read user value if given
	  read(paramtop, "./ReverseCheckFrequency", p.rev_check_frequency);
	}
      }

      if( paramtop.count("./MonitorForces") == 1 ) {
	read(paramtop, "./MonitorForces", p.monitorForcesP);
      }
      else { 
	p.monitorForcesP = true;
      }

      if( paramtop.count("./InlineMeasurements") == 0 ) {
	XMLBufferWriter dummy;
	push(dummy, "InlineMeasurements");
	pop(dummy); // InlineMeasurements
	p.inline_measurement_xml = dummy.printCurrentContext();
	
      }
      else {
	XMLReader measurements_xml(paramtop, "./InlineMeasurements");
	std::ostringstream inline_os;
	measurements_xml.print(inline_os);
	p.inline_measurement_xml = inline_os.str();
	QDPIO::cout << "InlineMeasurements are: " << std::endl;
	QDPIO::cout << p.inline_measurement_xml << std::endl;
      }

      
    }
    catch(const std::string& e ) { 
      QDPIO::cerr << "Caught Exception: " << e << std::endl;
      QDP_abort(1);
    }
    
    END_CODE();
  }

  void write(XMLWriter& xml, const std::string& path, const MCControl& p) 
  {
    START_CODE();
    
    try {
      push(xml, path);
      xml << p.cfg.xml;
      write(xml, "RNG", p.rng_seed);
      write(xml, "StartUpdateNum", p.start_update_num);
      write(xml, "NWarmUpUpdates", p.n_warm_up_updates);
      write(xml, "NProductionUpdates", p.n_production_updates);
      write(xml, "NUpdatesThisRun", p.n_updates_this_run);
      write(xml, "SaveInterval", p.save_interval);
      write(xml, "SavePrefix", p.save_prefix);
      write(xml, "SaveVolfmt", p.save_volfmt);
      { 
	bool pario = ( p.save_pario == QDPIO_PARALLEL );
	write(xml, "ParallelIO", pario);
      }
      write(xml, "ReproCheckP", p.repro_checkP);
      if( p.repro_checkP ) { 
	write(xml, "ReproCheckFrequency", p.repro_check_frequency);
      }
      write(xml, "ReverseCheckP", p.rev_checkP);
      if( p.rev_checkP ) { 
	write(xml, "ReverseCheckFrequency", p.rev_check_frequency);
      }
      write(xml, "MonitorForces", p.monitorForcesP);

      xml << p.inline_measurement_xml;
      
      pop(xml);
      
    }
    catch(const std::string& e ) { 
      QDPIO::cerr << "Caught Exception: " << e << std::endl;
      QDP_abort(1);
    }
    
    END_CODE();
  }


  struct HMCTrjParams 
  { 
    multi1d<int> nrow;
    
    // Polymorphic
    std::string Monomials_xml; // XML for the monomials
    std::string H_MC_xml;      // XML for the Hamiltonian
    std::string Integrator_xml; // XML for the Integrator
  };
  
  void write(XMLWriter& xml, const std::string& path, const HMCTrjParams& p)
  {
    START_CODE();
    
    try { 
      push(xml, path);
      write(xml, "nrow", p.nrow);
      xml << p.Monomials_xml;   // XML For the mon
      xml << p.H_MC_xml;
      xml << p.Integrator_xml;
      pop(xml);
    }
    catch(const std::string& e ) { 
      QDPIO::cerr << "Caught Exception: " << e << std::endl;
      QDP_abort(1);
    }
    
    END_CODE();
  }


  void read(XMLReader& xml, const std::string& path, HMCTrjParams& p) 
  {
    START_CODE();
    
    try {
      XMLReader paramtop(xml, path);
      
      read(paramtop, "./nrow", p.nrow);
      XMLReader Monomials_xml_reader(paramtop, "./Monomials");
      std::ostringstream os_Monomials;
      Monomials_xml_reader.print(os_Monomials);
      p.Monomials_xml = os_Monomials.str();
      QDPIO::cout << "Monomials xml is:" << std::endl;
      QDPIO::cout << p.Monomials_xml << std::endl;

      // Now the XML for the Hamiltonians
      XMLReader H_MC_xml(paramtop, "./Hamiltonian");
      std::ostringstream os_H_MC;
      H_MC_xml.print(os_H_MC);
      p.H_MC_xml = os_H_MC.str();
      
      QDPIO::cout << "HMC_xml is: " << std::endl;
      QDPIO::cout << p.H_MC_xml;
      
      
      // Read the Integrator parameters
      XMLReader MD_integrator_xml(paramtop, "./MDIntegrator");
      std::ostringstream os_integrator;
      MD_integrator_xml.print(os_integrator);
      p.Integrator_xml = os_integrator.str();

      QDPIO::cout << "Integrator XML is: " << std::endl;
      QDPIO::cout << p.Integrator_xml << std::endl;
    }
    catch( const std::string& e ) { 
      QDPIO::cerr << "Error reading XML : " << e << std::endl;
      QDP_abort(1);
    }
    
    END_CODE();
  }

  template<typename UpdateParams>
  void saveState(const UpdateParams& update_params, 
		 MCControl& mc_control,
		 unsigned long update_no,
		const multi1d<LatticeColorMatrix>& u) {
    // Do nothing
  }

  // Specialise
  template<>
  void saveState(const HMCTrjParams& update_params, 
		 MCControl& mc_control,
		 unsigned long update_no,
		 const multi1d<LatticeColorMatrix>& u)
  {
    START_CODE();
    
    // File names
    std::ostringstream restart_data_filename;
    restart_data_filename << mc_control.save_prefix << "_restart_" << update_no << ".xml" ;
    
    std::ostringstream restart_config_filename;
    restart_config_filename << mc_control.save_prefix << "_cfg_" << update_no << ".lime";
      
    XMLBufferWriter restart_data_buffer;

    
    // Copy old params
    MCControl p_new = mc_control;
    
    // Get Current RNG Seed
    QDP::RNG::savern(p_new.rng_seed);
   
    // Set the current traj number
    p_new.start_update_num = update_no;
    
    // Set the num_updates_this_run
    unsigned long total = mc_control.n_warm_up_updates 
      + mc_control.n_production_updates ;

    if ( total < mc_control.n_updates_this_run + update_no ) { 
      p_new.n_updates_this_run = total - update_no;
    }

    // Set the name and type of the config 
    {
      // Parse the cfg XML including the parallel IO part
      SZINQIOGaugeInitEnv::Params  cfg;

      // Reset the filename in it
      cfg.cfg_file = restart_config_filename.str();
      cfg.cfg_pario = mc_control.save_pario;

      // Prepare to write out
      p_new.cfg = SZINQIOGaugeInitEnv::createXMLGroup(cfg);
    }


    push(restart_data_buffer, "Params");
    write(restart_data_buffer, "MCControl", p_new);
    write(restart_data_buffer, "HMCTrj", update_params);
    pop(restart_data_buffer);


    // Save the config

    // some dummy header for the file
    XMLBufferWriter file_xml;
    push(file_xml, "HMC");
    proginfo(file_xml);
    pop(file_xml);


    // Save the config
    writeGauge(file_xml, 
	       restart_data_buffer,
	       u,
	       restart_config_filename.str(),
	       p_new.save_volfmt,
	       p_new.save_pario);    


    // Write a restart DATA file from the buffer XML 
    // Do this after the config is written, so that if the cfg
    // write fails, there is no restart file...
    //
    // production will then likely fall back to last good pair.

    XMLFileWriter restart_xml(restart_data_filename.str().c_str());
    restart_xml << restart_data_buffer;
    restart_xml.close();
    
    END_CODE();
  }


 
  // Predeclare this 
  bool checkReproducability( const multi1d<LatticeColorMatrix>& P_new, 
			     const multi1d<LatticeColorMatrix>& Q_new,
			     const QDP::Seed& seed_new,
			     const multi1d<LatticeColorMatrix>& P_old,
			     const multi1d<LatticeColorMatrix>& Q_old,
			     const QDP::Seed& seed_old );

  template<typename UpdateParams>
  void doHMC(multi1d<LatticeColorMatrix>& u,
	     AbsHMCTrj<multi1d<LatticeColorMatrix>,
	               multi1d<LatticeColorMatrix> >& theHMCTrj,
	     MCControl& mc_control, 
	     const UpdateParams& update_params,
	     multi1d< Handle<AbsInlineMeasurement> >& user_measurements) 
  {
    START_CODE();


    // Turn monitoring off/on
    QDPIO::cout << "Setting Force monitoring to " << mc_control.monitorForcesP  << std::endl;
    setForceMonitoring(mc_control.monitorForcesP) ;
    QDP::StopWatch swatch;

    XMLWriter& xml_out = TheXMLOutputWriter::Instance();
    XMLWriter& xml_log = TheXMLLogWriter::Instance();

    push(xml_out, "doHMC");
    push(xml_log, "doHMC");

    multi1d< Handle< AbsInlineMeasurement > > default_measurements(1);
    InlinePlaquetteEnv::Params plaq_params;
    plaq_params.frequency = 1;
    // It is a handle
    default_measurements[0] = new InlinePlaquetteEnv::InlineMeas(plaq_params);

    {
      // Initialise the RNG
      QDP::RNG::setrn(mc_control.rng_seed);
      
      // Fictitious momenta for now
      multi1d<LatticeColorMatrix> p(Nd);
      
      // Create a field state
      GaugeFieldState gauge_state(p,u);
      
      // Set the update number
      unsigned long cur_update=mc_control.start_update_num;
      
      // Compute how many updates to do
      unsigned long total_updates = mc_control.n_warm_up_updates
	+ mc_control.n_production_updates;
      
      unsigned long to_do = 0;
      if ( total_updates >= mc_control.n_updates_this_run + cur_update +1 ) {
	to_do = mc_control.n_updates_this_run;
      }
      else {
	to_do = total_updates - cur_update ;
      }
      
      QDPIO::cout << "MC Control: About to do " << to_do << " updates" << std::endl;

      // XML Output
      push(xml_out, "MCUpdates");
      push(xml_log, "MCUpdates");

      for(int i=0; i < to_do; i++) 
      {
	push(xml_out, "elem"); // Caller writes elem rule
	push(xml_log, "elem"); // Caller writes elem rule

	push(xml_out, "Update");
	push(xml_log, "Update");
	// Increase current update counter
	cur_update++;
	
	// Decide if the next update is a warm up or not
	bool warm_up_p = cur_update  <= mc_control.n_warm_up_updates;
	QDPIO::cout << "Doing Update: " << cur_update << " warm_up_p = " << warm_up_p << std::endl;

	// Log
	write(xml_out, "update_no", cur_update);
	write(xml_log, "update_no", cur_update);

	write(xml_out, "WarmUpP", warm_up_p);
	write(xml_log, "WarmUpP", warm_up_p);

	bool do_reverse = false;
	if( mc_control.rev_checkP 
	    && ( cur_update % mc_control.rev_check_frequency == 0 )) {
	  do_reverse = true;
	  QDPIO::cout << "Doing Reversibility Test this traj" << std::endl;
	}


	// Check if I need to do any reproducibility testing
	if( mc_control.repro_checkP 
	    && (cur_update % mc_control.repro_check_frequency == 0 ) 
	    ) { 

	  // Yes - reproducibility trajectory
	  // Save fields and RNG at start of trajectory
	  QDPIO::cout << "Saving start config and RNG seed for reproducability test" << std::endl;

	  GaugeFieldState repro_bkup_start( gauge_state.getP(), gauge_state.getQ());
	  QDP::Seed rng_seed_bkup_start;
	  QDP::RNG::savern(rng_seed_bkup_start);
	  
	  // DO the trajectory
	  QDPIO::cout << "Before HMC trajectory call" << std::endl;
	  swatch.reset(); 
	  swatch.start();

	  // This may do a reversibility check 
	  theHMCTrj( gauge_state, warm_up_p, do_reverse ); 
	  swatch.stop(); 
	  
	  QDPIO::cout << "After HMC trajectory call: time= "
		      << swatch.getTimeInSeconds() 
		      << " secs" << std::endl;
	  
	  write(xml_out, "seconds_for_trajectory", swatch.getTimeInSeconds());
	  write(xml_log, "seconds_for_trajectory", swatch.getTimeInSeconds());

	  // Save the fields and RNG at the end
	  QDPIO::cout << "Saving end config and RNG seed for reproducability test" << std::endl;
	  GaugeFieldState repro_bkup_end( gauge_state.getP(), gauge_state.getQ());
	  QDP::Seed rng_seed_bkup_end;
	  QDP::RNG::savern(rng_seed_bkup_end);

	  // Restore the starting field and the starting RNG
	  QDPIO::cout << "Restoring start config and RNG for reproducability test" << std::endl;

	  gauge_state.getP() = repro_bkup_start.getP(); 
	  gauge_state.getQ() = repro_bkup_start.getQ(); 
	  QDP::RNG::setrn(rng_seed_bkup_start); 

	  // Do the repro trajectory
	  QDPIO::cout << "Before HMC repro trajectory call" << std::endl;
	  swatch.reset(); 
	  swatch.start();
	  // Dont repeat the reversibility check in the repro test
	  theHMCTrj( gauge_state, warm_up_p, false ); 
	  swatch.stop(); 
	  
	  QDPIO::cout << "After HMC repro trajectory call: time= "
		      << swatch.getTimeInSeconds() 
		      << " secs" << std::endl;
	  
	  write(xml_out, "seconds_for_repro_trajectory", swatch.getTimeInSeconds());
	  write(xml_log, "seconds_for_repro_trajectory", swatch.getTimeInSeconds());
 
	  // Save seed at end of traj for comparison
	  QDP::Seed rng_seed_end2;
	  QDP::RNG::savern(rng_seed_end2);


	  // Check the reproducibility 
	  bool pass = checkReproducability( gauge_state.getP(), 
					    gauge_state.getQ(), 
					    rng_seed_end2,
					    repro_bkup_end.getP(), 
					    repro_bkup_end.getQ(), 
					    rng_seed_bkup_end);

	  
	  if( !pass ) { 
	    QDPIO::cout << "Reproducability check failed on update " << cur_update << std::endl;
	    QDPIO::cout << "Aborting" << std::endl;
	    write(xml_out, "ReproCheck", pass);
	    write(xml_log, "ReproCheck", pass);
	    QDP_abort(1);
	  }
	  else { 
	    QDPIO::cout << "Reproducability check passed on update " << cur_update << std::endl;
	    write(xml_out, "ReproCheck", pass);
	    write(xml_log, "ReproCheck", pass);
	  }


	}
	else { 

	  // Do the trajectory without accepting
	  QDPIO::cout << "Before HMC trajectory call" << std::endl;
	  swatch.reset();
	  swatch.start();
	  theHMCTrj( gauge_state, warm_up_p, do_reverse  );
	  swatch.stop();
	
	  QDPIO::cout << "After HMC trajectory call: time= "
		      << swatch.getTimeInSeconds() 
		      << " secs" << std::endl;
	  
	  write(xml_out, "seconds_for_trajectory", swatch.getTimeInSeconds());
	  write(xml_log, "seconds_for_trajectory", swatch.getTimeInSeconds());

	}
	swatch.reset();
	swatch.start();

	// Create a gauge header for inline measurements.
	// Since there are defaults always measured, we must always
	// create a header.
	//
	// NOTE: THIS HEADER STUFF NEEDS A LOT MORE THOUGHT
	//
	QDPIO::cout << "HMC: start inline measurements" << std::endl;
	{
	  XMLBufferWriter gauge_xml;
	  push(gauge_xml, "ChromaHMC");
	  write(gauge_xml, "update_no", cur_update);
	  write(gauge_xml, "HMCTrj", update_params);
	  pop(gauge_xml);

	  // Reset and set the default gauge field
	  QDPIO::cout << "HMC: initial resetting default gauge field" << std::endl;
	  InlineDefaultGaugeField::reset();
	  QDPIO::cout << "HMC: set default gauge field" << std::endl;
	  InlineDefaultGaugeField::set(gauge_state.getQ(), gauge_xml);
	  QDPIO::cout << "HMC: finished setting default gauge field" << std::endl;

	  // Measure inline observables 
	  push(xml_out, "InlineObservables");

	  // Always measure defaults
	  for(int m=0; m < default_measurements.size(); m++) 
	  {
	    QDPIO::cout << "HMC: do default measurement = " << m << std::endl;
	    QDPIO::cout << "HMC: dump named objects" << std::endl;
	    TheNamedObjMap::Instance().dump();

	    // Caller writes elem rule 
	    AbsInlineMeasurement& the_meas = *(default_measurements[m]);
	    push(xml_out, "elem");
	    the_meas(cur_update, xml_out);
	    pop(xml_out);

	    QDPIO::cout << "HMC: finished default measurement = " << m << std::endl;
	  }
	
	  // Only measure user measurements after warm up
	  if( ! warm_up_p ) 
	  {
	    QDPIO::cout << "Doing " << user_measurements.size() 
			<<" user measurements" << std::endl;
	    for(int m=0; m < user_measurements.size(); m++) 
	    {
	      QDPIO::cout << "HMC: considering user measurement number = " << m << std::endl;
	      AbsInlineMeasurement& the_meas = *(user_measurements[m]);
	      if( cur_update % the_meas.getFrequency() == 0 ) 
	      { 
		// Caller writes elem rule
		push(xml_out, "elem");
		QDPIO::cout << "HMC: calling user measurement number = " << m << std::endl;
		the_meas(cur_update, xml_out);
		QDPIO::cout << "HMC: finished user measurement number = " << m << std::endl;
		pop(xml_out); 
	      }
	    }
	    QDPIO::cout << "HMC: finished user measurements" << std::endl;
	  }
	  pop(xml_out); // pop("InlineObservables");

	  // Reset the default gauge field
	  QDPIO::cout << "HMC: final resetting default gauge field" << std::endl;
	  InlineDefaultGaugeField::reset();
	  QDPIO::cout << "HMC: finished final resetting default gauge field" << std::endl;
	}

	swatch.stop();
	QDPIO::cout << "After all measurements: time= "
		    << swatch.getTimeInSeconds() 
		    << " secs" << std::endl;

	write(xml_out, "seconds_for_measurements", swatch.getTimeInSeconds());
	write(xml_log, "seconds_for_measurements", swatch.getTimeInSeconds());

	if( cur_update % mc_control.save_interval == 0 ) 
	{
	  swatch.reset();
	  swatch.start();

	  // Save state
	  saveState<UpdateParams>(update_params, mc_control, cur_update, gauge_state.getQ());

	  swatch.stop();
	  QDPIO::cout << "After saving state: time= "
		      << swatch.getTimeInSeconds() 
		      << " secs" << std::endl;
	}

	pop(xml_log); // pop("Update");
	pop(xml_out); // pop("Update");

	pop(xml_log); // pop("elem");
	pop(xml_out); // pop("elem");
      }   
      
      // Save state
      saveState<UpdateParams>(update_params, mc_control, cur_update, gauge_state.getQ());
      
      pop(xml_log); // pop("MCUpdates")
      pop(xml_out); // pop("MCUpdates")
    }

    pop(xml_log); // pop("doHMC")
    pop(xml_out); // pop("doHMC")
    
    END_CODE();
  }
  
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

    // Gauge initialization
    foo &= GaugeInitEnv::registerAll();

    return foo;
  }
};

using namespace Chroma;

//! Hybrid Monte Carlo
/*! \defgroup hmcmain Hybrid Monte Carlo
 *  \ingroup main
 *
 * Main program for dynamical fermion generation
xml */

int main(int argc, char *argv[]) 
{
  Chroma::initialize(&argc, &argv);
  
  START_CODE();

  // Chroma Init stuff -- Open DATA and XMLDAT
  QDPIO::cout << "Linkage = " << linkageHack() << std::endl;

  StopWatch snoop;
  snoop.reset();
  snoop.start();

  XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();
  XMLFileWriter& xml_log = Chroma::getXMLLogInstance();

  push(xml_out, "hmc");
  push(xml_log, "hmc");

  HMCTrjParams trj_params;
  MCControl    mc_control;

  try
  {
    XMLReader xml_in(Chroma::getXMLInputFileName());

    XMLReader paramtop(xml_in, "/Params");
    read( paramtop, "./HMCTrj", trj_params);
    read( paramtop, "./MCControl", mc_control);
    
    // Write out the input
    write(xml_out, "Input", xml_in);
    write(xml_log, "Input", xml_in);
  }
  catch(const std::string& e) {
    QDPIO::cerr << "hmc: Caught Exception while reading file: " << e << std::endl;
    QDP_abort(1);
  }

  if (mc_control.start_update_num >= mc_control.n_production_updates)
  {
    QDPIO::cout << "hmc: run is finished" << std::endl;
    pop(xml_log);
    pop(xml_out);
    exit(0);
  }

  QDPIO::cout << "Call QDP create layout" << std::endl;
  Layout::setLattSize(trj_params.nrow);
  Layout::create();
  QDPIO::cout << "Finished with QDP create layout" << std::endl;

  proginfo(xml_out);    // Print out basic program info
  proginfo(xml_log);    // Print out basic program info

  // Start up the config
  multi1d<LatticeColorMatrix> u(Nd);
  try
  {
    XMLReader file_xml;
    XMLReader config_xml;
    
    QDPIO::cout << "Initialize gauge field" << std::endl;
    StopWatch swatch;
    swatch.reset();
    swatch.start();
    {
      std::istringstream  xml_c(mc_control.cfg.xml);
      XMLReader  cfgtop(xml_c);
      QDPIO::cout << "Gauge initialization: cfg_type = " << mc_control.cfg.id << std::endl;

      Handle< GaugeInit >
	gaugeInit(TheGaugeInitFactory::Instance().createObject(mc_control.cfg.id,
							       cfgtop,
							       mc_control.cfg.path));
      (*gaugeInit)(file_xml, config_xml, u);
    }
    swatch.stop();
    QDPIO::cout << "Gauge field successfully initialized: time= " 
		<< swatch.getTimeInSeconds() 
		<< " secs" << std::endl;

    swatch.reset();
    swatch.start();
    {
      for(int mu=0; mu < Nd; mu++) { 
	reunit(u[mu]);
      }
    }
    swatch.stop();
    QDPIO::cout << "Gauge field reunitarized: time="
		<< swatch.getTimeInSeconds()
		<< " secs" << std::endl;

    // Write out the config header
    write(xml_out, "Config_info", config_xml);
    write(xml_log, "Config_info", config_xml);
  }
  catch(std::bad_cast) 
  {
    QDPIO::cerr << "hmc: caught cast error" << std::endl;
    QDP_abort(1);
  }
  catch(const std::string& e) 
  {
    QDPIO::cerr << "hmc: Caught Exception: " << e << std::endl;
    QDP_abort(1);
  }
  catch(std::exception& e) 
  {
    QDPIO::cerr << "hmc: Caught standard library exception: " << e.what() << std::endl;
    QDP_abort(1);
  }
  catch(...)
  {
    QDPIO::cerr << "hmc: caught generic exception during measurement" << std::endl;
    QDP_abort(1);
  }

  
  // Set up the monomials
  try { 
    std::istringstream Monomial_is(trj_params.Monomials_xml);
    XMLReader monomial_reader(Monomial_is);
    readNamedMonomialArray(monomial_reader, "/Monomials");
  }
  catch(const std::string& e) { 
    QDPIO::cout << "Caught Exception reading Monomials" << std::endl;
    QDP_abort(1);
  }

  std::istringstream H_MC_is(trj_params.H_MC_xml);
  XMLReader H_MC_xml(H_MC_is);
  ExactHamiltonianParams ham_params(H_MC_xml, "/Hamiltonian");
    
  Handle< AbsHamiltonian< multi1d<LatticeColorMatrix>,     
    multi1d<LatticeColorMatrix> > > H_MC(new ExactHamiltonian(ham_params));
 

  std::istringstream MDInt_is(trj_params.Integrator_xml);
  XMLReader MDInt_xml(MDInt_is);
  LCMToplevelIntegratorParams int_par(MDInt_xml, "/MDIntegrator");
  Handle< AbsMDIntegrator< multi1d<LatticeColorMatrix>,
    multi1d<LatticeColorMatrix> > > Integrator(new LCMToplevelIntegrator(int_par));


  LatColMatHMCTrj theHMCTrj( H_MC, Integrator );

 
  multi1d < Handle< AbsInlineMeasurement > > the_measurements;

  // Get the measurements
  try 
  { 
    std::istringstream Measurements_is(mc_control.inline_measurement_xml);

    XMLReader MeasXML(Measurements_is);

    std::ostringstream os;
    MeasXML.print(os);
    QDPIO::cout << os.str() << std::endl << std::flush;

    read(MeasXML, "/InlineMeasurements", the_measurements);
  }
  catch(const std::string& e) { 
    QDPIO::cerr << "hmc: Caught exception while reading measurements: " << e << std::endl
		<< std::flush;

    QDP_abort(1);
  }

  QDPIO::cout << "There are " << the_measurements.size() << " user measurements " << std::endl;

  
  // Run
  try { 
    doHMC<HMCTrjParams>(u, theHMCTrj, mc_control, trj_params, the_measurements);
  } 
  catch(std::bad_cast) 
  {
    QDPIO::cerr << "HMC: caught cast error" << std::endl;
    QDP_abort(1);
  }
  catch(std::bad_alloc) 
  { 
    // This might happen on any node, so report it
    std::cerr << "HMC: caught bad memory allocation" << std::endl;
    QDP_abort(1);
  }
  catch(const std::string& e) 
  { 
    QDPIO::cerr << "HMC: Caught std::string exception: " << e << std::endl;
    QDP_abort(1);
  }
  catch(std::exception& e) 
  {
    QDPIO::cerr << "HMC: Caught standard library exception: " << e.what() << std::endl;
    QDP_abort(1);
  }
  catch(...) 
  {
    QDPIO::cerr << "HMC: Caught generic/unknown exception" << std::endl;
    QDP_abort(1);
  }

  pop(xml_log);  // hmc
  pop(xml_out);  // hmc

  snoop.stop();
  QDPIO::cout << "HMC: total time = "
	      << snoop.getTimeInSeconds() 
	      << " secs" << std::endl;

  END_CODE();

  Chroma::finalize();
  exit(0);
}


// Repro check
namespace Chroma { 

  bool 
  checkReproducability( const multi1d<LatticeColorMatrix>& P_new, 
			const multi1d<LatticeColorMatrix>& Q_new,
			const QDP::Seed& seed_new,
			const multi1d<LatticeColorMatrix>& P_old,
			const multi1d<LatticeColorMatrix>& Q_old,
			const QDP::Seed& seed_old ) 
  {
    
    // Compare local contributions of P
    int diffs_found = 0;
    if ( P_new.size() != P_old.size() ) { 
      // Something bad happened if P_old and P_new are not the same
      return false;
    }
    
    if ( Q_new.size() != Q_old.size() ) { 
      // Something bad happened if Q_old and Q_new are not the same
      return false;
    }
    
    // Count the number of bytes to compare
    int bytes=
      2*Nc*Nc*Layout::sitesOnNode()*sizeof(WordType< LatticeColorMatrix >::Type_t);
    
    // Compare P_s
    for(int mu=0; mu < P_new.size(); mu++) { 
      const unsigned char *p1 = (const unsigned char *)P_new[mu].getF();
      const unsigned char *p2 = (const unsigned char *)P_old[mu].getF();
      for(int b=0; b < bytes; b++) {
	unsigned char diff = *p1 - *p2; 
	if( diff != 0 ) diffs_found++;
	p1++; p2++;
      }
    }
    
    // Sum up the number of diffs found globally
    QDPInternal::globalSum(diffs_found);
    
    if( diffs_found != 0 ) { 
      QDPIO::cout << "Found " << diffs_found << " different bytes in momentum repro check" << std::endl;
      return false;
    }
    
    diffs_found = 0;
    for(int mu=0; mu < P_new.size(); mu++) { 
      const unsigned char *p1 = (const unsigned char *)Q_new[mu].getF();
      const unsigned char *p2 = (const unsigned char *)Q_old[mu].getF();
      for(int b=0; b < bytes; b++) {
	unsigned char diff = *p1 - *p2; 
	if( diff != 0 ) diffs_found++;
	p1++; p2++;
      }
    }
    
    // Sum up the number of diffs found globally
    QDPInternal::globalSum(diffs_found);
    
    if( diffs_found != 0 ) { 
      QDPIO::cout << "Found " << diffs_found << " different bytes in gauge repro check" << std::endl;
      return false;
    }
  
    if( ! toBool( seed_new == seed_old ) ) { 
      QDPIO::cout << "New and old RNG seeds do not match " << std::endl;
      return false;
    }
    
    return true;
  }

}
