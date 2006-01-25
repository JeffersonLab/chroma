// $Id: hmc.cc,v 2.3 2006-01-25 05:08:21 edwards Exp $
/*! \file
 *  \brief Main code for HMC with dynamical fermion generation
 */

#include "chroma.h"
#include <string>

using namespace Chroma;

namespace Chroma { 
  
  struct MCControl {
    Cfg_t cfg;
    QDP::Seed rng_seed;
    unsigned long start_update_num;
    unsigned long n_warm_up_updates;
    unsigned long n_production_updates;
    unsigned int  n_updates_this_run;
    unsigned int  save_interval;
    std::string   save_prefix;
    QDP_volfmt_t  save_volfmt;
    std::string   inline_measurement_xml;

  };

  void read(XMLReader& xml, const std::string& path, MCControl& p) {
    try { 
      XMLReader paramtop(xml, path);
      read(paramtop, "./Cfg", p.cfg);
      read(paramtop, "./RNG", p.rng_seed);
      read(paramtop, "./StartUpdateNum", p.start_update_num);
      read(paramtop, "./NWarmUpUpdates", p.n_warm_up_updates);
      read(paramtop, "./NProductionUpdates", p.n_production_updates);
      read(paramtop, "./NUpdatesThisRun", p.n_updates_this_run);
      read(paramtop, "./SaveInterval", p.save_interval);
      read(paramtop, "./SavePrefix", p.save_prefix);
      read(paramtop, "./SaveVolfmt", p.save_volfmt);

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
	QDPIO::cout << "InlineMeasurements are: " << endl;
	QDPIO::cout << p.inline_measurement_xml << endl;
      }
    }
    catch(const std::string& e ) { 
      QDPIO::cerr << "Caught Exception: " << e << endl;
      QDP_abort(1);
    }
  }

  void write(XMLWriter& xml, const std::string& path, const MCControl& p) {
    try {
      push(xml, path);
      write(xml, "Cfg", p.cfg);
      write(xml, "RNG", p.rng_seed);
      write(xml, "StartUpdateNum", p.start_update_num);
      write(xml, "NWarmUpUpdates", p.n_warm_up_updates);
      write(xml, "NProductionUpdates", p.n_production_updates);
      write(xml, "NUpdatesThisRun", p.n_updates_this_run);
      write(xml, "SaveInterval", p.save_interval);
      write(xml, "SavePrefix", p.save_prefix);
      write(xml, "SaveVolfmt", p.save_volfmt);
      xml << p.inline_measurement_xml;

      pop(xml);

    }
    catch(const std::string& e ) { 
      QDPIO::cerr << "Caught Exception: " << e << endl;
      QDP_abort(1);
    }
  }


  struct HMCTrjParams { 

    multi1d<int> nrow;

    // Polymorphic
    std::string H_MC_xml;
    std::string H_MD_xml;
    std::string Integrator_xml;

  };

  void write(XMLWriter& xml, const std::string& path, const HMCTrjParams& p)
  {
    try { 
      push(xml, path);
      write(xml, "nrow", p.nrow);
      xml << p.H_MC_xml;
      xml << p.H_MD_xml;
      xml << p.Integrator_xml;
      pop(xml);
    }
    catch(const std::string& e ) { 
      QDPIO::cerr << "Caught Exception: " << e << endl;
      QDP_abort(1);
    }
  }


  void read(XMLReader& xml, const std::string& path, HMCTrjParams& p) 
  {
    try {
      XMLReader paramtop(xml, path);
      
      read(paramtop, "./nrow", p.nrow);

      // Now the XML for the Hamiltonians
      XMLReader H_MC_xml(paramtop, "./MC_Hamiltonian");
      std::ostringstream os_H_MC;
      H_MC_xml.print(os_H_MC);
      p.H_MC_xml = os_H_MC.str();
      
      QDPIO::cout << "HMC_xml is: " << endl;
      QDPIO::cout << p.H_MC_xml;


      if( paramtop.count("./MD_Hamiltonian") == 1 ) { 
	// Separate MD Hamiltonian specified
	XMLReader H_MD_xml(paramtop, "./MD_Hamiltonian");
	std::ostringstream os_H_MD;
	H_MD_xml.print(os_H_MD);
	p.H_MD_xml = os_H_MD.str();
      }
      else { 
	QDPIO::cout << "Copying MC_Hamiltonian for MD_Hamiltonian" << endl;

	XMLReader H_MC_Monomials(paramtop, "./MC_Hamiltonian/Monomials");
	XMLBufferWriter H_MD_writer;
	push(H_MD_writer, "MD_Hamiltonian");
	write(H_MD_writer, "Monomials", H_MC_Monomials);
	pop(H_MD_writer);
	p.H_MD_xml = H_MD_writer.printCurrentContext();
      }
       
	
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

    // Set the name of the config 
    p_new.cfg.cfg_file = restart_config_filename.str();

    // Hijack this for now and assumes it means what I want it to mean
    p_new.cfg.cfg_type = CFG_TYPE_SZINQIO;


    push(restart_data_buffer, "Params");
    write(restart_data_buffer, "MCControl", p_new);
    write(restart_data_buffer, "HMCTrj", update_params);
    pop(restart_data_buffer);

    // Write a restart DATA file from the buffer XML
    
    XMLFileWriter restart_xml(restart_data_filename.str().c_str());
    restart_xml << restart_data_buffer;
    restart_xml.close();

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
	       QDPIO_SERIAL);    
  }

 

  template<typename UpdateParams>
  void doHMC(multi1d<LatticeColorMatrix>& u,
	     AbsHMCTrj<multi1d<LatticeColorMatrix>,
	               multi1d<LatticeColorMatrix> >& theHMCTrj,
	     MCControl& mc_control, 
	     const UpdateParams& update_params,
	     multi1d< Handle<AbsInlineMeasurement> >& user_measurements) {

    XMLWriter& xml_out = TheXMLOutputWriter::Instance();
    push(xml_out, "doHMC");

    multi1d< Handle< AbsInlineMeasurement > > default_measurements(1);
    InlinePlaquetteParams plaq_params;
    plaq_params.frequency = 1;
    // It is a handle
    default_measurements[0] = new InlinePlaquette(plaq_params);

    try {

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
      if ( total_updates > mc_control.n_updates_this_run + cur_update +1 ) {
	to_do = mc_control.n_updates_this_run;
      }
      else {
	to_do = total_updates - cur_update ;
      }
      
      QDPIO::cout << "MC Control: About to do " << to_do << " updates" << endl;

      // XML Output
      push(xml_out, "MCUpdates");

      for(int i=0; i < to_do; i++) {
	push(xml_out, "elem"); // Caller writes elem rule

	push(xml_out, "Update");
	// Increase current update counter
	cur_update++;
	
	// Decide if the next update is a warm up or not
	bool warm_up_p = cur_update  <= mc_control.n_warm_up_updates;
	QDPIO::cout << "Doing Update: " << cur_update << " warm_up_p = " << warm_up_p << endl;

	// Log
	write(xml_out, "update_no", cur_update);
	write(xml_out, "WarmUpP", warm_up_p);

        QDPIO::cout << "Before HMC trajectory call" << endl;

	// Do the trajectory without accepting 
	theHMCTrj( gauge_state, warm_up_p );

        QDPIO::cout << "After HMC trajectory call" << endl;

	// Create a gauge header for inline measurements.
	// Since there are defaults always measured, we must always
	// create a header.
	//
	// NOTE: THIS HEADER STUFF NEEDS A LOT MORE THOUGHT
	//
	{
	  XMLBufferWriter gauge_xml;
	  push(gauge_xml, "ChromaHMC");
	  write(gauge_xml, "update_no", cur_update);
	  write(gauge_xml, "HMCTrj", update_params);
	  pop(gauge_xml);

	  // Measure inline observables 
	  push(xml_out, "InlineObservables");

	  // Always measure defaults
	  for(int m=0; m < default_measurements.size(); m++) 
	  {
	    // Caller writes elem rule 
	    AbsInlineMeasurement& the_meas = *(default_measurements[m]);
	    push(xml_out, "elem");
	    the_meas( gauge_state.getQ(), gauge_xml, cur_update, xml_out);
	    pop(xml_out);
	  }
	
	  // Only measure user measurements after warm up
	  if( ! warm_up_p ) 
	  {
	    QDPIO::cout << "Doing " << user_measurements.size() 
			<<" user measurements" << endl;
	    for(int m=0; m < user_measurements.size(); m++) 
	    {
	      AbsInlineMeasurement& the_meas = *(user_measurements[m]);
	      if( cur_update % the_meas.getFrequency() == 0 ) 
	      { 
		// Caller writes elem rule
		push(xml_out, "elem");
		the_meas( gauge_state.getQ(), gauge_xml, cur_update, xml_out );
		pop(xml_out); 
	      }
	    }
	  }
	  pop(xml_out); // pop("InlineObservables");
	}

	if( cur_update % mc_control.save_interval == 0 ) {
	  // Save state
	  saveState<UpdateParams>(update_params, mc_control, cur_update, gauge_state.getQ());
	  
	}

	pop(xml_out); // pop("Update");
	pop(xml_out); // pop("elem");
      }   
      
      // Save state
      saveState<UpdateParams>(update_params, mc_control, cur_update, gauge_state.getQ());
      
      pop(xml_out); // pop("MCUpdates")
    }
    catch( const std::string& e) { 
      QDPIO::cerr << "Caught Exception: " << e << endl;
      QDP_abort(1);
    }

    pop(xml_out);
  }
  
  bool linkageHack(void)
  {
    bool foo = true;
    
    // Gauge Monomials
    foo &= GaugeMonomialEnv::registered;
    
    // Ferm Monomials
    foo &= WilsonTypeFermMonomialAggregrateEnv::registered;
    
    // MD Integrators
    foo &= LCMMDIntegratorAggregateEnv::registered;

    // Chrono predictor
    foo &= ChronoPredictorAggregrateEnv::registered;

    // Inline Measurements
    foo &= InlineAggregateEnv::registered;

    return foo;
  }
};

using namespace Chroma;

//! Hybrid Monte Carlo
/*! \defgroup hmcmain Hybrid Monte Carlo
 *  \ingroup main
 *
 * Main program for dynamical fermion generation
 */

int main(int argc, char *argv[]) 
{
  Chroma::initialize(&argc, &argv);
  
  // Chroma Init stuff -- Open DATA and XMLDAT
  linkageHack();

  XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();
  push(xml_out, "hmc");

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
  }
  catch(const std::string& e) {
    QDPIO::cerr << "hmc: Caught Exception while reading file: " << e << endl;
    QDP_abort(1);
  }

  QDPIO::cout << "Call QDP create layout" << endl;
  Layout::setLattSize(trj_params.nrow);
  Layout::create();
  QDPIO::cout << "Finished with QDP create layout" << endl;

  proginfo(xml_out);    // Print out basic program info

  // Start up the config
  multi1d<LatticeColorMatrix> u(Nd);
  {
    XMLReader file_xml;
    XMLReader config_xml;
    
    QDPIO::cout << "Initialize gauge field" << endl;
    gaugeStartup(file_xml, config_xml, u, mc_control.cfg);
    QDPIO::cout << "Finished initializing gauge field" << endl;

    // Write out the config header
    push(xml_out, "Config_info");
    write(xml_out, "file_xml", file_xml);
    write(xml_out, "config_xml", config_xml);
    pop(xml_out);
  }
  
  // Get the MC_Hamiltonian
  std::istringstream H_MC_is(trj_params.H_MC_xml);
  XMLReader H_MC_xml(H_MC_is);
  Handle< ExactAbsHamiltonian< multi1d<LatticeColorMatrix>,     
    multi1d<LatticeColorMatrix> > > H_MC(new ExactLatColMatHamiltonian(H_MC_xml, "/MC_Hamiltonian"));
 

  // Get the MD_Hamiltonian
  std::istringstream H_MD_is(trj_params.H_MD_xml);
  XMLReader H_MD_xml(H_MD_is);
  
  Handle< AbsHamiltonian< multi1d<LatticeColorMatrix>, 
    multi1d<LatticeColorMatrix> > > H_MD(new ExactLatColMatHamiltonian(H_MD_xml, "/MD_Hamiltonian"));


  std::istringstream Integrator_is(trj_params.Integrator_xml);
  XMLReader MD_xml(Integrator_is);
 
    
  // Get the Integrator
  std::string integrator_name;
  read(MD_xml, "/MDIntegrator/Name", integrator_name);
  
  // Get the Leapfrog 
  Handle< AbsMDIntegrator<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > MD( TheMDIntegratorFactory::Instance().createObject(integrator_name, MD_xml, "/MDIntegrator", H_MD) );
  
    
  // Get the HMC
  LatColMatHMCTrj theHMCTrj( H_MC, MD );

 
  multi1d < Handle< AbsInlineMeasurement > > the_measurements;

  // Get the measurements
  try { 
    std::istringstream Measurements_is(mc_control.inline_measurement_xml);

    XMLReader MeasXML(Measurements_is);

    std::ostringstream os;
    MeasXML.print(os);
    QDPIO::cout << os.str() << endl << flush;


    read(MeasXML, "/InlineMeasurements", the_measurements);

  }
  catch(const std::string& e) { 
    QDPIO::cerr << "hmc: Caught exception while reading measurements: " << e << endl
		<< flush;

    QDP_abort(1);
  }

  QDPIO::cout << "There are " << the_measurements.size() << " user measurements " << endl;

  
  // Run
  try { 
    doHMC<HMCTrjParams>(u, theHMCTrj, mc_control, trj_params, the_measurements);
  } 
  catch( const std::string& e ) { 
    QDPIO::cerr << "hmc: Caught string exception: " << e << endl;
    QDP_abort(1);
  }
  catch( std::exception& e ) {
    QDPIO::cerr << "hmc: Caught standard library exception: " << e.what() << endl;
    QDP_abort(1);
  }
  catch(...) {
    QDPIO::cerr << "hmc: Caught generic/unknown exception" << endl;
    QDP_abort(1);
  }

  pop(xml_out);

  Chroma::finalize();
  exit(0);
}

