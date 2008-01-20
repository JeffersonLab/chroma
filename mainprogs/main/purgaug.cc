// $Id: purgaug.cc,v 3.10 2008-01-20 17:46:43 edwards Exp $
/*! \file
 *  \brief Main code for pure gauge field generation
 */

#include "chroma.h"
#include "actions/gauge/gaugeacts/gaugeacts_aggregate.h"

using namespace Chroma;

namespace Chroma 
{

  //! Holds gauge action
  struct HBGauge
  {
    string  gauge_act;     /*!<  Holds gauge action xml */
  };


  //! Read the parameters
  void read(XMLReader& xml_in, const std::string& path, HBGauge& p)
  {
    try {
      // Read the inverter Parameters
      XMLReader xml_tmp(xml_in, "./GaugeAction");
      std::ostringstream os;
      xml_tmp.print(os);
      p.gauge_act = os.str();
    }
    catch(const string& s) {
      QDPIO::cerr << "Caught Exception while reading gauge action: " << s <<endl;
      QDP_abort(1);
    }

    QDPIO::cout << "Gauge action: read \n" << p.gauge_act << endl;
  }


  //! Writer
  void write(XMLWriter& xml, const std::string& path, const HBGauge& p)
  {
    xml << p.gauge_act;
  }


  //! Reader
  void read(XMLReader& xml, const std::string& path, HBParams& p)
  {
    try { 
      XMLReader paramtop(xml, path);
      read(paramtop, "NmaxHB", p.NmaxHB);
      read(paramtop, "nOver", p.nOver);
    }
    catch(const std::string& e ) { 
      QDPIO::cerr << "Caught Exception reading HBParams: " << e << endl;
      QDP_abort(1);
    }
  }

  //! Writer
  void write(XMLWriter& xml, const std::string& path, const HBParams& p)
  {
    push(xml, path);

    write(xml, "NmaxHB", p.NmaxHB);
    write(xml, "nOver", p.nOver);

    pop(xml);
  }


  //! Params controlling running of monte carlo
  struct MCControl 
  {
    QDP::Seed rng_seed;
    unsigned long start_update_num;
    unsigned long n_warm_up_updates;
    unsigned long n_production_updates;
    unsigned int  n_updates_this_run;
    unsigned int  save_interval;
    std::string   save_prefix;
    QDP_volfmt_t  save_volfmt;
  };

  void read(XMLReader& xml, const std::string& path, MCControl& p) 
  {
    try { 
      XMLReader paramtop(xml, path);
      read(paramtop, "./RNG", p.rng_seed);
      read(paramtop, "./StartUpdateNum", p.start_update_num);
      read(paramtop, "./NWarmUpUpdates", p.n_warm_up_updates);
      read(paramtop, "./NProductionUpdates", p.n_production_updates);
      read(paramtop, "./NUpdatesThisRun", p.n_updates_this_run);
      read(paramtop, "./SaveInterval", p.save_interval);
      read(paramtop, "./SavePrefix", p.save_prefix);
      read(paramtop, "./SaveVolfmt", p.save_volfmt);

      if (p.n_updates_this_run % p.save_interval != 0)
	throw string("UpdateThisRun not a multiple of SaveInterval");
    }
    catch(const std::string& e ) { 
      QDPIO::cerr << "Caught Exception reading MCControl: " << e << endl;
      QDP_abort(1);
    }
  }

  void write(XMLWriter& xml, const std::string& path, const MCControl& p) 
  {
    push(xml, path);

    write(xml, "RNG", p.rng_seed);
    write(xml, "StartUpdateNum", p.start_update_num);
    write(xml, "NWarmUpUpdates", p.n_warm_up_updates);
    write(xml, "NProductionUpdates", p.n_production_updates);
    write(xml, "NUpdatesThisRun", p.n_updates_this_run);
    write(xml, "SaveInterval", p.save_interval);
    write(xml, "SavePrefix", p.save_prefix);
    write(xml, "SaveVolfmt", p.save_volfmt);

    pop(xml);
  }


  //! Holds params for Heat-bath
  struct HBItrParams 
  { 
    multi1d<int> nrow;

    HBGauge   hb_gaugeact;    /*!< This is polymorphic */
    HBParams  hb_params;      /*!< Solely the HB bit */
  };

  void write(XMLWriter& xml, const std::string& path, const HBItrParams& p)
  {
    push(xml, path);
    write(xml, "nrow", p.nrow);
    write(xml, "GaugeAction", p.hb_gaugeact);
    write(xml, "HBParams", p.hb_params);
    pop(xml);
  }


  void read(XMLReader& xml, const std::string& path, HBItrParams& p) 
  {
    try {
      XMLReader paramtop(xml, path);
      
      read(paramtop, "nrow", p.nrow);
      read(paramtop, "GaugeAction", p.hb_gaugeact);
      read(paramtop, "HBParams", p.hb_params);
    }
    catch( const std::string& e ) { 
      QDPIO::cerr << "Error reading HBItrParams XML : " << e << endl;
      QDP_abort(1);
    }
  }

  //! Main struct from input params and output restarts
  struct HBControl 
  {
    HBItrParams   hbitr_params;
    MCControl     mc_control;
    Cfg_t         cfg;
    std::string   inline_measurement_xml;
  };


  //! Reader
  void read(XMLReader& xml_in, const std::string& path, HBControl& p) 
  {
    try {
      XMLReader paramtop(xml_in, path);

      read(paramtop, "HBItr", p.hbitr_params);
      read(paramtop, "MCControl", p.mc_control);
      read(paramtop, "Cfg", p.cfg);

      if( paramtop.count("./InlineMeasurements") == 0 ) {
	XMLBufferWriter dummy;
	push(dummy, "InlineMeasurements");
	pop(dummy); // InlineMeasurements
	p.inline_measurement_xml = dummy.printCurrentContext();
      }
      else 
      {
	XMLReader measurements_xml(paramtop, "./InlineMeasurements");
	std::ostringstream inline_os;
	measurements_xml.print(inline_os);
	p.inline_measurement_xml = inline_os.str();
	QDPIO::cout << "InlineMeasurements are: " << endl;
	QDPIO::cout << p.inline_measurement_xml << endl;
      }
    }
    catch(const std::string& e) {
      QDPIO::cerr << "Caught Exception reading HBControl: " << e << endl;
      QDP_abort(1);
    }
  }


  //! Writer
  void write(XMLWriter& xml, const std::string& path, const HBControl& p) 
  {
    push(xml, path);

    write(xml, "Cfg", p.cfg);
    write(xml, "MCControl", p.mc_control);
    xml << p.inline_measurement_xml;
    write(xml, "HBItr", p.hbitr_params);

    pop(xml);
  }



  //--------------------------------------------------------------------------
  // Specialise
  MCControl newMCHeader(const HBItrParams& update_params, 
			const MCControl& mc_control,
			unsigned long update_no)
  {
    START_CODE();

    // Copy old params
    MCControl p_new = mc_control;
    
    // Get Current RNG Seed
    QDP::RNG::savern(p_new.rng_seed);
   
    // Set the current traj number
    p_new.start_update_num = update_no;
    
    // Reset the warmups
    p_new.n_warm_up_updates = 0;
    
    // Set the num_updates_this_run
    unsigned long total = mc_control.n_production_updates;

    if ( total < mc_control.n_updates_this_run + update_no ) { 
      p_new.n_updates_this_run = total - update_no;
    }

    END_CODE();

    return p_new;
  }



  // Specialise
  void saveState(const HBItrParams& update_params, 
		 MCControl& mc_control,
		 unsigned long update_no,
		 const string& inline_measurement_xml,
		 const multi1d<LatticeColorMatrix>& u)
  {
    START_CODE();

    MCControl mc_new = newMCHeader(update_params, mc_control, update_no);

    // Files
    std::ostringstream restart_data_filename;
    std::ostringstream restart_config_filename;

    unsigned long save_num = update_no / mc_control.save_interval;
    restart_data_filename << mc_control.save_prefix << ".ini.xml" << save_num;
    restart_config_filename << mc_control.save_prefix << ".lime" << save_num;
    
    {
      HBControl hb;
      hb.hbitr_params = update_params;
      hb.mc_control = mc_new;
      hb.inline_measurement_xml = inline_measurement_xml;

      // Set the name of the restart file
      hb.cfg.cfg_file = restart_config_filename.str();

      // Hijack this for now and assumes it means what I want it to mean
      hb.cfg.cfg_type = CFG_TYPE_SZINQIO;

      // Write a restart DATA file from the buffer XML
      XMLFileWriter restart_xml(restart_data_filename.str().c_str());
      write(restart_xml, "purgaug", hb);
      restart_xml.close();
    }

    {
      // Save the config

      // some dummy header for the file
      XMLBufferWriter file_xml;
      push(file_xml, "HB");
      proginfo(file_xml);
      pop(file_xml);

      XMLBufferWriter config_xml;
      push(config_xml, "ChromaHB");
      write(config_xml, "MCControl", mc_new);
      write(config_xml, "HBItr", update_params);
      pop(config_xml);

      // Save the config
      writeGauge(file_xml, 
		 config_xml,
		 u,
		 restart_config_filename.str(),
		 mc_new.save_volfmt,
		 QDPIO_SERIAL);    
    }
    
    END_CODE();
  }


  //--------------------------------------------------------------------------
  void doMeas(XMLWriter& xml_out,
	      multi1d<LatticeColorMatrix>& u,
	      HBControl& hb_control, 
	      bool warm_up_p,
	      unsigned long cur_update,
	      const multi1d< Handle< AbsInlineMeasurement > >& default_measurements,
	      const multi1d< Handle<AbsInlineMeasurement> >& user_measurements) 
  {
    START_CODE();

    // Create a gauge header for inline measurements.
    // Since there are defaults always measured, we must always
    // create a header.
    //
    // NOTE: THIS HEADER STUFF NEEDS A LOT MORE THOUGHT
    //
    MCControl mc_new = newMCHeader(hb_control.hbitr_params, hb_control.mc_control, cur_update);

    XMLBufferWriter gauge_xml;
    push(gauge_xml, "ChromaHB");
    write(gauge_xml, "MCControl", mc_new);
    write(gauge_xml, "HBItr", hb_control.hbitr_params);
    pop(gauge_xml);

    // Reset and set the default gauge field
    InlineDefaultGaugeField::reset();
    InlineDefaultGaugeField::set(u, gauge_xml);

    // Measure inline observables 
    push(xml_out, "InlineObservables");

    // Always measure defaults
    for(int m=0; m < default_measurements.size(); m++) 
    {
      // Caller writes elem rule 
      AbsInlineMeasurement& the_meas = *(default_measurements[m]);
      push(xml_out, "elem");
      the_meas(cur_update, xml_out);
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
	  the_meas(cur_update, xml_out );
	  pop(xml_out); 
	}
      }
    }
    pop(xml_out); // pop("InlineObservables");

    // Reset the default gauge field
    InlineDefaultGaugeField::reset();
    
    END_CODE();
  }
  


  //--------------------------------------------------------------------------
  void doWarmUp(XMLWriter& xml_out,
		multi1d<LatticeColorMatrix>& u,
		const LinearGaugeAction& S_g,
		HBControl& hb_control,
		const multi1d< Handle< AbsInlineMeasurement > >& default_measurements,
		const multi1d< Handle<AbsInlineMeasurement> >& user_measurements) 
  {
    START_CODE();

    // Set the update number
    unsigned long cur_update = 0;
      
    // Compute how many updates to do
    unsigned long to_do = hb_control.mc_control.n_warm_up_updates;
      
    QDPIO::cout << "WarmUp Control: About to do " << to_do << " updates" << endl;

    // XML Output
    push(xml_out, "WarmUpdates");

    for(int i=0; i < to_do; i++) 
    {
      push(xml_out, "elem"); // Caller writes elem rule

      push(xml_out, "Update");
      // Increase current update counter
      cur_update++;
	
      // Log
      write(xml_out, "update_no", cur_update);
      write(xml_out, "WarmUpP", true);

      // Do the update, but with no measurements
      mciter(u, S_g, hb_control.hbitr_params.hb_params); //one hb sweep

      // Do measurements
      doMeas(xml_out, u, hb_control, true, cur_update,
	     default_measurements, user_measurements);

      pop(xml_out); // pop("Update");
      pop(xml_out); // pop("elem");
    }

    pop(xml_out); // pop("WarmUpdates")
    
    END_CODE();
  }
  

  //--------------------------------------------------------------------------
  void doProd(XMLWriter& xml_out,
	      multi1d<LatticeColorMatrix>& u,
	      const LinearGaugeAction& S_g,
	      HBControl& hb_control, 
	      const multi1d< Handle< AbsInlineMeasurement > >& default_measurements,
	      const multi1d< Handle<AbsInlineMeasurement> >& user_measurements) 
  {
    START_CODE();

    // Set the update number
    unsigned long cur_update = hb_control.mc_control.start_update_num;
      
    // Compute how many updates to do
    unsigned long total_updates = hb_control.mc_control.n_production_updates;
      
    unsigned long to_do = 0;
    if ( total_updates > hb_control.mc_control.n_updates_this_run + cur_update +1 ) {
      to_do = hb_control.mc_control.n_updates_this_run;
    }
    else {
      to_do = total_updates - cur_update ;
    }
      
    QDPIO::cout << "MC Control: About to do " << to_do << " updates" << endl;

    // XML Output
    push(xml_out, "MCUpdates");

    for(int i=0; i < to_do; i++) 
    {
      push(xml_out, "elem"); // Caller writes elem rule

      push(xml_out, "Update");
      // Increase current update counter
      cur_update++;
	
      // Decide if the next update is a warm up or not
      QDPIO::cout << "Doing Update: " << cur_update << " warm_up_p = " << false << endl;

      // Log
      write(xml_out, "update_no", cur_update);
      write(xml_out, "WarmUpP", false);

      // Do the update
      mciter(u, S_g, hb_control.hbitr_params.hb_params); //one hb sweep

      // Do measurements
      doMeas(xml_out, u, hb_control, false, cur_update,
	     default_measurements, user_measurements);

      // Save if needed
      if( cur_update % hb_control.mc_control.save_interval == 0 ) 
      {
	saveState(hb_control.hbitr_params, hb_control.mc_control, 
		  cur_update,
		  hb_control.inline_measurement_xml, u);
      }

      pop(xml_out); // pop("Update");
      pop(xml_out); // pop("elem");
    }

    pop(xml_out); // pop("MCUpdates")
    
    END_CODE();
  }
  

  //--------------------------------------------------------------------------
  void doHB(multi1d<LatticeColorMatrix>& u,
	    const LinearGaugeAction& S_g,
	    HBControl& hb_control, 
	    multi1d< Handle<AbsInlineMeasurement> >& user_measurements) 
  {
    START_CODE();

    XMLWriter& xml_out = TheXMLOutputWriter::Instance();
    push(xml_out, "doHB");

    multi1d< Handle< AbsInlineMeasurement > > default_measurements(1);
    InlinePlaquetteEnv::Params plaq_params;
    plaq_params.frequency = 1;
    // It is a handle
    default_measurements[0] = new InlinePlaquetteEnv::InlineMeas(plaq_params);

    try 
    {
      // Initialise the RNG
      QDP::RNG::setrn(hb_control.mc_control.rng_seed);
      
      // If warmups are required, do them first
      if (hb_control.mc_control.n_warm_up_updates > 0)
      {
	doWarmUp(xml_out, u, S_g, hb_control, default_measurements, user_measurements);
	hb_control.mc_control.n_warm_up_updates = 0;  // reset
      }

      // Do the production updates
      doProd(xml_out, u, S_g, hb_control, default_measurements, user_measurements);
    }
    catch( const std::string& e) { 
      QDPIO::cerr << "Caught Exception: " << e << endl;
      QDP_abort(1);
    }

    pop(xml_out);
    
    END_CODE();
  }
  
  bool linkageHack(void)
  {
    bool foo = true;
    
    // Gauge actions
    foo &= GaugeActsEnv::registerAll();

    // Inline Measurements
    foo &= InlineAggregateEnv::registerAll();

    return foo;
  }
}


using namespace Chroma;

//! Pure gauge field generation via heatbath
/*! \defgroup purgaug Heat-bath
 *  \ingroup main
 *
 * Main program for heat-bath generation of 
 */

int main(int argc, char *argv[]) 
{
  Chroma::initialize(&argc, &argv);
  
  START_CODE();

  // Chroma Init stuff -- Open DATA and XMLDAT
  linkageHack();

  XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();
  push(xml_out, "purgaug");

  HBControl  hb_control;

  try
  {
    XMLReader xml_in(Chroma::getXMLInputFileName());

    read(xml_in, "/purgaug", hb_control);

    // Write out the input
    write(xml_out, "Input", xml_in);
  }
  catch( const std::string& e ) {
    QDPIO::cerr << "Caught Exception reading input XML: " << e << endl;
    QDP_abort(1);
  }
  catch( std::exception& e ) {
    QDPIO::cerr << "Caught standard library exception: " << e.what() << endl;
    QDP_abort(1);
  }
  catch(...) {
    QDPIO::cerr << "Caught unknown exception " << endl;
    QDP_abort(1);
  }

  Layout::setLattSize(hb_control.hbitr_params.nrow);
  Layout::create();

  proginfo(xml_out);    // Print out basic program info

  // Start up the config
  multi1d<LatticeColorMatrix> u(Nd);
  {
    XMLReader file_xml;
    XMLReader config_xml;
    
    gaugeStartup(file_xml, config_xml, u, hb_control.cfg);

    // Write out the config header
    push(xml_out, "Config_info");
    write(xml_out, "file_xml", file_xml);
    write(xml_out, "config_xml", config_xml);
    pop(xml_out);
  }
  

  // Create the gauge action
  // This code is limited to only rb sets and subsets.
  // The main point is
  // The number of subsets within the staples, etc. and to get the subsets
  // straight
  Handle< LinearGaugeAction > S_g;
  try
  {
    std::istringstream is(hb_control.hbitr_params.hb_gaugeact.gauge_act);
    XMLReader gaugeact_reader(is);

    // Get the name of the gauge act
    std::string gaugeact_string;
    try { 
      read(gaugeact_reader, "/GaugeAction/Name", gaugeact_string);
    }
    catch( const std::string& e) 
    {
      QDPIO::cerr << "Error grepping the gaugeact name: " << e<<  endl;
      QDP_abort(1);
    }

    // Throw an exception if not found
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    GaugeAction<P,Q>* gaugeact = 
      TheGaugeActFactory::Instance().createObject(gaugeact_string, 
						  gaugeact_reader, 
						  "/GaugeAction");
    S_g = dynamic_cast<LinearGaugeAction*>(gaugeact);
  }
  catch(std::bad_cast) 
  {
    QDPIO::cerr << "PURGAUG: caught cast error" << endl;
    QDP_abort(1);
  }
  catch(std::bad_alloc) 
  { 
    // This might happen on any node, so report it
    cerr << "PURGAUG: caught bad memory allocation" << endl;
    QDP_abort(1);
  }
  catch(const std::string& e) 
  {
    QDPIO::cerr << "PURGAUG: Caught Exception: " << e << endl;
    QDP_abort(1);
  }
  catch(std::exception& e) 
  {
    QDPIO::cerr << "PURGAUG: Caught standard library exception: " << e.what() << endl;
    QDP_abort(1);
  }
  catch(...)
  {
    // This might happen on any node, so report it
    cerr << "PURGAUG: caught generic exception during measurement" << endl;
    QDP_abort(1);
  }


  // Get the measurements
  multi1d < Handle< AbsInlineMeasurement > > the_measurements;

  try { 
    std::istringstream Measurements_is(hb_control.inline_measurement_xml);

    XMLReader MeasXML(Measurements_is);

    std::ostringstream os;
    MeasXML.print(os);
    QDPIO::cout << os.str() << endl << flush;

    read(MeasXML, "/InlineMeasurements", the_measurements);
  }
  catch(const std::string& e) { 
    QDPIO::cerr << "Caugth exception while reading measurements: " << e << endl
		<< flush;

    QDP_abort(1);
  }

  QDPIO::cout << "There are " << the_measurements.size() << " user measurements " << endl;

  
  // Run
  try 
  {
    doHB(u, *S_g, hb_control, the_measurements);
  } 
  catch(std::bad_cast) 
  {
    QDPIO::cerr << "PURGAUG: caught cast error" << endl;
    QDP_abort(1);
  }
  catch(std::bad_alloc) 
  { 
    QDPIO::cerr << "PURGAUG: caught bad memory allocation" << endl;
    QDP_abort(1);
  }
  catch(const std::string& e) 
  { 
    QDPIO::cerr << "PURGAUG: Caught string exception: " << e << endl;
    QDP_abort(1);
  }
  catch(std::exception& e) 
  {
    QDPIO::cerr << "PURGAUG: Caught standard library exception: " << e.what() << endl;
    QDP_abort(1);
  }
  catch(...) 
  {
    QDPIO::cerr << "PURGAUG: Caught generic/unknown exception" << endl;
    QDP_abort(1);
  }

  pop(xml_out);

  END_CODE();

  Chroma::finalize();
  exit(0);
}

