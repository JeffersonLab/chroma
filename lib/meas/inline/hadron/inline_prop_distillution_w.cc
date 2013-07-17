/*! \file
 * \brief Compute propagators from distillution
 *
 * Propagator calculation in distillution
 */

#include "fermact.h"
#include "meas/inline/hadron/inline_prop_distillution_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "util/ferm/distillution_noise.h"
#include "meas/hadron/distillution_factory.h"
#include "qdp_map_obj.h"
#include "qdp_map_obj_disk.h"
#include "qdp_disk_map_slice.h"
#include "util/ferm/key_prop_distillution.h"
#include "util/ferm/transf.h"
#include "util/ferm/spin_rep.h"
#include "util/ferm/diractodr.h"
#include "util/ferm/twoquark_contract_ops.h"
#include "util/ft/time_slice_set.h"
#include "util/info/proginfo.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"
#include "meas/inline/make_xml_file.h"

#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{ 
  //----------------------------------------------------------------------------
  namespace InlinePropDistillutionEnv 
  {
    //! Propagator input
    void read(XMLReader& xml, const string& path, Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "save_solnP", input.save_solnP);
      read(inputtop, "gauge_id", input.gauge_id);
      read(inputtop, "distillution_id", input.distillution_id);
      read(inputtop, "src_file", input.src_file);
      read(inputtop, "soln_file", input.soln_file);
    }

    //! Propagator output
    void write(XMLWriter& xml, const string& path, const Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "save_solnP", input.save_solnP);
      write(xml, "gauge_id", input.gauge_id);
      write(xml, "distillution_id", input.distillution_id);
      write(xml, "src_file", input.src_file);
      write(xml, "soln_file", input.soln_file);

      pop(xml);
    }


    //! Propagator input
    void read(XMLReader& xml, const string& path, Params::Param_t::Contract_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "quark_lines", input.quark_lines);
      read(inputtop, "mass", input.mass);

      // Read quark-line parameters
      input.quark_line_xml = readXMLGroup(inputtop, "QuarkLine", "QuarkLineType");
    }

    //! Propagator output
    void write(XMLWriter& xml, const string& path, const Params::Param_t::Contract_t& input)
    {
      push(xml, path);

      write(xml, "quark_lines", input.quark_lines);
      write(xml, "mass", input.mass);
      xml << input.quark_line_xml.xml;

      pop(xml);
    }


    //! Propagator input
    void read(XMLReader& xml, const string& path, Params::Param_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "Propagator", input.prop);
      read(inputtop, "Contractions", input.contract);
    }

    //! Propagator output
    void write(XMLWriter& xml, const string& path, const Params::Param_t& input)
    {
      push(xml, path);

      write(xml, "Propagator", input.prop);
      write(xml, "Contractions", input.contract);

      pop(xml);
    }


    //! Propagator input
    void read(XMLReader& xml, const string& path, Params& input)
    {
      Params tmp(xml, path);
      input = tmp;
    }

    //! Propagator output
    void write(XMLWriter& xml, const string& path, const Params& input)
    {
      push(xml, path);
    
      write(xml, "Param", input.param);
      write(xml, "NamedObject", input.named_obj);

      pop(xml);
    }
  } // namespace InlinePropDistillutionEnv 


  //----------------------------------------------------------------------------
  namespace InlinePropDistillutionEnv 
  {
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineMeas(Params(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }
      
    const std::string name = "PROP_DISTILLUTION";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= DistillutionFactoryEnv::registerAll();
	success &= WilsonTypeFermActsEnv::registerAll();
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
    }


    //----------------------------------------------------------------------------
    // Param stuff
    Params::Params() { frequency = 0; }

    Params::Params(XMLReader& xml_in, const std::string& path) 
    {
      try 
      {
	XMLReader paramtop(xml_in, path);

	if (paramtop.count("Frequency") == 1)
	  read(paramtop, "Frequency", frequency);
	else
	  frequency = 1;

	// Parameters for source construction
	read(paramtop, "Param", param);

	// Read in the output propagator/source configuration info
	read(paramtop, "NamedObject", named_obj);

	// Possible alternate XML file pattern
	if (paramtop.count("xml_file") != 0) 
	{
	  read(paramtop, "xml_file", xml_file);
	}
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << endl;
	QDP_abort(1);
      }
    }



    // Function call
    void 
    InlineMeas::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
      // If xml file not empty, then use alternate
      if (params.xml_file != "")
      {
	string xml_file = makeXMLFileName(params.xml_file, update_no);

	push(xml_out, "PropDistillution");
	write(xml_out, "update_no", update_no);
	write(xml_out, "xml_file", xml_file);
	pop(xml_out);

	XMLFileWriter xml(xml_file);
	func(update_no, xml);
      }
      else
      {
	func(update_no, xml_out);
      }
    }


    // Real work done here
    void 
    InlineMeas::func(unsigned long update_no,
		     XMLWriter& xml_out) 
    {
      START_CODE();

      StopWatch snoop;
      snoop.reset();
      snoop.start();

      // Test and grab a reference to the gauge field
      multi1d<LatticeColorMatrix> u;
      XMLBufferWriter gauge_xml;
      try
      {
	u = TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);
	TheNamedObjMap::Instance().get(params.named_obj.gauge_id).getRecordXML(gauge_xml);
      }
      catch( std::bad_cast ) 
      {
	QDPIO::cerr << name << ": caught dynamic cast error" << endl;
	QDP_abort(1);
      }
      catch (const string& e) 
      {
	QDPIO::cerr << name << ": map call failed: " << e << endl;
	QDP_abort(1);
      }

      push(xml_out, "PropDistillution");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << name << ": propagator calculation" << endl;

      proginfo(xml_out);    // Print out basic program info

      // Write out the input
      write(xml_out, "Input", params);

      // Write out the config header
      write(xml_out, "Config_info", gauge_xml);

      push(xml_out, "Output_version");
      write(xml_out, "out_version", 1);
      pop(xml_out);

      // Calculate some gauge invariant observables just for info.
      MesPlq(xml_out, "Observables", u);

      //
      // Read in the source along with relevant information.
      // 
      QDPIO::cout << "Snarf the distillution factory from a named buffer" << endl;
      try
      {
	TheNamedObjMap::Instance().getData< Handle< DistillutionNoise > >(params.named_obj.distillution_id);
      }    
      catch (std::bad_cast) {
	QDPIO::cerr << name << ": caught dynamic cast error" << endl;
	QDP_abort(1);
      }
      catch (const string& e) {
	QDPIO::cerr << name << ": error extracting source_header: " << e << endl;
	QDP_abort(1);
      }
      catch( const char* e) {
	QDPIO::cerr << name <<": Caught some char* exception:" << endl;
	QDPIO::cerr << e << endl;
	QDPIO::cerr << "Rethrowing" << endl;
	throw;
      }

      // Cast should be valid now
      const DistillutionNoise& dist_noise_obj =
	*(TheNamedObjMap::Instance().getData< Handle<DistillutionNoise> >(params.named_obj.distillution_id));

      // Some diagnostics
      QDPIO::cout << "Distillution factory: ensemble= XX" << dist_noise_obj.getEnsemble() << "XX  "
		  << "sequence= XX" << dist_noise_obj.getSequence() << "XX\n"
		  << "t_origin= " << dist_noise_obj.getOrigin() << std::endl;

      // Will use TimeSliceSet-s a lot
      const int decay_dir = dist_noise_obj.getDecayDir();
      const int Lt        = Layout::lattSize()[decay_dir];

      // A sanity check
      if (decay_dir != Nd-1)
      {
	QDPIO::cerr << __func__ << ": TimeSliceIO only supports decay_dir= " << Nd-1 << "\n";
	QDP_abort(1);
      }

      // The time-slice set
      TimeSliceSet time_slice_set(decay_dir);


      //
      // Map-object-disk storage of the source file
      //
      QDP::MapObjectDisk<KeyPropDistillution_t, TimeSliceIO<LatticeColorVectorF> > source_obj;
      source_obj.setDebug(0);

      QDPIO::cout << "Open source file" << endl;

      if (! source_obj.fileExists(params.named_obj.src_file))
      {
	QDPIO::cerr << name << ": source file does not exist: src_file= " << params.named_obj.src_file << std::endl;
	QDP_abort(1);
      }
      else
      {
	source_obj.open(params.named_obj.src_file, std::ios_base::in);
      }

      QDPIO::cout << "Finished opening solution file" << endl;


      //
      // Map-object-disk storage
      //
      QDP::MapObjectDisk<KeyPropDistillution_t, TimeSliceIO<LatticeColorVectorF> > prop_obj;
      prop_obj.setDebug(0);

      if (params.named_obj.save_solnP)
      {
	QDPIO::cout << "Open solution file" << endl;

	if (! prop_obj.fileExists(params.named_obj.soln_file))
	{
	  XMLBufferWriter file_xml;

	  push(file_xml, "MODMetaData");
	  write(file_xml, "id", string("propDist"));
	  write(file_xml, "lattSize", QDP::Layout::lattSize());
	  write(file_xml, "decay_dir", decay_dir);
	  write(file_xml, "ensemble", dist_noise_obj.getEnsemble());
	  write(file_xml, "sequence", dist_noise_obj.getSequence());
	  write(file_xml, "t_origin", dist_noise_obj.getOrigin());
	  file_xml << params.param.contract.quark_line_xml.xml;
	  write(file_xml, "quark_line", params.param.contract.quark_lines[0]);
	  proginfo(file_xml);    // Print out basic program info
	  write(file_xml, "Params", params.param);
	  write(file_xml, "Config_info", gauge_xml);
	  pop(file_xml);

	  std::string file_str(file_xml.str());
	  
	  prop_obj.insertUserdata(file_xml.str());
	  prop_obj.open(params.named_obj.soln_file, std::ios_base::in | std::ios_base::out | std::ios_base::trunc);
	}
	else
	{
	  prop_obj.open(params.named_obj.soln_file);
	}

	QDPIO::cout << "Finished opening solution file" << endl;
      }


      // Total number of iterations
      int ncg_had = 0;


      // Rotation from DR to DP
      SpinMatrix diracToDRMat(DiracToDRMat());
      std::vector<MatrixSpinRep_t> diracToDrMatPlus = convertTwoQuarkSpin(diracToDRMat);
      std::vector<MatrixSpinRep_t> diracToDrMatMinus = convertTwoQuarkSpin(adj(diracToDRMat));


      //
      // Try the factories
      //
      try
      {
	StopWatch swatch;
	swatch.reset();
	QDPIO::cout << "Try the various factories" << endl;

	// Typedefs to save typing
	typedef LatticeFermion               T;
	typedef multi1d<LatticeColorMatrix>  P;
	typedef multi1d<LatticeColorMatrix>  Q;

	//
	// Initialize fermion action
	//
	std::istringstream  xml_s(params.param.prop.fermact.xml);
	XMLReader  fermacttop(xml_s);
	QDPIO::cout << "FermAct = " << params.param.prop.fermact.id << endl;

	// Generic Wilson-Type stuff
	Handle< FermionAction<T,P,Q> >
	  S_f(TheFermionActionFactory::Instance().createObject(params.param.prop.fermact.id,
							       fermacttop,
							       params.param.prop.fermact.path));

	Handle< FermState<T,P,Q> > state(S_f->createState(u));

	Handle< SystemSolver<LatticeFermion> > PP = S_f->qprop(state,
							       params.param.prop.invParam);
      
	QDPIO::cout << "Suitable factory found: compute all the quark props" << endl;
	swatch.start();


	//
	// Loop over the source color and spin, creating the source
	// and calling the relevant propagator routines.
	//
	// Loop over each quark-line
	for(std::vector<int>::const_iterator quark_line= params.param.contract.quark_lines.begin();
	    quark_line != params.param.contract.quark_lines.end();
	    ++quark_line)
	{
	  QDPIO::cout << "quark_line = " << *quark_line << "  type= " << params.param.contract.quark_line_xml.id << endl; 

	  //
	  // Factory for quark line
	  //
	  Handle<AbsQuarkLine> quark_line_fact;

	  try 
	  {
	    std::istringstream  xml_l(params.param.contract.quark_line_xml.xml);
	    XMLReader  linktop(xml_l);
	    QDPIO::cout << "Quark line type = " << params.param.contract.quark_line_xml.id << endl;
	    QDPIO::cout << "Quark line xml = XX" << params.param.contract.quark_line_xml.xml << "XX" << endl;
	    
	    QDPIO::cout << "Create quark-line factory" << std::endl;

	    quark_line_fact = 
	      TheQuarkLineFactory::Instance().createObject(params.param.contract.quark_line_xml.id,
							   linktop,
							   params.param.contract.quark_line_xml.path,
							   dist_noise_obj, 
							   source_obj, 
							   time_slice_set,
							   *quark_line,
							   params.param.contract.mass);

	    QDPIO::cout << "Factory created" << std::endl;
	  }
	  catch(const std::string& e) 
	  {
	    QDPIO::cerr << InlinePropDistillutionEnv::name << ": error creating quark-line factory: " << e << endl;
	    QDP_abort(1);
	  }
	  catch(...) 
	  {
	    QDPIO::cerr << InlinePropDistillutionEnv::name << ": generic exception in creating quark-line factory" << endl;
	    QDP_abort(1);
	  }


	  // Loop over source locations
	  std::vector<int> t_sources(quark_line_fact->getTimeSources());

	  for(int tt=0; tt < t_sources.size(); ++tt)
	  {
	    int t_source = t_sources[tt];  // This is the pretend time-slice. The actual value is shifted.
	    QDPIO::cout << "t_source = " << t_source << endl; 

	    // The space distillution loop
	    for(int dist_src=0; dist_src < quark_line_fact->getNumSpaceDils(); ++dist_src)
	    {
	      StopWatch sniss1;
	      sniss1.reset();
	      sniss1.start();
	      QDPIO::cout << "dist_src = " << dist_src << endl; 

	      // Prepare a distilluted source
	      LatticeColorVector vec_srce = quark_line_fact->getSrc(t_source, dist_src);

	      //
	      // Loop over each spin source and invert. 
	      // Use the same colorvector source. No spin dilution will be used.
	      //
	      multi2d<LatticeColorVector> ferm_out(Ns,Ns);

	      for(int spin_source=0; spin_source < Ns; ++spin_source)
	      {
		QDPIO::cout << "spin_source = " << spin_source << endl; 

		// Insert a ColorVector into spin index spin_source
		// This only overwrites sections, so need to initialize first
		LatticeFermion chi = zero;
		CvToFerm(vec_srce, chi, spin_source);

		LatticeFermion quark_soln = zero;

		// Do the propagator inversion
		SystemSolverResults_t res = (*PP)(quark_soln, chi);
		ncg_had = res.n_count;

		// Extract into the temporary output array
		for(int spin_sink=0; spin_sink < Ns; ++spin_sink)
		{
		  ferm_out(spin_sink,spin_source) = peekSpin(quark_soln, spin_sink);
		}
	      } // for spin_source


	      // Rotate from DeGrand-Rossi (DR) to Dirac-Pauli (DP)
	      {
		multi2d<LatticeColorVector> ferm_tmp;

		multiplyRep(ferm_tmp, diracToDrMatMinus, ferm_out);
		multiplyRep(ferm_out, ferm_tmp, diracToDrMatPlus);
	      }

	      sniss1.stop();
	      QDPIO::cout << "Time to assemble and transmogrify propagators for dist_src= " << dist_src << "  time = " 
			  << sniss1.getTimeInSeconds() 
			  << " secs" << endl;


	      // Write out each time-slice chunk of a lattice colorvec soln to disk
	      QDPIO::cout << "Potentially write propagator solutions to disk" << std::endl;
	      StopWatch sniss2;
	      sniss2.reset();
	      sniss2.start();

	      // Write the solutions
	      if (params.named_obj.save_solnP)
	      {
		QDPIO::cout << "Write propagator solution to disk" << std::endl;
		std::list<KeyPropDistillution_t> snk_keys(quark_line_fact->getSnkKeys(t_source, dist_src));

		for(std::list<KeyPropDistillution_t>::const_iterator key= snk_keys.begin();
		    key != snk_keys.end();
		    ++key)
		{
		  LatticeColorVectorF tmptmp = ferm_out(key->spin_snk,key->spin_src); 

		  prop_obj.insert(*key, TimeSliceIO<LatticeColorVectorF>(tmptmp, dist_noise_obj.getTime(key->t_slice)));
		} // for key
	      }

	      sniss2.stop();
	      QDPIO::cout << "Time to write propagators for dist_src= " << dist_src << "  time = " 
			  << sniss2.getTimeInSeconds() 
			  << " secs" << endl;

	    } // for dist_src
	  } // for tt
	} // for quark_line

	swatch.stop();
	QDPIO::cout << "Propagators computed: time= " 
		    << swatch.getTimeInSeconds() 
		    << " secs" << endl;
      }
      catch (const std::string& e) 
      {
	QDPIO::cout << name << ": caught exception around qprop: " << e << endl;
	QDP_abort(1);
      }

      push(xml_out,"Relaxation_Iterations");
      write(xml_out, "ncg_had", ncg_had);
      pop(xml_out);

      pop(xml_out);  // prop_dist

      snoop.stop();
      QDPIO::cout << name << ": total time = "
		  << snoop.getTimeInSeconds() 
		  << " secs" << endl;

      QDPIO::cout << name << ": ran successfully" << endl;

      END_CODE();
    } 

  }

} // namespace Chroma
