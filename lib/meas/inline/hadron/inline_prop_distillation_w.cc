/*! \file
 * \brief Compute propagators from distillation
 *
 * Propagator calculation in distillation
 */

#include "fermact.h"
#include "meas/inline/hadron/inline_prop_distillation_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "qdp_map_obj.h"
#include "qdp_map_obj_disk.h"
#include "qdp_disk_map_slice.h"
#include "util/ferm/key_prop_distillation.h"
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
  namespace InlinePropDistillationEnv 
  {
    //! Propagator input
    void read(XMLReader& xml, const string& path, Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_id", input.gauge_id);
      read(inputtop, "src_file", input.src_file);
      read(inputtop, "soln_file", input.soln_file);
    }

    //! Propagator output
    void write(XMLWriter& xml, const string& path, const Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "gauge_id", input.gauge_id);
      write(xml, "src_file", input.src_file);
      write(xml, "soln_file", input.soln_file);

      pop(xml);
    }


    //! Propagator input
    void read(XMLReader& xml, const string& path, Params::Param_t::Contract_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "num_vecs", input.num_vecs);
      read(inputtop, "t_sources", input.t_sources);
      read(inputtop, "decay_dir", input.decay_dir);
      read(inputtop, "Nt_forward", input.Nt_forward);
      read(inputtop, "Nt_backward", input.Nt_backward);
      read(inputtop, "mass_label", input.mass_label);
    }

    //! Propagator output
    void write(XMLWriter& xml, const string& path, const Params::Param_t::Contract_t& input)
    {
      push(xml, path);

      write(xml, "num_vecs", input.num_vecs);
      write(xml, "t_sources", input.t_sources);
      write(xml, "decay_dir", input.decay_dir);
      write(xml, "Nt_forward", input.Nt_forward);
      write(xml, "Nt_backward", input.Nt_backward);
      write(xml, "mass_label", input.mass_label);

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
  } // namespace InlinePropDistillationEnv 


  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  namespace InlinePropDistillationEnv 
  {
    // Anonymous namespace
    namespace
    {
      // Convenience type
      typedef QDP::MapObjectDisk<KeyPropDistillation_t, TimeSliceIO<LatticeColorVectorF> > MOD_t;

      //----------------------------------------------------------------------------
      //! Get source key
      KeyPropDistillation_t getSrcKey(int t_source, int colorvec_src)
      {
	KeyPropDistillation_t key;

	key.t_source     = t_source;
	key.t_slice      = t_source;
	key.colorvec_src = colorvec_src;
	key.spin_src     = -1;
	key.spin_snk     = -1;
	// key.mass   NOTE: no mass key in source

	return key;
      }

	
      //----------------------------------------------------------------------------
      //! Read a source vector
      LatticeColorVector getSrc(MOD_t& source_obj, int t_source, int colorvec_src)
      {
	QDPIO::cout << __func__ << ": on t_source= " << t_source << "  colorvec_src= " << colorvec_src << endl;

	// Get the source vector
	KeyPropDistillation_t src_key = getSrcKey(t_source, colorvec_src);
	LatticeColorVectorF vec_srce = zero;

	TimeSliceIO<LatticeColorVectorF> time_slice_io(vec_srce, t_source);
	source_obj.get(src_key, time_slice_io);

	return vec_srce;
      }


      //----------------------------------------------------------------------------
      //! Get active time-slices
      std::vector<bool> getActiveTSlices(int t_source, int Nt_forward, int Nt_backward)
      {
	// Initialize the active time slices
	const int decay_dir = Nd-1;
	const int Lt = Layout::lattSize()[decay_dir];

	std::vector<bool> active_t_slices(Lt);
	for(int t=0; t < Lt; ++t)
	{
	  active_t_slices[t] = false;
	}

	// Forward
	for(int dt=0; dt < Nt_forward; ++dt)
	{
	  int t = t_source + dt;
	  active_t_slices[t % Lt] = true;
	}

	// Backward
	for(int dt=0; dt < Nt_backward; ++dt)
	{
	  int t = t_source - dt;
	  while (t < 0) {t += Lt;} 

	  active_t_slices[t % Lt] = true;
	}

	return active_t_slices;
      }


      //----------------------------------------------------------------------------
      //! Get source keys
      std::list<KeyPropDistillation_t> getSrcKeys(int t_source, int colorvec_src)
      {
	std::list<KeyPropDistillation_t> keys;

	keys.push_back(getSrcKey(t_source,colorvec_src));

	return keys;
      }

	
      //----------------------------------------------------------------------------
      //! Get sink keys
      std::list<KeyPropDistillation_t> getSnkKeys(int t_source, int colorvec_src, int Nt_forward, int Nt_backward, const std::string mass)
      {
	std::list<KeyPropDistillation_t> keys;

	std::vector<bool> active_t_slices = getActiveTSlices(t_source, Nt_forward, Nt_backward);
	
	const int decay_dir = Nd-1;
	const int Lt = Layout::lattSize()[decay_dir];

	for(int spin_source=0; spin_source < Ns; ++spin_source)
	{
	  for(int spin_sink=0; spin_sink < Ns; ++spin_sink)
	  {
	    for(int t=0; t < Lt; ++t)
	    {
	      if (! active_t_slices[t]) {continue;}

	      KeyPropDistillation_t key;

	      key.t_source     = t_source;
	      key.t_slice      = t;
	      key.colorvec_src = colorvec_src;
	      key.spin_src     = spin_source;
	      key.spin_snk     = spin_sink;
	      key.mass         = mass;

	      //            QDPIO::cout << key << std::flush;

	      keys.push_back(key);
	    } // for t
	  } // for spin_sink
	} // for spin_source

	return keys;
      }
	
    } // end anonymous
  } // end namespace

	
  //----------------------------------------------------------------------------
  namespace InlinePropDistillationEnv 
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
      
    const std::string name = "PROP_DISTILLATION";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
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


    //----------------------------------------------------------------------------
    //----------------------------------------------------------------------------
    // Function call
    void 
    InlineMeas::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
      // If xml file not empty, then use alternate
      if (params.xml_file != "")
      {
	string xml_file = makeXMLFileName(params.xml_file, update_no);

	push(xml_out, "PropDistillation");
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

      push(xml_out, "PropDistillation");
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

      // Will use TimeSliceSet-s a lot
      const int decay_dir = params.param.contract.decay_dir;
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
      QDP::MapObjectDisk<KeyPropDistillation_t, TimeSliceIO<LatticeColorVectorF> > source_obj;
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
      QDP::MapObjectDisk<KeyPropDistillation_t, TimeSliceIO<LatticeColorVectorF> > prop_obj;
      prop_obj.setDebug(0);

      QDPIO::cout << "Open solution file" << endl;

      if (! prop_obj.fileExists(params.named_obj.soln_file))
      {
	XMLBufferWriter file_xml;

	push(file_xml, "MODMetaData");
	write(file_xml, "id", string("propDistillation"));
	write(file_xml, "lattSize", QDP::Layout::lattSize());
	write(file_xml, "decay_dir", decay_dir);
	write(file_xml, "num_vecs", params.param.contract.num_vecs);
	write(file_xml, "Nt_backward", params.param.contract.Nt_backward);
	write(file_xml, "Nt_backward", params.param.contract.Nt_backward);
	write(file_xml, "mass_label", params.param.contract.mass_label);
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


      // Total number of iterations
      int ncg_had = 0;


#if 0
      // NOTE: not doing this spin rotation, but leaving the code here for reference
      // Rotation from DR to DP
      SpinMatrix diracToDRMat(DiracToDRMat());
      std::vector<MatrixSpinRep_t> diracToDrMatPlus = convertTwoQuarkSpin(diracToDRMat);
      std::vector<MatrixSpinRep_t> diracToDrMatMinus = convertTwoQuarkSpin(adj(diracToDRMat));
#endif


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
	const int num_vecs            = params.param.contract.num_vecs;
	const multi1d<int>& t_sources = params.param.contract.t_sources;

	// Loop over each time source
	for(int tt=0; tt < t_sources.size(); ++tt)
	{
	  int t_source = t_sources[tt];  // This is the actual time-slice.
	  QDPIO::cout << "t_source = " << t_source << endl; 

	  // The space distillation loop
	  for(int colorvec_src=0; colorvec_src < num_vecs; ++colorvec_src)
	  {
	    StopWatch sniss1;
	    sniss1.reset();
	    sniss1.start();
	    QDPIO::cout << "colorvec_src = " << colorvec_src << endl; 

	    // Get the source vector
	    LatticeColorVector vec_srce = getSrc(source_obj, t_source, colorvec_src);

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


#if 0
	    // NOTE: not doing this spin rotation, but leaving the code here for reference
	    // Instead, leave the code in the DR basis - the convention for all the perambulators
	    // written

	    // Rotate from DeGrand-Rossi (DR) to Dirac-Pauli (DP)
	    {
	      multi2d<LatticeColorVector> ferm_tmp;

	      multiplyRep(ferm_tmp, diracToDrMatMinus, ferm_out);
	      multiplyRep(ferm_out, ferm_tmp, diracToDrMatPlus);
	    }
#endif

	    sniss1.stop();
	    QDPIO::cout << "Time to assemble and transmogrify propagators for colorvec_src= " << colorvec_src << "  time = " 
			<< sniss1.getTimeInSeconds() 
			<< " secs" << endl;


	    // Write out each time-slice chunk of a lattice colorvec soln to disk
	    QDPIO::cout << "Write propagator solutions to disk" << std::endl;
	    StopWatch sniss2;
	    sniss2.reset();
	    sniss2.start();

	    // Write the solutions
	    QDPIO::cout << "Write propagator solution to disk" << std::endl;
	    std::list<KeyPropDistillation_t> snk_keys(getSnkKeys(t_source, colorvec_src, 
								 params.param.contract.Nt_forward,
								 params.param.contract.Nt_backward,
								 params.param.contract.mass_label));

	    for(std::list<KeyPropDistillation_t>::const_iterator key= snk_keys.begin();
		key != snk_keys.end();
		++key)
	    {
	      LatticeColorVectorF tmptmp = ferm_out(key->spin_snk,key->spin_src); 

	      prop_obj.insert(*key, TimeSliceIO<LatticeColorVectorF>(tmptmp, key->t_slice));
	    } // for key

	    sniss2.stop();
	    QDPIO::cout << "Time to write propagators for colorvec_src= " << colorvec_src << "  time = " 
			<< sniss2.getTimeInSeconds() 
			<< " secs" << endl;
	    
	  } // for colorvec_src
	} // for tt

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
