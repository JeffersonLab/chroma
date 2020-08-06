/*! \file
 * \brief Compute propagators from distillation
 *
 * Propagator calculation in distillation
 */

#include "fermact.h"
#include "meas/inline/hadron/inline_prop_distillation_stoch_w.h"
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
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"
#include "meas/inline/make_xml_file.h"

#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{ 
  //----------------------------------------------------------------------------------
  // Utility functions
  namespace
  {
    //! Error output
    StandardOutputStream& operator<<(StandardOutputStream& os, const multi1d<int>& d)
    {
      if (d.size() > 0)
      {
	os << d[0];

	for(int i=1; i < d.size(); ++i)
	  os << " " << d[i];
      }

      return os;
    }
  }

  //----------------------------------------------------------------------------
  namespace InlinePropDistillationStochEnv 
  {
    //! Propagator input
    void read(XMLReader& xml, const std::string& path, Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_id", input.gauge_id);
      read(inputtop, "src_file", input.src_file);
      read(inputtop, "soln_file", input.soln_file);
      read(inputtop, "prop_file", input.prop_file);
    }

    //! Propagator output
    void write(XMLWriter& xml, const std::string& path, const Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "gauge_id", input.gauge_id);
      write(xml, "src_file", input.src_file);
      write(xml, "soln_file", input.soln_file);
      write(xml, "prop_file", input.prop_file);

      pop(xml);
    }


    //! Propagator input
    void read(XMLReader& xml, const std::string& path, Params::Param_t::Contract_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "num_vecs", input.num_vecs);
      read(inputtop, "t_sources", input.t_sources);
      read(inputtop, "decay_dir", input.decay_dir);
      read(inputtop, "Nt_forward", input.Nt_forward);
      read(inputtop, "Nt_backward", input.Nt_backward);
      read(inputtop, "mass_label", input.mass_label);
      read(inputtop, "spatial_mask_size", input.spatial_mask_size);
    }

    //! Propagator output
    void write(XMLWriter& xml, const std::string& path, const Params::Param_t::Contract_t& input)
    {
      push(xml, path);

      write(xml, "num_vecs", input.num_vecs);
      write(xml, "t_sources", input.t_sources);
      write(xml, "decay_dir", input.decay_dir);
      write(xml, "Nt_forward", input.Nt_forward);
      write(xml, "Nt_backward", input.Nt_backward);
      write(xml, "mass_label", input.mass_label);
      write(xml, "spatial_mask_size", input.spatial_mask_size);

      pop(xml);
    }


    //! Propagator input
    void read(XMLReader& xml, const std::string& path, Params::Param_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "Propagator", input.prop);
      read(inputtop, "Contractions", input.contract);
    }

    //! Propagator output
    void write(XMLWriter& xml, const std::string& path, const Params::Param_t& input)
    {
      push(xml, path);

      write(xml, "Propagator", input.prop);
      write(xml, "Contractions", input.contract);

      pop(xml);
    }


    //! Propagator input
    void read(XMLReader& xml, const std::string& path, Params& input)
    {
      Params tmp(xml, path);
      input = tmp;
    }

    //! Propagator output
    void write(XMLWriter& xml, const std::string& path, const Params& input)
    {
      push(xml, path);
    
      write(xml, "Param", input.param);
      write(xml, "NamedObject", input.named_obj);

      pop(xml);
    }
  } // namespace InlinePropDistillationStochEnv 


  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  namespace InlinePropDistillationStochEnv 
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
      //! Read a source std::vector
      LatticeColorVector getSrc(MOD_t& source_obj, int t_source, int colorvec_src)
      {
	QDPIO::cout << __func__ << ": on t_source= " << t_source << "  colorvec_src= " << colorvec_src << std::endl;

	// Get the source std::vector
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
	
      //----------------------------------------------------------------------------
      //! Get soln
      multi2d<LatticeColorVector> getSoln(MOD_t& source_obj, int t_source, int colorvec_src, const std::string mass)
      {
	multi2d<LatticeColorVector> ferm_out(Ns,Ns);

	const int decay_dir = Nd-1;
	const int Lt = Layout::lattSize()[decay_dir];

	for(int spin_source=0; spin_source < Ns; ++spin_source)
	{
	  for(int spin_sink=0; spin_sink < Ns; ++spin_sink)
	  {
	    KeyPropDistillation_t key;

	    key.t_source     = t_source;
	    key.t_slice      = t_source;
	    key.colorvec_src = colorvec_src;
	    key.spin_src     = spin_source;
	    key.spin_snk     = spin_sink;
	    key.mass         = mass;

	    // Get the source std::vector
	    LatticeColorVectorF vec_srce = zero;

	    TimeSliceIO<LatticeColorVectorF> time_slice_io(vec_srce, t_source);
	    if (source_obj.get(key, time_slice_io) != 0)
	    {
	      QDPIO::cerr << __func__ << ": error retrieving t_source= " << t_source << std::endl;
	      QDP_abort(1);
	    }

	    ferm_out(spin_sink, spin_source) = vec_srce;
	  } // for spin_sink
	} // for spin_source

	return ferm_out;
      }


      //----------------------------------------------------------------------------
      //! Get soln
      void doTrace(multi2d<ComplexD>& trace_mom,
		   const LatticeColorVector& vec_srce,
		   const multi2d<LatticeColorVector>& fred_out, int t_source, int colorvec_src,
		   const std::vector<MatrixSpinRep_t>& diracToDrMatPlus,
		   const std::vector<MatrixSpinRep_t>& diracToDrMatMinus,
		   const SftMom& phases,
		   int num_vecs)
      {
	// Check
	if (1)
	{
	  Double ddd = norm2(vec_srce);   // check for debugging
	  QDPIO::cout << __func__ << ": colorvec_src= " << colorvec_src << "  norm(left_vec)= " << ddd << "\n";
	}

	// Rotate from DeGrand-Rossi (DR) to Dirac-Pauli (DP)
	multi2d<LatticeColorVector> ferm_out;
	{
	  multi2d<LatticeColorVector> ferm_tmp;
	    
	  multiplyRep(ferm_tmp, diracToDrMatMinus, fred_out);
	  multiplyRep(ferm_out, ferm_tmp, diracToDrMatPlus);
	}
	
	for(int spin_source = 0; spin_source < Ns; ++spin_source)
	{
	  LatticeComplex lop = localInnerProduct(vec_srce, ferm_out(spin_source,spin_source));
	  
	  // Slow fourier-transform
	  for(int mom_num = 0; mom_num < phases.numMom(); ++mom_num)
	  {
	    multi1d<ComplexD> op_sum = sumMulti(phases[mom_num] * lop, phases.getSet());

#if 0
	    for(int ttt = 0; ttt < op_sum.size(); ++ttt)
	    {
	      QDPIO::cout << "OP_SUM: "
			  << "  t_source= " << ttt
			  << "  colorvec_src= " << colorvec_src
			  << "  spin_source= " << spin_source
			  << "  mom_num= " << mom_num
			  << "  mom= " << phases.numToMom(mom_num)
			  << "  op_sum[ " << ttt << " ]= " << op_sum[ttt]
			  << std::endl;
	    }
#endif

	    trace_mom(t_source, mom_num) += op_sum[t_source];
	  }
	}
	  
	for(int mom_num = 0 ; mom_num < phases.numMom() ; ++mom_num)
	{
	  QDPIO::cout << "TRACE: "
		      << "  t_source= " << t_source
		      << "  colorvec_src= " << colorvec_src
		      << "  num_vecs= " << num_vecs
		      << "  mom= " << phases.numToMom(mom_num)
		      << "  trace( " << mom_num << " )= " << trace_mom(t_source,mom_num)
		      << std::endl;
	}
      }

    } // end anonymous
  } // end namespace

	
  //----------------------------------------------------------------------------
  namespace InlinePropDistillationStochEnv 
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
      
    const std::string name = "PROP_DISTILLATION_STOCH";

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
	QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << std::endl;
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
	std::string xml_file = makeXMLFileName(params.xml_file, update_no);

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
	QDPIO::cerr << name << ": caught dynamic cast error" << std::endl;
	QDP_abort(1);
      }
      catch (const std::string& e) 
      {
	QDPIO::cerr << name << ": std::map call failed: " << e << std::endl;
	QDP_abort(1);
      }

      push(xml_out, "PropDistillation");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << name << ": propagator calculation" << std::endl;

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
      // Sanity checks
      //
      if (params.param.contract.spatial_mask_size.size() != Nd-1)
      {
	QDPIO::cerr << name << ": param.contract.spatial mask size incorrect 1" << std::endl;
	QDP_abort(1);
      }
      
      for(int mu=0; mu < Nd-1; ++mu)
      {
	if (params.param.contract.spatial_mask_size[mu] == 0)
	{
	  QDPIO::cerr << name << ": not valid param.contract.spatial mask" << std::endl;
	  QDP_abort(1);
	}

	if ((QDP::Layout::lattSize()[mu] % params.param.contract.spatial_mask_size[mu]) != 0)
	{
	  QDPIO::cerr << name << ": lattice size not divisible by param.contract.spatial mask size" << std::endl;
	  QDP_abort(1);
	}
      }
      

      //
      // Build a list of all the colorvec sources
      //
      std::list<int> source_list;

      if (1)
      {
#if 0
	// The distillation sources
	for(int colorvec_src = 0; colorvec_src < num_vecs; ++colorvec_src)
	  source_list.push_back(colorvec_src);
#endif

	// The number of spatial noises
	int num_sites = 1;
	for(int mu = 0; mu < params.param.contract.spatial_mask_size.size(); ++mu)
	{
	  num_sites *= params.param.contract.spatial_mask_size[mu];
	}

	// Loop over all sites with a spatial mask
	// These will form the sources
#if 0
        for(int ipos=1; ipos <= num_sites; ++ipos)
#endif
	for(int ipos = num_sites; ipos > 0; --ipos)
	  source_list.push_back(-ipos);
      }


      //
      // Map-object-disk storage of the source file
      //
      QDP::MapObjectDisk<KeyPropDistillation_t, TimeSliceIO<LatticeColorVectorF> > source_obj;
      source_obj.setDebug(0);

      QDPIO::cout << "Open source file" << std::endl;

      if (! source_obj.fileExists(params.named_obj.src_file))
      {
	QDPIO::cerr << name << ": source file does not exist: src_file= " << params.named_obj.src_file << std::endl;
	QDP_abort(1);
      }
      else
      {
	source_obj.open(params.named_obj.src_file, std::ios_base::in);
      }

      QDPIO::cout << "Finished opening solution file" << std::endl;


      //
      // Map-object-disk storage
      //
      QDP::MapObjectDisk<KeyPropDistillation_t, TimeSliceIO<LatticeColorVectorF> > soln_obj;
      soln_obj.setDebug(0);

      QDPIO::cout << "Open solution file" << std::endl;

      if (! soln_obj.fileExists(params.named_obj.soln_file))
      {
	QDPIO::cerr << name << ": soln file does not exist: soln_file= " << params.named_obj.soln_file << std::endl;
	QDP_abort(1);
      }
      else
      {
	soln_obj.open(params.named_obj.soln_file, std::ios_base::in);
      }
      
      QDPIO::cout << "Finished opening solution file" << std::endl;


      //
      // Map-object-disk storage
      //
      QDP::MapObjectDisk<KeyPropDistillation_t, TimeSliceIO<LatticeColorVectorF> > prop_obj;
      prop_obj.setDebug(0);

      QDPIO::cout << "Open solution file" << std::endl;

      if (! prop_obj.fileExists(params.named_obj.prop_file))
      {
	XMLBufferWriter file_xml;

	push(file_xml, "MODMetaData");
	write(file_xml, "id", std::string("propDistillation"));
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
	prop_obj.open(params.named_obj.prop_file, std::ios_base::in | std::ios_base::out | std::ios_base::trunc);
      }
      else
      {
	prop_obj.open(params.named_obj.prop_file);
      }
      
      QDPIO::cout << "Finished opening prop file" << std::endl;


      // Total number of iterations
      int ncg_had = 0;


#if 1
      // NOTE: not doing this spin rotation, but leaving the code here for reference
      // Rotation from DR to DP
      SpinMatrix diracToDRMat(DiracToDRMat());
      std::vector<MatrixSpinRep_t> diracToDrMatPlus = convertTwoQuarkSpin(diracToDRMat);
      std::vector<MatrixSpinRep_t> diracToDrMatMinus = convertTwoQuarkSpin(adj(diracToDRMat));
#endif


#if 1
      SftMom phases(2, false, decay_dir);
      multi2d<ComplexD> trace_mom(phases.numSubsets(), phases.numMom());
      trace_mom = zero;
#endif

      //
      // Try the factories
      //
      try
      {
	StopWatch swatch;
	swatch.reset();
	QDPIO::cout << "Try the various factories" << std::endl;

	// Typedefs to save typing
	typedef LatticeFermion               T;
	typedef multi1d<LatticeColorMatrix>  P;
	typedef multi1d<LatticeColorMatrix>  Q;

	//
	// Initialize fermion action
	//
	std::istringstream  xml_s(params.param.prop.fermact.xml);
	XMLReader  fermacttop(xml_s);
	QDPIO::cout << "FermAct = " << params.param.prop.fermact.id << std::endl;

	// Generic Wilson-Type stuff
	Handle< FermionAction<T,P,Q> >
	  S_f(TheFermionActionFactory::Instance().createObject(params.param.prop.fermact.id,
							       fermacttop,
							       params.param.prop.fermact.path));

	Handle< FermState<T,P,Q> > state(S_f->createState(u));

	Handle< SystemSolver<LatticeFermion> > PP = S_f->qprop(state,
							       params.param.prop.invParam);
      
	QDPIO::cout << "Suitable factory found: compute all the quark props" << std::endl;
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
	  QDPIO::cout << "t_source = " << t_source << std::endl; 

#if 1
	  if (1)
	  {
	    for(int colorvec_src = 0; colorvec_src < num_vecs; ++colorvec_src)
	    {
	      doTrace(trace_mom,
		      getSrc(source_obj, t_source, colorvec_src),
		      getSoln(soln_obj, t_source, colorvec_src, params.param.contract.mass_label),
		      t_source, colorvec_src,
		      diracToDrMatPlus,
		      diracToDrMatMinus,
		      phases,
		      num_vecs);

	    }
	  }
#endif

	  //
	  // The space distillation loop
	  //
	  for(auto col_src = source_list.begin(); col_src != source_list.end(); ++col_src)
	  {
	    const int colorvec_src = *col_src;
	    
	    StopWatch sniss1;
	    sniss1.reset();
	    sniss1.start();
	    QDPIO::cout << "colorvec_src = " << colorvec_src << std::endl; 

	    // Get the source std::vector
	    LatticeColorVector vec_srce = getSrc(source_obj, t_source, colorvec_src);

	    // Check
	    if (1)
	    {
	      Double ddd = norm2(vec_srce);   // check for debugging
	      QDPIO::cout << __func__ << ": colorvec_src= " << colorvec_src << "  norm(source_vec)= " << ddd << "\n";
	    }

	    //
	    // Loop over each spin source and invert. 
	    // Use the same colorstd::vector source. No spin dilution will be used.
	    //
	    multi2d<LatticeColorVector> ferm_out(Ns,Ns);

	    for(int spin_source=0; spin_source < Ns; ++spin_source)
	    {
	      QDPIO::cout << "spin_source = " << spin_source << std::endl; 

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
			<< " secs" << std::endl;


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


#if 1
	    //
	    // Compute some time-sliced 1pt traces
	    //
	    if (1)
	    {
	      // Get the source std::vector
	      if (colorvec_src < 0)
	      {
		int orig_vec = (-colorvec_src) | (1 << 15);
		vec_srce = getSrc(source_obj, t_source, -orig_vec);
	      }

	      doTrace(trace_mom,
		      vec_srce,
		      ferm_out,
		      t_source, colorvec_src,
		      diracToDrMatPlus,
		      diracToDrMatMinus,
		      phases,
		      num_vecs);
	    }
#endif

	    sniss2.stop();
	    QDPIO::cout << "Time to write propagators for colorvec_src= " << colorvec_src << "  time = " 
			<< sniss2.getTimeInSeconds() 
			<< " secs" << std::endl;
	    
	  } // for colorvec_src
	} // for tt

	swatch.stop();
	QDPIO::cout << "Propagators computed: time= " 
		    << swatch.getTimeInSeconds() 
		    << " secs" << std::endl;
      }
      catch (const std::string& e) 
      {
	QDPIO::cout << name << ": caught exception around qprop: " << e << std::endl;
	QDP_abort(1);
      }

      push(xml_out,"Relaxation_Iterations");
      write(xml_out, "ncg_had", ncg_had);
      pop(xml_out);

      pop(xml_out);  // prop_dist

      snoop.stop();
      QDPIO::cout << name << ": total time = "
		  << snoop.getTimeInSeconds() 
		  << " secs" << std::endl;

      QDPIO::cout << name << ": ran successfully" << std::endl;

      END_CODE();
    } 

  }

} // namespace Chroma
