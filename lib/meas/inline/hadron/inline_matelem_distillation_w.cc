/*! \file
 * \brief Compute propagators from distillation
 *
 * Propagator calculation in distillation
 */

#include "qdp.h"
#include "fermact.h"
#include "meas/inline/hadron/inline_matelem_distillation_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "qdp_map_obj.h"
#include "qdp_map_obj_disk.h"
#include "qdp_map_obj_disk_multiple.h"
#include "qdp_map_obj_memory.h"
#include "qdp_disk_map_slice.h"
#include "util/ferm/subset_vectors.h"
#include "util/ferm/key_prop_distillation.h"
#include "util/ferm/key_timeslice_colorvec.h"
#include "util/ferm/key_prop_colorvec.h"
#include "util/ferm/key_prop_matelem.h"
#include "util/ferm/key_val_db.h"
#include "util/ferm/transf.h"
#include "util/ferm/spin_rep.h"
#include "util/ferm/diractodr.h"
#include "util/ferm/twoquark_contract_ops.h"
#include "util/ft/sftmom.h"
#include "util/ft/time_slice_set.h"
#include "util/info/proginfo.h"
#include "util/info/proginfo.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"
#include "meas/inline/make_xml_file.h"

#include "meas/inline/io/named_objmap.h"

#include "chroma_config.h"

#if ! defined (QDP_IS_QDPJIT2)

#ifndef QDP_IS_QDPJIT_NO_NVPTX

#ifdef BUILD_JIT_CONTRACTION_KERNELS
#include "custom_kernels/custom_kernels.h"
#endif

namespace Chroma 
{ 
  //----------------------------------------------------------------------------
  // Utilities
  namespace
  {
    
    multi1d<SubsetVectorWeight_t> readEigVals(const std::string& meta)
    {    
      std::istringstream  xml_l(meta);
      XMLReader eigtop(xml_l);

      std::string pat("/MODMetaData/Weights");
      //      multi1d<SubsetVectorWeight_t> eigenvalues;
      multi1d< multi1d<Real> > eigens;
      try
      {
	read(eigtop, pat, eigens);
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << __func__ << ": Caught Exception reading meta= XX" << meta << "XX   with path= " << pat << "   error= " << e << std::endl;
	QDP_abort(1);
      }

      eigtop.close();

      multi1d<SubsetVectorWeight_t> eigenvalues(eigens.size());

      for(int i=0; i < eigens.size(); ++i)
	eigenvalues[i].weights = eigens[i];

      return eigenvalues;
    }
  } // end anonymous namespace



  //----------------------------------------------------------------------------
  namespace InlineMatElemDistillationEnv 
  {
    //! Propagator input
    void read(XMLReader& xml, const std::string& path, Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "colorvec_files", input.colorvec_files);
      read(inputtop, "prop_file", input.prop_file);
      read(inputtop, "peram_file", input.peram_file);
    }

    //! Propagator output
    void write(XMLWriter& xml, const std::string& path, const Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "colorvec_files", input.colorvec_files);
      write(xml, "prop_file", input.prop_file);
      write(xml, "peram_file", input.peram_file);

      pop(xml);
    }


    //! Propagator input
    void read(XMLReader& xml, const std::string& path, Params::Param_t::Contract_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "num_vecs", input.num_vecs);
      read(inputtop, "lt_orig", input.lt_orig);
      read(inputtop, "t_offset", input.t_offset);
      read(inputtop, "t_sources", input.t_sources);
      read(inputtop, "decay_dir", input.decay_dir);
      read(inputtop, "Nt_forward", input.Nt_forward);
      read(inputtop, "Nt_backward", input.Nt_backward);
      read(inputtop, "mass_label", input.mass_label);

      input.zero_colorvecs = false;
      if( inputtop.count("zero_colorvecs") == 1 ) {
	read(inputtop, "zero_colorvecs", input.zero_colorvecs );
	if (input.zero_colorvecs)
	  {
	    QDPIO::cout << "zero_colorvecs found, *** timing mode activated ***\n";
	  }
      }

      input.fuse_timeloop = false;
      if( inputtop.count("fuse_timeloop") == 1 ) {
	read(inputtop, "fuse_timeloop", input.fuse_timeloop );
	if (input.fuse_timeloop)
	  {
	    QDPIO::cout << "fuse_timeloop found!\n";
	  }
      }

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

      pop(xml);
    }


    //! Propagator input
    void read(XMLReader& xml, const std::string& path, Params::Param_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "Contractions", input.contract);
    }

    //! Propagator output
    void write(XMLWriter& xml, const std::string& path, const Params::Param_t& input)
    {
      push(xml, path);

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
  } // namespace InlinePropDistillationEnv 


  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  namespace InlineMatElemDistillationEnv 
  {
    //----------------------------------------------------------------------------
    // Convenience type
    typedef QDP::MapObjectDisk<KeyTimeSliceColorVec_t, TimeSliceIO<LatticeColorVectorF> > MOD_t;

    // Convenience type
    typedef QDP::MapObjectDiskMultiple<KeyTimeSliceColorVec_t, TimeSliceIO<LatticeColorVectorF> > MODS_t;

    // Convenience type
    typedef QDP::MapObjectMemory<KeyTimeSliceColorVec_t, SubLatticeColorVectorF> SUB_MOD_t;


    
    // Anonymous namespace
    namespace
    {
      //----------------------------------------------------------------------------
      //! Class to hold std::map of eigenvectors
      class SubEigenMap
      {
      public:
	//! Constructor
	SubEigenMap(MODS_t& eigen_source_, int decay_dir, bool zero_colorvecs, int t_offset, int lt_orig) : eigen_source(eigen_source_), time_slice_set(decay_dir), zero_colorvecs(zero_colorvecs), t_offset(t_offset), lt_orig(lt_orig) {}

	//! Getter
	const SubLatticeColorVectorF& getVec(int t_source, int colorvec_src) const;

	//! The set to be used in sumMulti
	const Set& getSet() const {return time_slice_set.getSet();}

      private:
	//! Eigenvectors
	MODS_t& eigen_source;

	// The time-slice set
	TimeSliceSet time_slice_set;

      private:
	//! Where we store the sublattice versions
	mutable SUB_MOD_t sub_eigen;
	bool zero_colorvecs;
	int t_offset;
	int lt_orig;
      };

      //----------------------------------------------------------------------------
      //! Getter
      const SubLatticeColorVectorF& SubEigenMap::getVec(int t_source, int colorvec_src) const
      {
	// The key
	KeyTimeSliceColorVec_t src_key( (t_source + t_offset ) % lt_orig , colorvec_src);

	// If item does not exist, read from original std::map and put in memory std::map
	if (! sub_eigen.exist(src_key))
	{
	  QDPIO::cout << __func__ << ": on t_source= " << t_source << "  colorvec_src= " << colorvec_src << std::endl;

	  // No need to initialize with 'zero' - we are returning a subtype.
	  LatticeColorVectorF vec_srce;

	  if (!zero_colorvecs)
	    {
	      TimeSliceIO<LatticeColorVectorF> time_slice_io(vec_srce, t_source );
	      eigen_source.get(src_key, time_slice_io);
	    }
	  
	  SubLatticeColorVectorF tmp(getSet()[t_source], vec_srce);

	  sub_eigen.insert(src_key, tmp);
	}
	return sub_eigen[src_key];
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
      //! Get sink keys
      std::list<KeyPropElementalOperator_t> getSnkKeys(int t_source, int spin_source, int Nt_forward, int Nt_backward, const std::string mass)
      {
	std::list<KeyPropElementalOperator_t> keys;

	std::vector<bool> active_t_slices = getActiveTSlices(t_source, Nt_forward, Nt_backward);
	
	const int decay_dir = Nd-1;
	const int Lt = Layout::lattSize()[decay_dir];
	
	for(int spin_sink=0; spin_sink < Ns; ++spin_sink)
	{
	  for(int t=0; t < Lt; ++t)
	  {
	    if (! active_t_slices[t]) {continue;}

	    KeyPropElementalOperator_t key;

	    key.t_source     = t_source;
	    key.t_slice      = t;
	    key.spin_src     = spin_source;
	    key.spin_snk     = spin_sink;
	    key.mass_label   = mass;

	    keys.push_back(key);
	  } // for t
	} // for spin_sink

	return keys;
      }
	
    } // end anonymous
  } // end namespace

	
  //----------------------------------------------------------------------------
  namespace InlineMatElemDistillationEnv 
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
      
    const std::string name = "MATELEM_DISTILLATION";

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

      push(xml_out, "PropDistillation");
      write(xml_out, "update_no", update_no);

      proginfo(xml_out);    // Print out basic program info

      // Write out the input
      write(xml_out, "Input", params);

      push(xml_out, "Output_version");
      write(xml_out, "out_version", 1);
      pop(xml_out);

      // Will use TimeSliceSet-s a lot
      const int decay_dir = params.param.contract.decay_dir;
      const int Lt        = Layout::lattSize()[decay_dir];
      const int lt_orig   = params.param.contract.lt_orig;
      const int t_offset  = params.param.contract.t_offset;

      // A sanity check
      if (decay_dir != Nd-1)
      {
	QDPIO::cerr << name << ": TimeSliceIO only supports decay_dir= " << Nd-1 << "\n";
	QDP_abort(1);
      }


      QDP::MapObjectDisk<KeyPropDistillation_t, TimeSliceIO<LatticeColorVectorF> > prop_obj;
      prop_obj.setDebug(0);

      QDPIO::cout << "Open source prop file" << std::endl;

      if (! prop_obj.fileExists(params.named_obj.prop_file))
      {
	QDPIO::cerr << name << ": source file does not exist: src_file= " << params.named_obj.prop_file << std::endl;
	QDP_abort(1);
      }
      else
      {
	prop_obj.open(params.named_obj.prop_file, std::ios_base::in);
      }

      QDPIO::cout << "Finished opening solution file" << std::endl;

      
      //
      // Read in the source along with relevant information.
      // 
      QDPIO::cout << "Snarf the source from a std::map object disk file" << std::endl;

      MODS_t eigen_source;
      eigen_source.setDebug(0);

      std::string eigen_meta_data;   // holds the eigenvalues

      if (!params.param.contract.zero_colorvecs)
	{
	  try
	    {
	      // Open
	      QDPIO::cout << "Open file= " << params.named_obj.colorvec_files[0] << std::endl;
	      eigen_source.open(params.named_obj.colorvec_files);

	      // Snarf the source info. 
	      QDPIO::cout << "Get user data" << std::endl;
	      eigen_source.getUserdata(eigen_meta_data);
	      //	QDPIO::cout << "User data= " << eigen_meta_data << std::endl;

	      // Write it
	      //	QDPIO::cout << "Write to an xml file" << std::endl;
	      //	XMLBufferWriter xml_buf(eigen_meta_data);
	      //	write(xml_out, "Source_info", xml_buf);
	    }    
	  catch (std::bad_cast) {
	    QDPIO::cerr << name << ": caught dynamic cast error" << std::endl;
	    QDP_abort(1);
	  }
	  catch (const std::string& e) {
	    QDPIO::cerr << name << ": error extracting source_header: " << e << std::endl;
	    QDP_abort(1);
	  }
	  catch( const char* e) {
	    QDPIO::cerr << name <<": Caught some char* exception:" << std::endl;
	    QDPIO::cerr << e << std::endl;
	    QDPIO::cerr << "Rethrowing" << std::endl;
	    throw;
	  }

	  QDPIO::cout << "Source successfully read and parsed" << std::endl;
	}
      
#if 0
      // Sanity check
      if (params.param.contract.num_vecs > eigen_source.size())
      {
	QDPIO::cerr << name << ": number of available eigenvectors is too small\n";
	QDP_abort(1);
      }
#endif

      QDPIO::cout << "Number of vecs available is large enough" << std::endl;

      // The sub-lattice eigenstd::vector std::map
      QDPIO::cout << "Initialize sub-lattice std::map" << std::endl;
      SubEigenMap sub_eigen_map(eigen_source, decay_dir, params.param.contract.zero_colorvecs, t_offset , lt_orig );
      QDPIO::cout << "Finished initializing sub-lattice std::map" << std::endl;


      //
      // DB storage
      //
      BinaryStoreDB< SerialDBKey<KeyPropElementalOperator_t>, SerialDBData<ValPropElementalOperator_t> > qdp_db;

      if (!params.param.contract.zero_colorvecs)
	{
	  // Open the file, and write the meta-data and the binary for this operator
	  if (! qdp_db.fileExists(params.named_obj.peram_file))
	    {
	      XMLBufferWriter file_xml;
	      push(file_xml, "DBMetaData");
	      write(file_xml, "id", std::string("propElemOp"));

	      auto tmp = QDP::Layout::lattSize();
	      tmp[Nd-1] = lt_orig;
	      
	      write(file_xml, "lattSize", tmp );
	      write(file_xml, "decay_dir", params.param.contract.decay_dir);
	      proginfo(file_xml);    // Print out basic program info
	      write(file_xml, "Params", params.param);
	      if (!params.param.contract.zero_colorvecs)
		write(file_xml, "Weights", readEigVals(eigen_meta_data));
	      pop(file_xml);

	      std::string file_str(file_xml.str());
	      qdp_db.setMaxUserInfoLen(file_str.size());

	      qdp_db.open(params.named_obj.peram_file, O_RDWR | O_CREAT, 0664);

	      qdp_db.insertUserdata(file_str);
	    }
	  else
	    {
	      qdp_db.open(params.named_obj.peram_file, O_RDWR, 0664);
	    }

	  QDPIO::cout << "Finished opening peram file" << std::endl;
	}

      // Total number of iterations
      int ncg_had = 0;

      //
      // Try the factories
      //
      try
      {
	StopWatch swatch;
	swatch.reset();
	QDPIO::cout << "Try the various factories" << std::endl;

	swatch.start();

#ifdef BUILD_JIT_CONTRACTION_KERNELS
	if ( params.param.contract.fuse_timeloop )
	  QDPIO::cout << "Using JIT contraction kernels (fused timeloop)\n";
	else
	  QDPIO::cout << "Using JIT contraction kernels (non-fused timeloop)\n";
#endif

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

	  t_source -= t_offset;
	  QDPIO::cout << "new t_source in current lattice = " << t_source << std::endl;
	  
	  // Loop over each spin source
	  for(int spin_source=0; spin_source < Ns; ++spin_source)
	  {
	    QDPIO::cout << "spin_source = " << spin_source << std::endl; 

	    // These are the common parts of a perambulator that are needed for this time source
	    std::list<KeyPropElementalOperator_t> snk_keys(getSnkKeys(t_source,
								      spin_source,
								      params.param.contract.Nt_forward,
								      params.param.contract.Nt_backward,
								      params.param.contract.mass_label));

	    if (1) {
	      // The final perambulator
	      QDP::MapObjectMemory<KeyPropElementalOperator_t, ValPropElementalOperator_t> peram;

	      // Initialize
	      for(std::list<KeyPropElementalOperator_t>::const_iterator key = snk_keys.begin();
		  key != snk_keys.end();
		  ++key)
		{
		  // The perambulator value
		  ValPropElementalOperator_t tmp;
		  peram.insert(*key, tmp);
	      
		  peram[*key].mat.resize(num_vecs,num_vecs);
		  peram[*key].mat = zero;
		} // key
	      QDPIO::cout << "peram initialized! " << std::endl; 

	      //
	      // The space distillation loop
	      //
	      for(int colorvec_src=0; colorvec_src < num_vecs; ++colorvec_src)
		{
		  StopWatch sniss1;
		  sniss1.reset();
		  sniss1.start();

		  StopWatch snarss1;
		  snarss1.reset();
		  snarss1.start();
		  QDPIO::cout << "Do spin_source= " << spin_source << "  colorvec_src= " << colorvec_src << std::endl; 


		  //
		  // Loop over each spin source and invert. 
		  // Use the same colorstd::vector source. No spin dilution will be used.
		  //
		  multi1d<LatticeColorVector> ferm_out(Ns);

		  if (!params.param.contract.zero_colorvecs)
		    {
		      // Extract into the temporary output array
		      for(int spin_sink=0; spin_sink < Ns; ++spin_sink)
			{
			  for ( int t = 0 ; t < Lt ; ++t )
			    {
			      KeyPropDistillation_t key;

			      key.t_source = ( t_source + t_offset ) % lt_orig;
			      key.t_slice = ( t + t_offset ) % lt_orig;
			      key.colorvec_src = colorvec_src;
			      key.spin_src = spin_source;
			      key.spin_snk = spin_sink;
			      key.mass = params.param.contract.mass_label;

			      LatticeColorVectorF tmp;

			      TimeSliceIO<LatticeColorVectorF> time_slice_io(tmp, t);

			      prop_obj.get( key , time_slice_io );
			      
			      ferm_out(spin_sink)[ sub_eigen_map.getSet()[t] ] = tmp;
			    }
			}

		    } // zero_colorvecs ??
		  else
		    {
		      for ( int s = 0 ; s < Ns ; ++s )
			{
			  zero_rep( ferm_out[s] );
			}
		    }
	      

		  snarss1.stop();
		  QDPIO::cout << "Time to read prop for spin_source= " << spin_source << "  colorvec_src= " << colorvec_src << "  time = " 
			      << snarss1.getTimeInSeconds()
			      << " secs" << std::endl;

		  // The perambulator part
		  // Loop over time

#ifndef BUILD_JIT_CONTRACTION_KERNELS
		  for(int t_slice = 0; t_slice < Lt; ++t_slice)
		    {
		      // Loop over all the keys
		      for(std::list<KeyPropElementalOperator_t>::const_iterator key= snk_keys.begin();
			  key != snk_keys.end();
			  ++key)
			{
			  if (key->t_slice != t_slice) {continue;}

			  // Loop over the sink colorvec, form the innerproduct and the resulting perambulator
			  for(int colorvec_sink=0; colorvec_sink < num_vecs; ++colorvec_sink)
			    {
			      peram[*key].mat(colorvec_sink,colorvec_src) = innerProduct(sub_eigen_map.getVec(t_slice, colorvec_sink), 
											 ferm_out(key->spin_snk));
			    } // for colorvec_sink
			} // for key
		    } // for t_slice
#else
		  if ( params.param.contract.fuse_timeloop)
		    {
		      //
		      // fused timeloop
		      //
		      for(int spin_snk = 0; spin_snk < Ns; ++spin_snk)
			{
			  int count = 0;
			  for(std::list<KeyPropElementalOperator_t>::const_iterator key= snk_keys.begin();  key != snk_keys.end(); ++key)
			    {
			      if (key->spin_snk != spin_snk) { continue;}
			      count += num_vecs;
			    }

			  //QDPIO::cout << "spin_snk = " << spin_snk << ", count = " << count << "\n";
		      
			  multi1d<SubLatticeColorVectorF*> vec_ptr( count );
			  multi1d<ComplexD*> contr_ptr( count );

			  int run_count = 0;
		      
			  // Loop over all the keys
			  for(std::list<KeyPropElementalOperator_t>::const_iterator key= snk_keys.begin();
			      key != snk_keys.end();
			      ++key)
			    {
			      if (key->spin_snk != spin_snk) {continue;}

			      //
			      // Pack pointers to the vectors and matrix elements
			      //
			      for (int i=0 ; i < num_vecs ; ++i ) {
				vec_ptr[run_count] = const_cast<SubLatticeColorVectorF*>( &sub_eigen_map.getVec( key->t_slice , i ) );
				contr_ptr[run_count] = &peram[*key].mat( i , colorvec_src );
				++run_count;
			      }
			    }
		      
			  //
			  // Big call-out 
			  //
			  multi_innerProduct( contr_ptr , vec_ptr , ferm_out(spin_snk) );
		  
			} // for spin_snk
		    }
		  else
		    {
		      //
		      // non-fused timeloop
		      //
		      for(int t_slice = 0; t_slice < Lt; ++t_slice)
			{
			  // Loop over all the keys
			  for(std::list<KeyPropElementalOperator_t>::const_iterator key= snk_keys.begin();
			      key != snk_keys.end();
			      ++key)
			    {
			      if (key->t_slice != t_slice) {continue;}
			      //
			      // Pack pointers to the vectors and matrix elements
			      //
			      multi1d<SubLatticeColorVectorF*> vec_ptr( num_vecs );
			      multi1d<ComplexD*> contr_ptr( num_vecs );
			      for (int i=0 ; i < num_vecs ; ++i ) {
				vec_ptr[i] = const_cast<SubLatticeColorVectorF*>( &sub_eigen_map.getVec( t_slice , i ) );
				contr_ptr[i] = &peram[*key].mat( i , colorvec_src );
			      }
		  
			      //
			      // Big call-out 
			      //
			      multi_innerProduct( contr_ptr , vec_ptr , ferm_out(key->spin_snk) );
			    } // for key
			} // for t_slice
		    }
#endif

		  sniss1.stop();
		  QDPIO::cout << "Time to compute and assemble peram for spin_source= " << spin_source << "  colorvec_src= " << colorvec_src << "  time = " 
			      << sniss1.getTimeInSeconds()
			      << " secs"
			      << " (for contraction: " << sniss1.getTimeInSeconds() - snarss1.getTimeInSeconds() << ")"
			      << std::endl;

		} // for colorvec_src

	      if (!params.param.contract.zero_colorvecs)
		{
		  // Write out each time-slice chunk of a lattice colorvec soln to disk
		  QDPIO::cout << "Write perambulator for spin_source= " << spin_source << "  to disk" << std::endl;
		  StopWatch sniss2;
		  sniss2.reset();
		  sniss2.start();

		  // The perambulator is complete. Write it.
		  for(std::list<KeyPropElementalOperator_t>::const_iterator key= snk_keys.begin();
		      key != snk_keys.end();
		      ++key)
		    {
		      KeyPropElementalOperator_t tmp = *key;
		      tmp.t_source = ( tmp.t_source + t_offset ) % lt_orig;
		      tmp.t_slice  = ( tmp.t_slice + t_offset ) % lt_orig;
		      
		      // Insert/write to disk
		      qdp_db.insert( tmp , peram[*key]);

		    } // for key
		
		  sniss2.stop();
		  QDPIO::cout << "Time to write perambulators for spin_src= " << spin_source << "  time = " 
			      << sniss2.getTimeInSeconds() 
			      << " secs" << std::endl;
		}

	    } // this ends the lifetime of qdp_db.
	  } // for spin_src
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

#endif
#endif
