/*! \file
 * \brief Inline measurement that construct unsmeared hadron nodes using distillation
 */

#include "meas/inline/hadron/inline_unsmeared_hadron_node_distillation_harom_opt_w.h"

#include "qdp_map_obj_memory.h"
#include "qdp_map_obj_disk.h"
#include "qdp_map_obj_disk_multiple.h"
#include "qdp_disk_map_slice.h"
#include "handle.h"

#include "meas/glue/mesplq.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "util/ferm/key_timeslice_colorvec.h"
#include "util/ferm/key_val_db.h"
#include "util/info/proginfo.h"
#include "io/xml_group_reader.h"
#include "meas/inline/make_xml_file.h"

#include "util/ferm/transf.h"
#include "util/ferm/key_val_db.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"
#include "meas/inline/make_xml_file.h"

#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{ 
  /*!
   * \ingroup hadron
   *
   * @{
   */
  namespace InlineUnsmearedHadronNodeDistillationHaromOptEnv 
  { 
    //! Propagator input
    void read(XMLReader& xml, const std::string& path, Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_id", input.gauge_id);
      read(inputtop, "colorvec_files", input.colorvec_files);
      read(inputtop, "dist_op_file", input.dist_op_file);
    }

    //! Propagator output
    void write(XMLWriter& xml, const std::string& path, const Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "gauge_id", input.gauge_id);
      write(xml, "colorvec_files", input.colorvec_files);
      write(xml, "dist_op_file", input.dist_op_file);

      pop(xml);
    }


    //! Propagator input
    void read(XMLReader& xml, const std::string& path, Params::Param_t::Contract_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "num_vecs", input.num_vecs);
      read(inputtop, "t_start", input.t_start);
      read(inputtop, "Nt_forward", input.Nt_forward);
      read(inputtop, "decay_dir", input.decay_dir);
      read(inputtop, "displacement_length", input.displacement_length);
      read(inputtop, "mass_label", input.mass_label);
      read(inputtop, "num_tries", input.num_tries);
    }

    //! Propagator output
    void write(XMLWriter& xml, const std::string& path, const Params::Param_t::Contract_t& input)
    {
      push(xml, path);

      write(xml, "num_vecs", input.num_vecs);
      write(xml, "t_start", input.t_start);
      write(xml, "Nt_forward", input.Nt_forward);
      write(xml, "decay_dir", input.decay_dir);
      write(xml, "displacement_length", input.displacement_length);
      write(xml, "mass_label", input.mass_label);
      write(xml, "num_tries", input.num_tries);

      pop(xml);
    }


    //! Propagator input
    void read(XMLReader& xml, const std::string& path, Params::Param_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "Propagator", input.prop);
      read(inputtop, "PropSources", input.prop_sources);
      read(inputtop, "Contractions", input.contract);
    }

    //! Propagator output
    void write(XMLWriter& xml, const std::string& path, const Params::Param_t& input)
    {
      push(xml, path);

      write(xml, "Propagator", input.prop);
      write(xml, "PropSources", input.prop_sources);
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
  } // end namespace



  //-------------------------------------------------------------------------------
  //-------------------------------------------------------------------------------
  namespace InlineUnsmearedHadronNodeDistillationHaromOptEnv
  { 
    //----------------------------------------------------------------------------
    // Convenience type
    typedef QDP::MapObjectDisk<KeyTimeSliceColorVec_t, TimeSliceIO<LatticeColorVectorF> > MOD_t;

    // Convenience type
    typedef QDP::MapObjectDiskMultiple<KeyTimeSliceColorVec_t, TimeSliceIO<LatticeColorVectorF> > MODS_t;

    // Convenience type
    typedef QDP::MapObjectMemory<KeyTimeSliceColorVec_t, SubLatticeColorVectorF> SUB_MOD_t;

  } // end namespace


  //-------------------------------------------------------------------------------
  namespace InlineUnsmearedHadronNodeDistillationHaromOptEnv
  { 
    // Anonymous namespace for registration
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

    const std::string name = "UNSMEARED_HADRON_NODE_DISTILLATION_HAROM_OPT";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
    }


    //! Anonymous namespace
    /*! Diagnostic stuff */
    namespace
    {
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

      StandardOutputStream& operator<<(StandardOutputStream& os, const std::vector<int>& d)
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


    //-------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------
    //! Invert off of each source, and do all the checking
    LatticeFermion doInversion(const SystemSolver<LatticeFermion>& PP, const LatticeColorVector& vec_srce, int spin_source, int num_tries)
    {
      //
      // Loop over each spin source and invert. 
      //
      LatticeFermion chi = zero;
      CvToFerm(vec_srce, chi, spin_source);

      LatticeFermion quark_soln = zero;
      int ncg_had;

      // Do the propagator inversion
      // Check if bad things are happening
      bool badP = true;
      for(int nn = 1; nn <= num_tries; ++nn)
      {	
	// Reset
	quark_soln = zero;
	badP = false;
	      
#if 0
	// Solve for the solution std::vector
	SystemSolverResults_t res = PP(quark_soln, chi);
#else
	SystemSolverResults_t res;
        res.n_count = 1;	
        res.resid = 1.0e-7;
#endif
	ncg_had = res.n_count;

	// Some sanity checks
	if (toDouble(res.resid) > 1.0e-3)
	{
	  QDPIO::cerr << name << ": have a resid > 1.0e-3. That is unacceptable" << std::endl;
	  QDP_abort(1);
	}

	// Check for finite values - neither NaN nor Inf
	if (isfinite(quark_soln))
	{
	  // Okay
	  break;
	}
	else
	{
	  QDPIO::cerr << name << ": WARNING - found something not finite, may retry\n";
	  badP = true;
	}
      }

      // Sanity check
      if (badP)
      {
	QDPIO::cerr << name << ": this is bad - did not get a finite solution vector after num_tries= " << num_tries << std::endl;
	QDP_abort(1);
      }

      return quark_soln;
    }
    

    //-------------------------------------------------------------------------------
    class SourcePropFactory
    {
    public:
      SourcePropFactory(const multi1d<LatticeColorMatrix>& u,
			const ChromaProp_t& prop, MODS_t& eigs, int num_tries);

      //! The solution
      void getSoln(LatticeColorVectorSpinMatrix& soln, int t_source, int colorvec_src);

    private:
      //! Eigenvectors
      MODS_t& eigen_source;

      //! Put it here
      int num_tries;

      //! qprop
      Handle< SystemSolver<LatticeFermion> > PP;
    };


    //-------------------------------------------------------------------------------
    SourcePropFactory::SourcePropFactory(const multi1d<LatticeColorMatrix>& u,
					 const ChromaProp_t& prop, MODS_t& eigs_, int num_tries_) : eigen_source(eigs_), num_tries(num_tries_)
    {
      StopWatch swatch;
      swatch.reset();
      swatch.start();
      QDPIO::cout << __func__ << ": initialize the prop cache" << std::endl;

      // Typedefs to save typing
      typedef LatticeFermion               T;
      typedef multi1d<LatticeColorMatrix>  P;
      typedef multi1d<LatticeColorMatrix>  Q;

      //
      // Initialize fermion action
      //
      std::istringstream  xml_s(prop.fermact.xml);
      XMLReader  fermacttop(xml_s);
      QDPIO::cout << "FermAct = " << prop.fermact.id << std::endl;

      // Generic Wilson-Type stuff
      Handle< FermionAction<T,P,Q> >
	S_f(TheFermionActionFactory::Instance().createObject(prop.fermact.id,
							     fermacttop,
							     prop.fermact.path));

      Handle< FermState<T,P,Q> > state(S_f->createState(u));

      PP = S_f->qprop(state, prop.invParam);
      
      swatch.stop(); 
      QDPIO::cout << __func__ << ": finished initializing the prop cache"
		    << "  time = " << swatch.getTimeInSeconds() 
		    << " secs" << std::endl;
    }

  
    //-------------------------------------------------------------------------------
    //! New t-slice
    void SourcePropFactory::getSoln(LatticeColorVectorSpinMatrix& soln, int t_slice, int colorvec_ind)
    {
      // Get the source vector
      StopWatch swatch;
      swatch.reset();
      swatch.start();

      LatticeColorVectorF vec_srce = zero;
      KeyTimeSliceColorVec_t src_key(t_slice, colorvec_ind);
      TimeSliceIO<LatticeColorVectorF> time_slice_io(vec_srce, t_slice);
      if (eigen_source.get(src_key, time_slice_io) != 0)
      {
	QDPIO::cerr << __func__ << ": error retrieving colorvec=" << colorvec_ind << "  t_slice= " << t_slice << std::endl;
	QDP_abort(1);
      }
      
      // Loop over each spin source
      for(int spin_ind=0; spin_ind < Ns; ++spin_ind)
      {	
	QDPIO::cout << __func__ << ": do t_slice= " << t_slice << "  spin_ind= " << spin_ind << "  colorvec_ind= " << colorvec_ind << std::endl; 
	    
	StopWatch snarss1;
	snarss1.reset();
	snarss1.start();

	//
	// Loop over each spin source and invert. 
	// Use the same color vector source. No spin dilution will be used.
	//
	LatticeFermion tmp = doInversion(*PP, vec_srce, spin_ind, num_tries);

	// shove a fermion into a colorvec-spinmatrix 
	FermToProp(tmp, soln, spin_ind);

	snarss1.stop();
	QDPIO::cout << __func__ << ": time to compute prop for t_slice= " << t_slice
		    << "  spin_ind= " << spin_ind
		    << "  colorvec_ind= " << colorvec_ind
		    << "  time = " << snarss1.getTimeInSeconds() 
		    << " secs" << std::endl;
      } // for spin_src

      swatch.stop(); 
      QDPIO::cout << __func__ << ": FINISHED: total time for t_slice = " << t_slice << "  colorvec_ind = " << colorvec_ind
		  << "  time= " << swatch.getTimeInSeconds() << " secs" << std::endl;
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
    //! Generalized propagator operator
    struct KeyGenPropElementalOperator_t
    {
      int                t_slice;       /*!< Propagator time slice */
      int                t_source;      /*!< Source time slice */
      int                t_sink;        /*!< Source time slice */
      int                spin_l;        /*!< Four spin indices */
      int                ins_l;         /*!< Four spin indices */
      int                ins_r;         /*!< Four spin indices */
      int                spin_r;        /*!< Four spin indices */
      std::vector<int>   displacement;  /*!< Displacement dirs of right colorstd::vector */
      multi1d<int>       mom;           /*!< D-1 momentum of this operator */
      std::string        mass;          /*!< A mass label */
    };

    //! Generalized propagator operator
    struct ValGenPropElementalOperator_t
    {
      multi2d<ComplexD>  op;              /*!< Colorstd::vector source and sink with momentum projection */
    };

    //----------------------------------------------------------------------------
    //! Holds key and value as temporaries
    struct KeyValGenPropElementalOperator_t
    {
      KeyGenPropElementalOperator_t  key;
      ValGenPropElementalOperator_t  val;
    };



    //----------------------------------------------------------------------------
    StandardOutputStream& operator<<(StandardOutputStream& os, const KeyGenPropElementalOperator_t& d)
    {
      os << "GenProp:"
	 << " t_sink= " << d.t_sink
	 << " t_slice= " << d.t_slice
	 << " t_source= " << d.t_source
	 << " spin_l= " << d.spin_l
	 << " ins_l= " << d.ins_l
	 << " ins_r= " << d.ins_r
	 << " spin_r= " << d.spin_r
	 << " displacement= " << d.displacement
	 << " mom= " << d.mom
	 << " mass= " << d.mass;

      return os;
    }
    

    //----------------------------------------------------------------------------
    //! KeyGenPropElementalOperator reader
    void read(BinaryReader& bin, KeyGenPropElementalOperator_t& param)
    {
      read(bin, param.t_slice);
      read(bin, param.t_source);
      read(bin, param.t_sink);
      read(bin, param.spin_l);
      read(bin, param.ins_l);
      read(bin, param.ins_r);
      read(bin, param.spin_r);
      read(bin, param.displacement);
      read(bin, param.mom);
      readDesc(bin, param.mass);
    }

    //! GenPropElementalOperator write
    void write(BinaryWriter& bin, const KeyGenPropElementalOperator_t& param)
    {
      write(bin, param.t_slice);
      write(bin, param.t_source);
      write(bin, param.t_sink);
      write(bin, param.spin_l);
      write(bin, param.ins_l);
      write(bin, param.ins_r);
      write(bin, param.spin_r);
      write(bin, param.displacement);
      write(bin, param.mom);
      writeDesc(bin, param.mass);
    }


    //----------------------------------------------------------------------------
    //! PropElementalOperator reader
    void read(BinaryReader& bin, ValGenPropElementalOperator_t& param)
    {
      int n;
      read(bin, n);    // the size is always written, even if 0
      param.op.resize(n,n);
  
      for(int i=0; i < n; ++i)
      {
	for(int j=0; j < n; ++j)
	{
	  read(bin, param.op(i,j));
	}
      }
    }

    //! GenPropElementalOperator write
    void write(BinaryWriter& bin, const ValGenPropElementalOperator_t& param)
    {
      int n = param.op.size1();  // all sizes the same
      write(bin, n);
      for(int i=0; i < n; ++i)
      {
	for(int j=0; j < n; ++j)
	{
	  write(bin, param.op(i,j));
	}
      }
    }

 
    //-------------------------------------------------------------------------------
    // Function call
    void 
    InlineMeas::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
      // If xml file not empty, then use alternate
      if (params.xml_file != "")
      {
	std::string xml_file = makeXMLFileName(params.xml_file, update_no);

	push(xml_out, "HadronNode");
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


    //-------------------------------------------------------------------------------
    // Function call
    void 
    InlineMeas::func(unsigned long update_no,
		     XMLWriter& xml_out) 
    {
      START_CODE();

      StopWatch snoop;
      snoop.reset();
      snoop.start();

      QDPIO::cout << name << ": construct unsmeared hadron nodes via distillation" << std::endl;
      
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

      push(xml_out, "UnsmearedHadronNode");
      write(xml_out, "update_no", update_no);

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
      const int num_vecs  = params.param.contract.num_vecs;
      const int decay_dir = params.param.contract.decay_dir;
      const int Lt        = Layout::lattSize()[decay_dir];

      // A sanity check
      if (decay_dir != Nd-1)
      {
	QDPIO::cerr << name << ": TimeSliceIO only supports decay_dir= " << Nd-1 << "\n";
	QDP_abort(1);
      }

      // Reset
      if (params.param.contract.num_tries <= 0)
      {
	params.param.contract.num_tries = 1;
      }


      //
      // Read in the source along with relevant information.
      // 
      QDPIO::cout << "Snarf the source from a std::map object disk file" << std::endl;

      MODS_t eigen_source;
      eigen_source.setDebug(0);

      std::string eigen_meta_data;   // holds the eigenvalues

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

      //
      // DB storage
      //
      QDPIO::cout << name << ": initialize sdb" << std::endl;
      BinaryStoreDB< SerialDBKey<KeyGenPropElementalOperator_t>, SerialDBData<ValGenPropElementalOperator_t> > qdp_db;

      // Open the file, and write the meta-data and the binary for this operator
      if (! qdp_db.fileExists(params.named_obj.dist_op_file))
      {
	XMLBufferWriter file_xml;

	push(file_xml, "DBMetaData");
	write(file_xml, "id", std::string("unsmearedGenpropElemOp"));
	write(file_xml, "lattSize", QDP::Layout::lattSize());
	write(file_xml, "decay_dir", decay_dir);
	proginfo(file_xml);    // Print out basic program info
	write(file_xml, "Config_info", gauge_xml);
	pop(file_xml);

	std::string file_str(file_xml.str());
	qdp_db.setMaxUserInfoLen(file_str.size());

	qdp_db.open(params.named_obj.dist_op_file, O_RDWR | O_CREAT, 0664);

	qdp_db.insertUserdata(file_str);
      }
      else
      {
	qdp_db.open(params.named_obj.dist_op_file, O_RDWR, 0664);
      }

      QDPIO::cout << "Finished opening distillation file" << std::endl;

      
      // 
      // Build active t_slice list
      //
      std::vector<bool> active_t_slices(Lt);

      if (1)
      {
	// Initialize the active time slices
	for(int t=0; t < Lt; ++t)
	{
	  active_t_slices[t] = false;
	}

	// Forward
	for(int dt=0; dt < params.param.contract.Nt_forward; ++dt)
	{
	  int t = (params.param.contract.t_start + dt) % Lt;
	  active_t_slices[t] = true;
	}
      }



      //
      // Try the factories
      //
      try
      {
	StopWatch swatch;
	swatch.reset();

	// Cache manager
	QDPIO::cout << name << ": initialize the prop factory" << std::endl;
	SourcePropFactory prop_factory(u, params.param.prop, eigen_source, params.param.contract.num_tries);

	//
	// Generate all the solution vectors up front and send them to harom
	//
	for(auto t_source_ptr = params.param.prop_sources.begin(); t_source_ptr != params.param.prop_sources.end(); ++t_source_ptr)
	{
	  int t_source       = *t_source_ptr;

	  QDPIO::cout << "\n\nt_source= " << t_source << std::endl; 
	  swatch.reset(); swatch.start();

	  // Sources
	  for(int colorvec_src=0; colorvec_src < num_vecs; ++colorvec_src)
	  {
	    QDPIO::cout << "Generate SOURCE solutions: colorvec_src = " << colorvec_src << std::endl;
	    
	    StopWatch snarss1;
	    snarss1.reset(); snarss1.start();

	    //
	    // Loop over each spin source and invert. 
	    // Use the same color vector source. No spin dilution will be used.
	    //
	    LatticeColorVectorSpinMatrix soln_srce;
	    prop_factory.getSoln(soln_srce, t_source, colorvec_src);

	    QDPIO::cout << "HACK - send soln to harom: t_source= " << t_source << " colorvec= " << colorvec_src << std::endl;

	    snarss1.stop();
	    QDPIO::cout << "Time to generate SOURCE solutions for colorvec_src= " << colorvec_src << "  time = " << snarss1.getTimeInSeconds() << " secs " <<std::endl;
	  }
	}

	swatch.stop();
	QDPIO::cout << "Time to generate all solutions: time = " << swatch.getTimeInSeconds() << " secs " <<std::endl;


	//
	// Receive genprops from harom
	//
	swatch.reset(); swatch.start();
	
	int t_comm = 0;
	while (t_comm >= 0)
	{
	  KeyGenPropElementalOperator_t  key;
	  ValGenPropElementalOperator_t  val;

	  // Write perambulators
	  //QDPIO::cout << "Write perambulators to disk" << std::endl;
	  qdp_db.insert(key, val);

	  // HACK
	  QDPIO::cout << "HACK - receive genprops from harom\n";
	  t_comm = -1;
	}
	  
	swatch.stop(); 
	QDPIO::cout << "Time to receive and write all genprops: time= " << swatch.getTimeInSeconds() << " secs" <<std::endl;
      }
      catch (const std::string& e) 
      {
	QDPIO::cout << name << ": caught exception around qprop: " << e << std::endl;
	QDP_abort(1);
      }

      // Close db
      qdp_db.close();

      // Close the xml output file
      pop(xml_out);     // UnsmearedHadronNode

      snoop.stop();
      QDPIO::cout << name << ": total time = " 
		  << snoop.getTimeInSeconds() 
		  << " secs" << std::endl;

      QDPIO::cout << name << ": ran successfully" << std::endl;

      END_CODE();
    } // func
  } // namespace InlineUnsmearedHadronNodeDistillationEnv

  /*! @} */  // end of group hadron

} // namespace Chroma

