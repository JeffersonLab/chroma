/*! \file
 * \brief Compute propagators from distillation
 *
 * Propagator calculation in distillation
 */

#include "qdp.h"
#include "fermact.h"
#include "meas/inline/hadron/inline_prop_and_matelem_distillation_harom_w.h"
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

#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <fstream>
#include <iostream>

#include "meas/inline/io/named_objmap.h"

#include "chroma_config.h"

#include "util/info/ts_comms.h"

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
#if 0
void test_mpi_speed()
{
  QDPIO::cout << "beginning test...\n";
  InlinePropAndMatElemDistillationHaromEnv::Comms comms;

  int s;

  int m1=1572864;
  int m2=23592960;
  int m3=14155776;
  int m4=11010048;
  int m5=25165824;

  switch ( Layout::nodeNumber() )
    {
    case 0:
      s = m1;
      comms.add_send_to( 8 , s ); 
      
      s = m2;
      comms.add_receive_from( 13 , s );

      s = m3;
      comms.add_receive_from( 14 , s );
      break;

    case 1:
      break;
    case 2:
      break;
    case 3:
      break;
    case 4:
      break;
    case 5:
      break;
    case 6:
      break;
    case 7:
      break;
    case 8:
      s = m1;
      comms.add_receive_from( 0 , s );

      s = m4;
      comms.add_receive_from( 14 , s );

      s = m5;
      comms.add_receive_from( 15 , s );
      break;

    case 9:
      break;
    case 10:
      break;
    case 11:
      break;
    case 12:
      break;
    case 13:
      s = m2;
      comms.add_send_to( 0 , s );
      break;

    case 14:
      s = m3;
      comms.add_send_to( 0 , s );

      s = m4;
      comms.add_send_to( 8 , s );
      break;

    case 15:
      s = m5;
      comms.add_send_to( 8 , s );
      break;

    }

  comms.finishSetup();

  QDPIO::cout << "comms setup!\n";  

  QMP_barrier();

  StopWatch sniss1;
  sniss1.reset();
  sniss1.start();


  comms.send_receive();
  comms.qmp_wait();

  QMP_barrier();

  sniss1.stop();
  QDPIO::cout << "Time = " << sniss1.getTimeInSeconds() << std::endl;

  int size = m1+m2+m3+m4+m5;
  QDPIO::cout << "MB/s = " << (double)size/1024./1024/(double)sniss1.getTimeInSeconds() << std::endl;

  
}
#endif




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
  }
  

  //----------------------------------------------------------------------------
  namespace InlinePropAndMatElemDistillationHaromEnv 
  {
    //! Propagator input
    void read(XMLReader& xml, const std::string& path, Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_id", input.gauge_id);
      read(inputtop, "colorvec_files", input.colorvec_files);
      read(inputtop, "prop_op_file", input.prop_op_file);
    }

    //! Propagator output
    void write(XMLWriter& xml, const std::string& path, const Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "gauge_id", input.gauge_id);
      write(xml, "colorvec_files", input.colorvec_files);
      write(xml, "prop_op_file", input.prop_op_file);

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
      //read(inputtop, "Nt_backward", input.Nt_backward);
      input.Nt_backward = 0;
      read(inputtop, "mass_label", input.mass_label);
      read(inputtop, "num_tries", input.num_tries);


      input.do_inversions = true;
      if( inputtop.count("do_inversions") == 1 ) {
      	read(inputtop, "do_inversions", input.do_inversions );
      	if (!input.do_inversions)
      	  {
      	    QDPIO::cout << "do_inversions found as false, *** no inverter calls ***\n";
      	  }
      }

      input.check_results = false;
      if( inputtop.count("check_results") == 1 ) {
      	read(inputtop, "check_results", input.check_results );
      	if (input.check_results)
      	  {
      	    QDPIO::cout << "check_results found as true, *** checking results with chroma contractions ***\n";
      	  }
      }

      input.cache_eigs = true;
      if( inputtop.count("cache_eigs") == 1 ) {
      	read(inputtop, "cache_eigs", input.cache_eigs );
      	if (!input.cache_eigs)
      	  {
      	    QDPIO::cout << "cache_eigs found as false, *** no memory caching of eigenvectors. Vectors read twice in total ***\n";
      	  }
      }


      read(inputtop, "fifo", input.fifo );
      read(inputtop, "nodes_per_cn", input.nodes_per_cn );
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
      write(xml, "num_tries", input.num_tries);

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
  } // namespace InlinePropDistillationEnv 


  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  namespace InlinePropAndMatElemDistillationHaromEnv 
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
	SubEigenMap(MODS_t& eigen_source_, int decay_dir ) : eigen_source(eigen_source_), time_slice_set(decay_dir) {}

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
      };

      //----------------------------------------------------------------------------
      //! Getter
      const SubLatticeColorVectorF& SubEigenMap::getVec(int t_source, int colorvec_src) const
      {
	// The key
	KeyTimeSliceColorVec_t src_key(t_source, colorvec_src);
	// If item does not exist, read from original std::map and put in memory std::map
	if (! sub_eigen.exist(src_key))
	{
	  QDPIO::cout << __func__ << ": on t_source= " << t_source << "  colorvec_src= " << colorvec_src << std::endl;

	  // No need to initialize with 'zero' - we are returning a subtype.
	  LatticeColorVectorF vec_srce;

	  TimeSliceIO<LatticeColorVectorF> time_slice_io(vec_srce, t_source);
	  eigen_source.get(src_key, time_slice_io);
	  
	  SubLatticeColorVectorF tmp(getSet()[t_source], vec_srce);

	  sub_eigen.insert(src_key, tmp);
	}
	return sub_eigen[src_key];
      }


      class SubEigenGetter
      {
      public:
	//! Constructor
	SubEigenGetter(MODS_t& eigen_source_, int decay_dir ) : eigen_source(eigen_source_), time_slice_set(decay_dir) {}

	//! The set to be used in sumMulti
	const Set& getSet() const {return time_slice_set.getSet();}

	SubLatticeColorVectorF get( int t_source, int colorvec_src ) const;

      private:
	//! Eigenvectors
	MODS_t& eigen_source;

	// The time-slice set
	TimeSliceSet time_slice_set;
      };

      
      SubLatticeColorVectorF SubEigenGetter::get( int t_source, int colorvec_src ) const
      {
	KeyTimeSliceColorVec_t src_key(t_source, colorvec_src);

	LatticeColorVectorF vec_srce;

	TimeSliceIO<LatticeColorVectorF> time_slice_io(vec_srce, t_source);
	eigen_source.get(src_key, time_slice_io);
	  
	SubLatticeColorVectorF tmp(getSet()[t_source], vec_srce);

	return tmp;
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
  namespace InlinePropAndMatElemDistillationHaromEnv 
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
      
    const std::string name = "PROP_AND_MATELEM_DISTILLATION_HAROM";

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


#if 0    
    struct ts_comms_t
    {
      std::string fifo_send_name;
      std::string fifo_recv_name;

      int         fifo_send_fd;
      int         fifo_recv_fd;
      
      std::string shm_name;
      int         shm_fd;

      void*       ts_buf;
    };
#endif

#if 0    
    namespace
    {
      int local_site(const multi1d<int>& coord, const multi1d<int>& latt_size)
      {
	int order = 0;

	for(int mmu=latt_size.size()-1; mmu >= 1; --mmu)
	  order = latt_size[mmu-1]*(coord[mmu] + order);

	order += coord[0];

	return order;
      }
    }
#endif



    namespace
    {
      std::ostream& operator<< (std::ostream& stream, const KeyPropElementalOperator_t& k) {
	stream << "t_slice = " << k.t_slice 
	       << "\nt_source = " << k.t_source 
	       << "\nspin_src = " << k.spin_src
	       << "\nspin_snk = " << k.spin_snk
	       << "\nmass_label = " << k.mass_label
	       << "\n";
	return stream;
      }

      template<class T>
      ValPropElementalOperator_t& peram_access( T& peram, const KeyPropElementalOperator_t& k )
      {
	return peram[k.t_source][k.t_slice][k.spin_src][k.spin_snk];
      }
    }


    // Real work done here
    void 
    InlineMeas::func(unsigned long update_no,
		     XMLWriter& xml_out) 
    {
      START_CODE();

      //test_mpi_speed();
      // Time = 0.018653
      // MB/s = 3859.96890580604


      StopWatch snoop;
      snoop.reset();
      snoop.start();

      const int ts_per_node = params.param.contract.fifo.size();
      const int nodes_per_cn = params.param.contract.nodes_per_cn;
      const int Nt_forward = params.param.contract.Nt_forward;
      const int Nt = Layout::lattSize()[3];

#if 1
      multi1d<int> ts_lattsize(3);
      ts_lattsize[0] = Layout::lattSize()[0];
      ts_lattsize[1] = Layout::lattSize()[1];
      ts_lattsize[2] = Layout::lattSize()[2];
      
      const int    ts_vol  = Layout::lattSize()[0] * Layout::lattSize()[1] * Layout::lattSize()[2];
      const size_t ts_sizebytes = ts_vol * sizeof( typename TSCollect<LatticeColorVector>::TypeSend_t );
			
      if (Layout::nodeNumber() % nodes_per_cn == 0) {
	int ret = ts_comms_setup( params.param.contract.fifo , ts_sizebytes );
	if (ret < 0)
	  QDP_error_exit("Error during ts comms setup");

      }
#endif


      //test_mpi_speed();


      
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
      SubEigenMap sub_eigen_map(eigen_source, decay_dir);
      QDPIO::cout << "Finished initializing sub-lattice std::map" << std::endl;

      QDPIO::cout << "Initialize sub-eigen getter std::map" << std::endl;
      SubEigenGetter sub_eigen_getter(eigen_source, decay_dir);
      QDPIO::cout << "Finished initializing sub-eigen getter" << std::endl;


      //
      // DB storage
      //
      BinaryStoreDB< SerialDBKey<KeyPropElementalOperator_t>, SerialDBData<ValPropElementalOperator_t> > qdp_db;

      // Open the file, and write the meta-data and the binary for this operator
      if (! qdp_db.fileExists(params.named_obj.prop_op_file))
	{
	  XMLBufferWriter file_xml;
	  push(file_xml, "DBMetaData");
	  write(file_xml, "id", std::string("propElemOp"));
	  write(file_xml, "lattSize", QDP::Layout::lattSize());
	  write(file_xml, "decay_dir", params.param.contract.decay_dir);
	  proginfo(file_xml);    // Print out basic program info
	  write(file_xml, "Params", params.param);
	  write(file_xml, "Config_info", gauge_xml);
	  write(file_xml, "Weights", readEigVals(eigen_meta_data));
	  pop(file_xml);

	  std::string file_str(file_xml.str());
	  qdp_db.setMaxUserInfoLen(file_str.size());

	  qdp_db.open(params.named_obj.prop_op_file, O_RDWR | O_CREAT, 0664);

	  qdp_db.insertUserdata(file_str);
	}
      else
	{
	  qdp_db.open(params.named_obj.prop_op_file, O_RDWR, 0664);
	}

      QDPIO::cout << "Finished opening peram file" << std::endl;

      

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
	  //
	  // First, send the eigenvectors to harom.
	  //
	  if ( Layout::nodeNumber() % nodes_per_cn == 0 )
	    {
	      for ( int i = 0 ; i < ts_per_node ; ++i )
		{
		  ts_comms_send( i , num_vecs );
		}
	    }

	  {
	    TSCollect<LatticeColorVector> ts_eig_collect;
	    ts_eig_collect.exec_prepare( ts_per_node , t_source , Nt_forward , nodes_per_cn );

	    //
	    StopWatch sniss_send;
	    sniss_send.reset();
	    sniss_send.start();

	    for(int colorvec_n=0; colorvec_n < num_vecs; ++colorvec_n)
	      {
		QDPIO::cout << "Sending eig number " << colorvec_n << "\n";
	      
		// sub_eigen_map contains QDPSubTypes, need to collect the pieces first
		//
		StopWatch sniss4;
		sniss4.reset();
		sniss4.start();

		LatticeColorVector eig;
		for ( int t0 = 0 ; t0 < Nt_forward ; ++t0 )
		  {
		    int t = ( t0 + t_source ) % Nt;
		    if (params.param.contract.cache_eigs)
		      eig = sub_eigen_map.getVec( t , colorvec_n );
		    else
		      eig = sub_eigen_getter.get( t , colorvec_n );
		  }

		QMP_barrier();
		sniss4.stop();
		QDPIO::cout << "Time for constructing eigen vector = " << sniss4.getTimeInSeconds() << std::endl;


		std::vector< typename TSCollect<LatticeColorVector>::TypeSend_t * > buf_eig( ts_per_node );
		if (Layout::nodeNumber() % nodes_per_cn == 0)
		  {
		    for ( int i = 0 ; i < ts_per_node ; ++i ) {
		      buf_eig[ i ] = reinterpret_cast< typename TSCollect<LatticeColorVector>::TypeSend_t * >( ts_comms_get_shm( i ) );
		    }
		  }

		ts_eig_collect.exec( eig , buf_eig );
		//auto& recv_eig = ts_eig_collect.getRecv();

		if (Layout::nodeNumber() % nodes_per_cn == 0)
		  {
		    for (int i = 0 ; i < ts_per_node ; ++i )
		      {
			//auto buf_eig = reinterpret_cast< typename TSCollect<LatticeColorVector>::TypeSend_t * >( ts_comms_get_shm( i ) );

			int ts = ( t_source + ( Layout::nodeNumber() / nodes_per_cn ) * ts_per_node + i ) % Nt;

			// for(int site=0; site < ts_vol; ++site)
			//   {
			//     multi1d<int> coord = crtesn(site, ts_lattsize);
			//     int lsite = local_site( coord , ts_lattsize );
			//     buf_eig[ lsite ] = recv_eig[ i ]( coord[0] , coord[1] , coord[2] );
			//   }

			// send the t_slice info, which lets harom know the eig is in shm
			//
			//QDPIO::cout << "sending to slot " << i << "  timeslice info " << ts << "\n";
			ts_comms_send( i , ts );

			// optimize, can wait for ack in separate loop
			int ack = ts_comms_recv( i ); // wait for harom to grab the eig
			assert( ack == 23 );
		      }
		  } // % nodes_per_cn
	      
		QMP_barrier();

	      } // colorvec_n

	    sniss_send.stop();
	    QDPIO::cout << "Time for sending eigen vectors = " << sniss_send.getTimeInSeconds() << std::endl;

	  } // scope to end TSCollect's life
#endif

	  

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

	    std::map< int , std::map< int , std::map< int , std::map< int , ValPropElementalOperator_t > > > > peram;   // [t_source][t_slice][spin_src][spin_snk]
	    std::map< int , std::map< int , std::map< int , std::map< int , ValPropElementalOperator_t > > > > peram2;  // [t_source][t_slice][spin_src][spin_snk]


	    // Initialize
	    for(std::list<KeyPropElementalOperator_t>::const_iterator key = snk_keys.begin();
		key != snk_keys.end();
		++key)
	      {
		peram_access(peram,*key).mat.resize(num_vecs,num_vecs);
		peram_access(peram,*key).mat = zero;

		if (params.param.contract.check_results)
		  {
		    peram_access(peram2,*key).mat.resize(num_vecs,num_vecs);
		    peram_access(peram2,*key).mat = zero;
		  }
	      } // key
	    QDPIO::cout << "peram initialized! " << std::endl; 


	    TSCollect<LatticeColorVector> tscollect;
	    tscollect.exec_prepare( ts_per_node , t_source , Nt_forward , nodes_per_cn );
	    
	    //
	    // The space distillation loop
	    //
	    for(int colorvec_src=0; colorvec_src < num_vecs; ++colorvec_src)
	      {
		StopWatch timer_prep;
		timer_prep.reset();
		timer_prep.start();

		StopWatch sniss1;
		sniss1.reset();
		sniss1.start();

		StopWatch snarss1;
		snarss1.reset();
		snarss1.start();
		QDPIO::cout << "Do spin_source= " << spin_source << "  colorvec_src= " << colorvec_src << std::endl; 

		// Get the source std::vector
		LatticeColorVector vec_srce = zero;
	      
		if (params.param.contract.cache_eigs)
		  vec_srce = sub_eigen_map.getVec( t_source, colorvec_src );
		else
		  vec_srce = sub_eigen_getter.get( t_source, colorvec_src );

		//vec_srce = sub_eigen_map.getVec(t_source, colorvec_src);

		//
		// Loop over each spin source and invert. 
		// Use the same colorstd::vector source. No spin dilution will be used.
		//
		multi1d<LatticeColorVector> ferm_out(Ns);
		LatticeFermion quark_soln = zero;


		timer_prep.stop();
		QDPIO::cout << "Time to prepare source = " << timer_prep.getTimeInSeconds() << "\n";


		if (params.param.contract.do_inversions)
		  {
		    // Insert a ColorVector into spin index spin_source
		    // This only overwrites sections, so need to initialize first
		    LatticeFermion chi = zero;
		    CvToFerm(vec_srce, chi, spin_source);

		    // Do the propagator inversion
		    // Check if bad things are happening
		    bool badP = true;
		    for(int nn = 1; nn <= params.param.contract.num_tries; ++nn)
		      {	
			// Reset
			quark_soln = zero;
			badP = false;
	      
			// Solve for the solution std::vector
			SystemSolverResults_t res = (*PP)(quark_soln, chi);

			ncg_had += res.n_count;

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
			QDPIO::cerr << name << ": this is bad - did not get a finite solution std::vector after num_tries= " 
				    << params.param.contract.num_tries << std::endl;
			QDP_abort(1);
		      }
		  }
		else
		  {
		    //
		    // Not doing inversions, ** timing mode **
		    //
		    for(int spin_sink=0; spin_sink < Ns; ++spin_sink)
		      zero_rep( ferm_out(spin_sink));
		  }
	      
		// Extract into the temporary output array
		for(int spin_sink=0; spin_sink < Ns; ++spin_sink)
		  {
		    ferm_out(spin_sink) = peekSpin(quark_soln, spin_sink);

#if 0
		    for(int site=0; site < Layout::vol(); ++site)
		      {
			multi1d<int> coord = crtesn(site, Layout::lattSize());
			int node	 = Layout::nodeNumber(coord);
			int linear = Layout::linearSiteIndex(coord);
			if (Layout::nodeNumber() == node)
			  {
			    ferm_out(spin_sink).elem(linear).elem().elem(0).real() = coord[0];
			    ferm_out(spin_sink).elem(linear).elem().elem(1).real() = coord[1];
			    ferm_out(spin_sink).elem(linear).elem().elem(2).real() = coord[2];

			    ferm_out(spin_sink).elem(linear).elem().elem(0).imag() = coord[3];
			    ferm_out(spin_sink).elem(linear).elem().elem(1).imag() = 0;
			    ferm_out(spin_sink).elem(linear).elem().elem(2).imag() = 0;
			  }
		      }
#endif

		  }

		snarss1.stop();
		QDPIO::cout << "Time to compute prop for spin_source= " << spin_source << "  colorvec_src= " << colorvec_src << "  time = " 
			    << snarss1.getTimeInSeconds() 
			    << " secs" << std::endl;



		
#if 1
		if (params.param.contract.check_results)
		  {
		    QDPIO::cout << "Doing chroma contractions (as part of checking the harom results)...\n";
		    // the local contraction, for comparison
		    //
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
				for(int site=0; site < Layout::vol(); ++site)
				  {
				    multi1d<int> coord = crtesn(site, Layout::lattSize());
				    int node	 = Layout::nodeNumber(coord);
				    int linear = Layout::linearSiteIndex(coord);
				    if (coord[3] != t_slice) continue;
				  }

				if (params.param.contract.cache_eigs)
				  peram_access(peram2,*key).mat(colorvec_sink,colorvec_src) = innerProduct(sub_eigen_map.getVec(t_slice, colorvec_sink) , ferm_out(key->spin_snk));
				else
				  peram_access(peram2,*key).mat(colorvec_sink,colorvec_src) = innerProduct(sub_eigen_getter.get(t_slice, colorvec_sink) , ferm_out(key->spin_snk));
				

			      } // for colorvec_sink
			  } // for key
		      } // for t_slice
		    QDPIO::cout << "..done\n";
		  }
#endif


		
#if 1
		for(int spin_sink=0; spin_sink < Ns; ++spin_sink)
		  {
		    // 1) gather TS from ferm_out(spin_sink)
		    // 2) TS -> shared mem
		    // 3) t_slice info
		    // 4) recv fifo: done
		      		  
		    // 1) gather TS from ferm_out(spin_sink)

		    QMP_barrier();

		    StopWatch sniss3;
		    sniss3.reset();
		    sniss3.start();

#if 1
		    std::vector< typename TSCollect<LatticeColorVector>::TypeSend_t * > buf_shm( ts_per_node );
		    if (Layout::nodeNumber() % nodes_per_cn == 0)
		      {
			for ( int i = 0 ; i < ts_per_node ; ++i )
			  buf_shm[ i ] = reinterpret_cast< typename TSCollect<LatticeColorVector>::TypeSend_t * >( ts_comms_get_shm( i ) );
		      }

		    tscollect.exec( ferm_out(spin_sink) , buf_shm );


		    //auto& recv = tscollect.getRecv();
#else
		    multi1d< multi3d< typename TSCollect<LatticeColorVector>::TypeSend_t > > recv(ts_per_node);
		    for ( int i = 0 ; i < ts_per_node ; ++i )
		      recv[i].resize( Layout::lattSize()[0],Layout::lattSize()[1],Layout::lattSize()[2] );
#endif


		    QMP_barrier();

		    sniss3.stop();
		    QDPIO::cout << "tscollect time = " << sniss3.getTimeInSeconds() << "\n";


		    if (Layout::nodeNumber() % nodes_per_cn == 0)
		      {
			for (int i = 0 ; i < ts_per_node ; ++i )
			  {
			    int ts = ( t_source + (Layout::nodeNumber() / nodes_per_cn) * ts_per_node + i ) % Nt;
			    ts_comms_send( i , ts );
			  }

			// StopWatch sniss5;
			// sniss5.reset();
			// sniss5.start();

			// Wait for harom instance to finish
			//
			for (int i = 0 ; i < ts_per_node ; ++i )
			  {
			    //QDPIO::cout << "waiting for contraction to be done on slot " << i << "\n";
			    int rcv = ts_comms_recv( i );
			    assert(rcv==23);
			  }

			// sniss5.stop();
			// printf("time waiting for harom to finish = %d secs\n", sniss5.getTimeInSeconds() );


			// StopWatch sniss6;
			// sniss6.reset();
			// sniss6.start();

			for (int i = 0 ; i < ts_per_node ; ++i )
			  {
			    int ts = ( t_source + ( Layout::nodeNumber() / nodes_per_cn ) * ts_per_node + i ) % Nt;

			    KeyPropElementalOperator_t key;
			    key.t_slice    = ts;
			    key.t_source   = t_source;
			    key.spin_src   = spin_source;
			    key.spin_snk   = spin_sink;
			    key.mass_label = params.param.contract.mass_label;

			    auto buf_complex = reinterpret_cast< ComplexD * >( ts_comms_get_shm( i ) );
			    
			    for(int colorvec_sink=0; colorvec_sink < num_vecs; ++colorvec_sink)
			      {
				peram_access(peram,key).mat(colorvec_sink,colorvec_src) = buf_complex[ colorvec_sink ];
			      }
			  }

			// sniss6.stop();
			// printf("time waiting to copy peram = %d secs\n", sniss6.getTimeInSeconds() );


			
		      }  // node % nodes_per_cn

		    QMP_barrier();
		    
		  } //spin_sink   
#endif
	      
		sniss1.stop();
		QDPIO::cout << "Time to compute and assemble peram for spin_source= " << spin_source << "  colorvec_src= " << colorvec_src << "  time = " 
			    << sniss1.getTimeInSeconds()
			    << " secs"
			    << " (for contraction: " << sniss1.getTimeInSeconds() - snarss1.getTimeInSeconds() << ")"
			    << std::endl;

	      } // for colorvec_src
	    

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
#if 1
		QMP_barrier();
		// Make sure the perambulator is available on the primary node.
		// If not, make it available.

		int ts = key->t_slice;

		//QDPIO::cout << "t_slice = " << ts << "\n";

		int tcorr = ( Nt + key->t_slice - t_source ) % Nt;
		int ts_node = tcorr / ts_per_node * nodes_per_cn;


		// if (ts_node == 0)
		//   QDPIO::cout << "primary node\n";
		// else
		//   QDPIO::cout << "on node " << ts_node << ", have to send it\n";

		if (ts_node != 0)
		  {
		    int size = peram_access(peram,*key).mat.size1() * peram_access(peram,*key).mat.size2() * sizeof( ComplexD );

		    // QDPIO::cout << "size1 = " << peram_access(peram,*key).mat.size1() 
		    // 	    << ", size2 = " << peram_access(peram,*key).mat.size2()
		    // 	    << ", bytes = " << size << "\n";


		    if (Layout::nodeNumber() == ts_node)
		      {
			QDPInternal::sendToWait( const_cast<ComplexD*>(peram_access( peram , *key ).mat.slice(0)) , 0 , size );
		      }
		    if (Layout::nodeNumber() == 0)
		      {
			QDPInternal::recvFromWait( const_cast<ComplexD*>(peram_access( peram , *key ).mat.slice(0)) , ts_node , size );
		      }
		  }
		QMP_barrier();
#endif
		    
		// Insert/write to disk
		qdp_db.insert(*key, peram_access(peram,*key));

#if 1
		if (params.param.contract.check_results)
		  {
		    QDPIO::cout << "Comparing chroma and harom peram.. " << std::endl;
		    if (Layout::primaryNode())
		      {
			for (int q=0;q<num_vecs;++q)
			  {
			    for (int w=0;w<num_vecs;++w)
			      {
				if ( toBool(norm2(peram_access(peram,*key).mat(q,w) - peram_access(peram2,*key).mat(q,w)) > 0.0001) )
				  {
				    QDPIO::cout << "mismatch: " << q << "," << w << " \t"
						<< peram_access(peram,*key).mat(q,w) << " vs "
						<< peram_access(peram2,*key).mat(q,w) << "\n";
				  }
			      }
			  }
		      }
		  }
#endif
		    
	      } // for key
		
	    sniss2.stop();
	    QDPIO::cout << "Time to write perambulators for spin_src= " << spin_source << "  time = " 
			<< sniss2.getTimeInSeconds() 
			<< " secs" << std::endl;

	  } // for spin_src

		
	} // for t_source

	//
	// Send harom instances the message to quit
	//
	if (Layout::nodeNumber() % nodes_per_cn == 0)
	  {
	    for (int i = 0 ; i < ts_per_node ; ++i )
	      ts_comms_send( i , -1 );

	    ts_comms_done();
	  }

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







    void* Comms::get_sendbuf( int node )
    {
      assert( map_sendbuf.count(node) == 1 );
      return map_sendbuf.at(node);
    }

    bool Comms::exists_recvbuf( int node )
    {
      return map_recvbuf.count(node) == 1;
    }

    void* Comms::get_recvbuf( int node )
    {
      if ( map_recvbuf.count(node) != 1 )
	{
	  printf("my node %d    node %d count %d \n", Layout::nodeNumber(), node, map_recvbuf.count(node) );
	  // std::cout << "count=" << map_recvbuf.count(node) << "\n";
	  // std::cout << "node=" << node << "\n";
	  // std::cout << "my node=" << Layout::nodeNumber() << "\n";
	}
      assert( map_recvbuf.count(node) == 1 );
      return map_recvbuf.at(node).first;
    }

    int Comms::get_recvbuf_size( int node )
    {
      assert( map_recvbuf.count(node) == 1 );
      return map_recvbuf.at(node).second;
    }

    void Comms::add_receive_from( int node , int size )
    {
      setup_started=true;
  
      void* tmp = static_cast<void*>(new(std::nothrow) char[size]);
      if( tmp == 0x0 ) 
	QDP_error_exit("Unable to allocate recv_buf\n");

      assert( map_recvbuf.count(node) == 0 );
      map_recvbuf[node] = std::make_pair( tmp , size );
  
      buffers.push_back( std::make_pair( tmp , size ) );

      msgmem.push_back( QMP_declare_msgmem( buffers.back().first , size ) );
      if( msgmem.back() == (QMP_msgmem_t)NULL ) { 
	QDP_error_exit("QMP_declare_msgmem for msg[0] failed\n");
      }

      msghandle.push_back( QMP_declare_receive_from( msgmem.back() , node , 0) );
      if( msghandle.back() == (QMP_msghandle_t)NULL ) {
	QDP_error_exit("QMP_declare_receive_from failed\n");
      }
    }

    void Comms::add_send_to( int node , int size )
    {
      setup_started=true;

      void* tmp = static_cast<void*>(new(std::nothrow) char[size]);
      if( tmp == 0x0 )
	QDP_error_exit("Unable to allocate send_buf\n");

      assert( map_sendbuf.count(node) == 0 );
      map_sendbuf[node] = tmp;

      buffers.push_back( std::make_pair( tmp , size ) );

      msgmem.push_back( QMP_declare_msgmem( buffers.back().first , size ) );
      if( msgmem.back() == (QMP_msgmem_t)NULL ) { 
	QDP_error_exit("QMP_declare_msgmem for msg[0] failed\n");
      }

      msghandle.push_back( QMP_declare_send_to( msgmem.back() , node , 0) );
      if( msghandle.back() == (QMP_msghandle_t)NULL ) {
	QDP_error_exit("QMP_declare_receive_from failed\n");
      }
    }

    void Comms::finishSetup()
    {
      setup_finished = true;

      if (!setup_started) {
	//std::cout << "node " << Layout::nodeNumber() << " has nothing to send/receive\n";
	//return;
      }

      mh = QMP_declare_multiple( msghandle.data() , msghandle.size() );
      if( mh == (QMP_msghandle_t)NULL ) { 
	QDP_error_exit("QMP_declare_multiple for mh failed in Map::operator()\n");
      }
    }


    Comms::~Comms()
    {
      cleanup();
    }


    void Comms::cleanup()
    {
      if ( setup_started ) {
	QMP_free_msghandle(mh);
	// QMP_free_msghandle(mh_a[1]);
	// QMP_free_msghandle(mh_a[0]);
	for ( auto& i : msgmem )
	  QMP_free_msgmem(i);
      }

      for ( auto& i : buffers )
	free( i.first );
    }

    void Comms::qmp_wait()
    {
      if (!setup_started) {
	//std::cout << "qmp_wait: node " << Layout::nodeNumber() << " has nothing to send/receive\n";
	//return;
      }

      QMP_status_t err;
      if ((err = QMP_wait(mh)) != QMP_SUCCESS)
	QDP_error_exit(QMP_error_string(err));
    }

    void Comms::fake_comms()
    {
      assert( map_sendbuf.count(0) == 1 );
      assert( map_recvbuf.count(0) == 1 );
  
      float * it = static_cast<float*>(map_sendbuf.at(0));
  
      memcpy( map_recvbuf.at(0).first , map_sendbuf.at(0) , map_recvbuf.at(0).second );
    }

    void Comms::send_receive()
    {
      assert(setup_finished);
      if (!setup_started) {
	//std::cout << "send receive: node " << Layout::nodeNumber() << " has nothing to send/receive\n";
	//return;
      }
  
      QMP_status_t err;
#if QDP_DEBUG >= 3
      QDP_info("Map: send = 0x%x  recv = 0x%x",send_buf,recv_buf);
      QDP_info("Map: establish send=%d recv=%d",destnodes[0],srcenodes[0]);
      {
	const multi1d<int>& me = Layout::nodeCoord();
	multi1d<int> scrd = Layout::getLogicalCoordFrom(destnodes[0]);
	multi1d<int> rcrd = Layout::getLogicalCoordFrom(srcenodes[0]);

	QDP_info("Map: establish-info   my_crds=[%d,%d,%d,%d]",me[0],me[1],me[2],me[3]);
	QDP_info("Map: establish-info send_crds=[%d,%d,%d,%d]",scrd[0],scrd[1],scrd[2],scrd[3]);
	QDP_info("Map: establish-info recv_crds=[%d,%d,%d,%d]",rcrd[0],rcrd[1],rcrd[2],rcrd[3]);
      }
#endif

#if QDP_DEBUG >= 3
      QDP_info("Map: calling start send=%d recv=%d",destnodes[0],srcenodes[0]);
#endif

#ifdef GPU_DEBUG_DEEP
      QDP_info("D2H %d bytes receive buffers",dstnum);
#endif

      if ((err = QMP_start(mh)) != QMP_SUCCESS)
	QDP_error_exit(QMP_error_string(err));
    }


    

  }

} // namespace Chroma

#endif
