/*! \file
 * \brief Inline measurement that construct unsmeared hadron nodes using distillation
 */

#include "meas/inline/hadron/inline_unsmeared_hadron_node_distillation_harom_opt2_w.h"

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
#include "meas/smear/link_smearing_aggregate.h"
#include "meas/smear/link_smearing_factory.h"

#include "util/ferm/transf.h"
#include "util/ferm/key_val_db.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"
#include "meas/inline/make_xml_file.h"

#include "util/info/ts_comms.h"
//#include "util/ft/sftmom.h"

#include "meas/inline/io/named_objmap.h"
#include <algorithm>

namespace Chroma 
{ 
  /*!
   * \ingroup hadron
   *
   * @{
   */
  namespace InlineUnsmearedHadronNodeDistillationHaromOpt2Env 
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

      read(inputtop, "fifo", input.fifo );
      read(inputtop, "nodes_per_cn", input.nodes_per_cn );
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

      write(xml, "fifo", input.fifo );
      write(xml, "nodes_per_cn", input.nodes_per_cn );

      pop(xml);
    }


    //! Propagator input
    void read(XMLReader& xml, const std::string& path, Params::Param_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "Propagator", input.prop);
      read(inputtop, "PropSources", input.prop_sources);
      read(inputtop, "Contractions", input.contract);
      read(inputtop, "Displacements", input.displacements);
      read(inputtop, "Moms", input.moms);
      read(inputtop, "SinkSources", input.sink_sources);

      input.link_smearing  = readXMLGroup(inputtop, "LinkSmearing", "LinkSmearingType");

    }

    //! Propagator output
    void write(XMLWriter& xml, const std::string& path, const Params::Param_t& input)
    {
      push(xml, path);

      write(xml, "Propagator", input.prop);
      write(xml, "PropSources", input.prop_sources);
      write(xml, "Contractions", input.contract);
      write(xml, "Displacements", input.displacements);
      write(xml, "Moms", input.moms);
      write(xml, "SinkSourcePairs", input.sink_sources);

      xml << input.link_smearing.xml;

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
  namespace InlineUnsmearedHadronNodeDistillationHaromOpt2Env
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
  namespace InlineUnsmearedHadronNodeDistillationHaromOpt2Env
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

    const std::string name = "UNSMEARED_HADRON_NODE_DISTILLATION_HAROM_OPT2";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= LinkSmearingEnv::registerAll();
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
	      
#if 1
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
      eigen_source.get(src_key, time_slice_io);

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

      XMLReader xml_tmp(xml_in, path);
      std::ostringstream os;
      xml_tmp.print(os);
      xml_str = os.str();
    }


    //----------------------------------------------------------------------------
    //! Generalized propagator operator
    struct KeyGenProp4ElementalOperator_t
    {
      int                t_sink;        /*!< Source time slice */
      int                t_slice;       /*!< Propagator time slice */
      int                t_source;      /*!< Source time slice */
      int                g;             /*!< DR gamma */
      std::vector<int>   displacement;  /*!< Displacement dirs of right colorstd::vector */
      multi1d<int>       mom;           /*!< D-1 momentum of this operator */
      std::string        mass;          /*!< A mass label */
    };

    //! Generalized propagator operator
    struct ValGenProp4ElementalOperator_t
    {
      ValGenProp4ElementalOperator_t() {}
      ValGenProp4ElementalOperator_t( ComplexD* F , int size ): op(F,size,size, Ns,Ns ) {}
      ValGenProp4ElementalOperator_t( int size ): op(size,size, Ns,Ns ) {}
      multi4d<ComplexD>  op;              /*!< Colorstd::vector source and sink with momentum projection */
    };

    //----------------------------------------------------------------------------
    //! Holds key and value as temporaries
    struct KeyValGenProp4ElementalOperator_t
    {
      KeyGenProp4ElementalOperator_t  key;
      ValGenProp4ElementalOperator_t  val;
    };



    //----------------------------------------------------------------------------
    StandardOutputStream& operator<<(StandardOutputStream& os, const KeyGenProp4ElementalOperator_t& d)
    {
      os << "GenProp4:"
	 << " t_sink= " << d.t_sink
	 << " t_slice= " << d.t_slice
	 << " t_source= " << d.t_source
	 << " g= " << d.g
	 << " displacement= " << d.displacement
	 << " mom= " << d.mom
	 << " mass= " << d.mass;

      return os;
    }
    

    //----------------------------------------------------------------------------
    //! KeyGenProp4ElementalOperator reader
    void read(BinaryReader& bin, KeyGenProp4ElementalOperator_t& param)
    {
      read(bin, param.t_sink);
      read(bin, param.t_slice);
      read(bin, param.t_source);
      read(bin, param.g);
      read(bin, param.displacement);
      read(bin, param.mom);
      readDesc(bin, param.mass);
    }

    //! GenProp4ElementalOperator write
    void write(BinaryWriter& bin, const KeyGenProp4ElementalOperator_t& param)
    {
      write(bin, param.t_sink);
      write(bin, param.t_slice);
      write(bin, param.t_source);
      write(bin, param.g);
      write(bin, param.displacement);
      write(bin, param.mom);
      writeDesc(bin, param.mass);
    }


    //----------------------------------------------------------------------------
    //! PropElementalOperator reader
    void read(BinaryReader& bin, ValGenProp4ElementalOperator_t& param)
    {
      read(bin, param.op);
    }

    //! GenProp4ElementalOperator write
    void write(BinaryWriter& bin, const ValGenProp4ElementalOperator_t& param)
    {
      write(bin, param.op);
    }




    //
    //  class Comms -------------------------------------
    //

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



        //----------------------------------------------------------------------------
    //! Normalize just one displacement array
    std::vector<int> normDisp(const std::vector<int>& orig)
    {
      START_CODE();

      std::vector<int> disp;
      std::vector<int> empty; 
      std::vector<int> no_disp(1); no_disp[0] = 0;

      // NOTE: a no-displacement is recorded as a zero-length array
      // Convert a length one array with no displacement into a no-displacement array
      if (orig.size() == 1)
	{
	  if (orig == no_disp)
	    disp = empty;
	  else
	    disp = orig;
	}
      else
	{
	  disp = orig;
	}

      END_CODE();

      return disp;
    } // void normDisp




    void recv_array( int ts , std::vector<int>& array )
    {
      int size = ts_comms_recv( ts );
      array.resize(size);
      for ( int i = 0 ; i < array.size() ; ++i )
	array[i] = ts_comms_recv( ts );
    }

    void recv_array( int ts , multi1d<int>& array )
    {
      int size = ts_comms_recv( ts );
      array.resize(size);
      for ( int i = 0 ; i < array.size() ; ++i )
	array[i] = ts_comms_recv( ts );
    }
    

    void recv_key( int ts , KeyGenProp4ElementalOperator_t& key )
    {
      key.t_sink = ts_comms_recv( ts );
      key.t_source = ts_comms_recv( ts );
      key.g = ts_comms_recv( ts );
      recv_array( ts , key.displacement );
      recv_array( ts , key.mom );
      key.mass = ts_comms_recv_str( ts );
      key.t_slice = ts_comms_recv( ts );
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

      const int ts_per_node = params.param.contract.fifo.size();
      const int nodes_per_cn = params.param.contract.nodes_per_cn;
      const int t_start = params.param.contract.t_start;
      const int Nt_forward = params.param.contract.Nt_forward;

      const int num_vecs  = params.param.contract.num_vecs;
      const int decay_dir = params.param.contract.decay_dir;
      const int Lt        = Layout::lattSize()[decay_dir];
      
      
#if 1
      multi1d<int> ts_lattsize(3);
      ts_lattsize[0] = Layout::lattSize()[0];
      ts_lattsize[1] = Layout::lattSize()[1];
      ts_lattsize[2] = Layout::lattSize()[2];
      
      const int    ts_vol  = Layout::lattSize()[0] * Layout::lattSize()[1] * Layout::lattSize()[2];
      const size_t ts_size_bytes = ts_vol * sizeof( typename TSCollect<LatticeColorVectorSpinMatrix>::TypeSend_t );
      const size_t tensor_size_bytes = sizeof( ComplexD ) * num_vecs * num_vecs * Ns * Ns;
      
      if (Layout::nodeNumber() % nodes_per_cn == 0) {
	int ret = ts_comms_setup( params.param.contract.fifo , std::max( ts_size_bytes , tensor_size_bytes ) );
	if (ret < 0)
	  QDP_error_exit("Error during ts comms setup");

	for ( int i = 0 ; i < ts_per_node ; ++i )
	  {	
	    int base_harom = ts_comms_recv( i ); // read base precision of harom
	    if (base_harom != sizeof(Real) )
	      {
		QDPIO::cout << "harom base precision = " << base_harom << "\n";
		QDPIO::cout << "chroma base precision = " << sizeof(Real) << "\n";
		QDP_error_exit("harom/chroma base precision mismatch");
	      }
	  }
      }
      QMP_barrier();
#endif


      if (Layout::nodeNumber() % nodes_per_cn == 0)
	{
	  for ( int i = 0 ; i < ts_per_node ; ++i )
	    {
	      ts_comms_send( i , params.xml_str );

	      int t_slice = ( t_start + ( Layout::nodeNumber() / nodes_per_cn ) * ts_per_node + i ) % Lt;

	      ts_comms_send( i , t_slice );
	    }
	}
      QMP_barrier();




      
      
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
      BinaryStoreDB< SerialDBKey<KeyGenProp4ElementalOperator_t>, SerialDBData<ValGenProp4ElementalOperator_t> > qdp_db;

      // Open the file, and write the meta-data and the binary for this operator
      if (! qdp_db.fileExists(params.named_obj.dist_op_file))
      {
	XMLBufferWriter file_xml;

	push(file_xml, "DBMetaData");
	write(file_xml, "id", std::string("unsmearedGenprop4ElemOp"));
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
      // Build momentum phases
      //
      QDPIO::cout << "Parse momentum list" << std::endl;
      
      // Possible momenta
      multi2d<int> moms;

      if (params.param.moms.size() > 0)
	{
	  int num_mom = params.param.moms.size();
	  int mom_size = params.param.moms[0].size();
	  moms.resize(num_mom,mom_size);
	  for(int i = 0; i < num_mom; i++)
	    {
	      moms[i] = params.param.moms[i];
	    }
	}
      else
	{
	  QDPIO::cerr << "Error - you have an empty mom list" << std::endl;
	  QDP_abort(1);
	}

#if 0
      // Initialize the slow Fourier transform phases
      QDPIO::cout << "Initialize phases" << std::endl;
      SftMom phases(moms,3);
#endif

      
      

      //
      // Setup the gauge smearing
      //
      QDPIO::cout << "Initalize link smearing" << std::endl;
      multi1d<LatticeColorMatrix> u_smr = u;

      try
      {
	std::istringstream  xml_l(params.param.link_smearing.xml);
	XMLReader  linktop(xml_l);
	QDPIO::cout << "Link smearing type = " << params.param.link_smearing.id << std::endl;
	
	Handle< LinkSmearing >
	  linkSmearing(TheLinkSmearingFactory::Instance().createObject(params.param.link_smearing.id,
								       linktop, 
								       params.param.link_smearing.path));

	QDPIO::cout << "Do the link smearing" << std::endl;
	(*linkSmearing)(u_smr);
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << name << ": Caught Exception link smearing: " << e << std::endl;
	QDP_abort(1);
      }

      QDPIO::cout << "Finished link smearing" << std::endl;

      
#if 0
      for (int i = 0 ; i < Nt_forward ; ++i )
	{
	  int t = ( i + t_start ) % Lt;

	  for (int q = 0 ; q < 3 ; ++q )
	    {
	      LatticeColorMatrix tmp = zero;
	      tmp[ phases.getSet()[t] ] = u_smr[q];
	      QDPIO::cout << "t_slice " << t << "  dir = " << q << " norm2(real(phase * field)) = " << norm2(real( phases[1] * tmp )) << "\n";
	    }
	}
#endif

      //
      // Sending gauge field to harom
      //
      TSCollect<LatticeColorMatrix> ts_collect_u;
      ts_collect_u.exec_prepare( ts_per_node , t_start , Nt_forward , nodes_per_cn );

      std::vector< typename TSCollect<LatticeColorMatrix>::TypeSend_t * > buf_eig( ts_per_node );
      if (Layout::nodeNumber() % nodes_per_cn == 0)
	{
	  for ( int i = 0 ; i < ts_per_node ; ++i )
	    {
	      buf_eig[ i ] = reinterpret_cast< typename TSCollect<LatticeColorMatrix>::TypeSend_t * >( ts_comms_get_shm( i ) );
	    }
	}
      QMP_barrier();
      
      for( int q = 0 ; q < 3 ; ++q )
	{
	  QDPIO::cout << "Sending field, dir = " << q << ", norm = " << norm2( u_smr[q] ) << "\n";
	  
	  ts_collect_u.exec( u_smr[q] , buf_eig );
	      
	  if (Layout::nodeNumber() % nodes_per_cn == 0)
	    {
	      for (int i = 0 ; i < ts_per_node ; ++i )
		{
		  ts_comms_send( i , 24 );
			
		  int ack = ts_comms_recv( i ); // wait for harom to be ready
		  assert( ack == 23 );
		}
	    } // % nodes_per_cn
	  QMP_barrier();
	}


      //return;

      QDPIO::cout << "Check displacement list" << std::endl;
      std::vector< std::vector<int> > disps;

      if (params.param.displacements.size() != 0)
	{
	  for(auto ins = params.param.displacements.begin(); ins != params.param.displacements.end(); ++ins)
	    {	
	      // Make sure displacement is something sensible
	      std::vector<int> disp = normDisp(*ins);
	  
	      disps.push_back(disp);
	    }
	}


      //
      // Try the factories
      //
      try
      {
	StopWatch swatch;
	swatch.reset();


	StopWatch swatch2;
	swatch2.reset();

	// Cache manager
	QDPIO::cout << name << ": initialize the prop factory" << std::endl;
	SourcePropFactory prop_factory(u, params.param.prop, eigen_source, params.param.contract.num_tries);


	TSCollect<LatticeColorVectorSpinMatrix> ts_collect_soln;
	ts_collect_soln.exec_prepare( ts_per_node , t_start , Nt_forward , nodes_per_cn );

	
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

	    snarss1.stop();
	    QDPIO::cout << "Time to generate SOURCE solutions for colorvec_src= " << colorvec_src << "  time = " << snarss1.getTimeInSeconds() << " secs " <<std::endl;


	    QDPIO::cout << "Sending soln to harom: t_source= " << t_source << " colorvec= " << colorvec_src << std::endl;

	    std::vector< typename TSCollect<LatticeColorVectorSpinMatrix>::TypeSend_t * > buf_eig( ts_per_node );
	    if (Layout::nodeNumber() % nodes_per_cn == 0)
	      {
		for ( int i = 0 ; i < ts_per_node ; ++i ) {
		  buf_eig[ i ] = reinterpret_cast< typename TSCollect<LatticeColorVectorSpinMatrix>::TypeSend_t * >( ts_comms_get_shm( i ) );
		}
	      }
	    
	    ts_collect_soln.exec( soln_srce , buf_eig );

	    //
	    // Tell harom the solution vector is ready
	    //
	    if (Layout::nodeNumber() % nodes_per_cn == 0)
	      {
		for (int i = 0 ; i < ts_per_node ; ++i )
		  {
		    ts_comms_send( i , 24 );
		    
		    // optimize, can wait for ack in separate loop
		    int ack = ts_comms_recv( i ); // wait for harom to be ready
		    assert( ack == 23 );
		  }
	      } // % nodes_per_cn
	    QMP_barrier();
	    
	  } // colorvec_src
	  
	} // t_source

	
	swatch.stop();
	QDPIO::cout << "Time to generate all solutions: time = " << swatch.getTimeInSeconds() << " secs " <<std::endl;

	QDPIO::cout << "Receiving tensors from genprop/harom:\n";

	{
	  std::vector<int> do_recv(ts_per_node);

	  if (Layout::nodeNumber() % nodes_per_cn == 0)
	    {
	      for (int i = 0 ; i < ts_per_node ; ++i )
		{	
		  do_recv[i] = ts_comms_recv( i );
		}
	    }
	  QMP_barrier();
	    
	  Integer do_recv_global = do_recv[0];
	  QDPInternal::broadcast( do_recv_global );

	  //printf("node = %d, do_recv_global = %d\n",Layout::nodeNumber(), toInt(do_recv_global) );

	  // Arbitrarily choose the first harom node to determine stopping condition
	  //
	  swatch2.reset(); swatch2.start();

	  while ( toInt(do_recv_global) == 23 )
	    {
	      std::vector< KeyGenProp4ElementalOperator_t >  key( Nt_forward );
	      std::vector< ValGenProp4ElementalOperator_t* >  val( Nt_forward );

	      for (int tcorr = 0 ; tcorr < Nt_forward ; ++tcorr )
		val[ tcorr ] = new ValGenProp4ElementalOperator_t;

	      int key_ref = -1;

	      //
	      // Now wait until all tensors are available in shared memory
	      //				    
	      if (Layout::nodeNumber() % nodes_per_cn == 0)
		{
		  for (int i = 0 ; i < ts_per_node ; ++i )
		    {
		      int tcorr = ( Layout::nodeNumber() / nodes_per_cn ) * ts_per_node + i;

		      recv_key( i , key[ tcorr ] );

		      if ( key_ref == -1 && Layout::nodeNumber() == 0 )
			key_ref = tcorr;

		      //QDPIO::cout << "genprop " << i << ", got key: " << key[i] << "\n";

		      delete val[ tcorr ];
		      val[ tcorr ] = new ValGenProp4ElementalOperator_t( reinterpret_cast< ComplexD * >( ts_comms_get_shm( i ) ) , num_vecs );
		    }
		}
	      QMP_barrier();
	      swatch2.stop();
	      QDPIO::cout << "Time for genprop:                    time = " << swatch2.getTimeInSeconds() << " secs " <<std::endl;


	      swatch2.reset(); swatch2.start();
	      if (Layout::nodeNumber() % nodes_per_cn == 0)
		{
		  for (int tcorr = 0 ; tcorr < Nt_forward ; ++tcorr )
		    {
		      int ts_node = tcorr / ts_per_node * nodes_per_cn;

		      if (ts_node != 0)
			{
			  int size = sizeof( ComplexD ) * num_vecs * num_vecs * Ns * Ns;
			  if (Layout::nodeNumber() == ts_node)
			    {
			      QDPInternal::sendToWait( const_cast<ComplexD*>(val[tcorr]->op.slice(0,0,0)) , 0 , size );
			    }
			  if (Layout::nodeNumber() == 0)
			    {
			      delete val[ tcorr ];
			      val[ tcorr ] = new ValGenProp4ElementalOperator_t (num_vecs );

			      QDPInternal::recvFromWait( const_cast<ComplexD*>(val[tcorr]->op.slice(0,0,0)) , ts_node , size );

			      key[ tcorr ] = key[ key_ref ];
			      key[ tcorr ].t_slice = ( t_start + tcorr ) % Lt;
			    }
			}
		      if ( tcorr == 0 || tcorr == Nt_forward-1 )
			QDPIO::cout << key[ tcorr ] << "\n";
		    }
		}
	      QMP_barrier();
	      swatch2.stop();
	      QDPIO::cout << "Time to copy tensor to primary node: time = " << swatch2.getTimeInSeconds() << " secs " <<std::endl;

	      //
	      // Now store the tensors
	      //				    
	      swatch2.reset(); swatch2.start();
	      for (int tcorr = 0 ; tcorr < Nt_forward ; ++tcorr )
		{
		  qdp_db.insert(key[ tcorr ], *val[ tcorr ]);
		  delete val[ tcorr ];
		}
	      swatch2.stop();
	      QDPIO::cout << "Time to store tensors:               time = " << swatch2.getTimeInSeconds() << " secs " <<std::endl;

		    
	      //
	      // Tell harom we're done with using the shared memory
	      //
	      swatch2.reset(); swatch2.start();
	      if (Layout::nodeNumber() % nodes_per_cn == 0)
		{
		  for (int i = 0 ; i < ts_per_node ; ++i )
		    {
		      ts_comms_send( i , 24 );
		    }
		}
	      QMP_barrier();


	      //
	      // Now make sure the tensors are on the primary node
	      //				    
	      if (Layout::nodeNumber() % nodes_per_cn == 0)
		{
		  for (int i = 0 ; i < ts_per_node ; ++i )
		    {	
		      do_recv[i] = ts_comms_recv( i );
		    }
		}
	      QMP_barrier();
	    
	      do_recv_global = do_recv[0];
	      QDPInternal::broadcast( do_recv_global );
	      
	    } // while

	  swatch2.stop();
	    
	} // node %
	QMP_barrier();
	
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

