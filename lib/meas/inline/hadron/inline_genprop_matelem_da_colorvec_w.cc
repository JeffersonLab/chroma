//
/*! \file
 * \brief Compute the matrix element of  LatticeColorVector*M^-1*Gamma*displace*M^-1**LatticeColorVector
 *
 * Generalized propagator calculation on a colorvector for meson distribution amplitudes and parton dstribution functions
 */

#include "handle.h"
#include "inline_genprop_matelem_da_colorvec_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/smear/link_smearing_aggregate.h"
#include "meas/smear/link_smearing_factory.h"
#include "meas/smear/displace.h"
#include "meas/glue/mesplq.h"
#include "qdp_map_obj.h"
#include "util/ferm/key_val_db.h"
#include "util/ferm/key_prop_colorvec.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "meas/inline/make_xml_file.h"

#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{ 
  /*!
   * \ingroup hadron
   *
   * @{
   */
  namespace InlineGenPropMatElemDAColorVecEnv  
  { 
    // Reader for input parameters
    // Reader for input parameters
    void read(XMLReader& xml, const string& path, Params::Param_t& param)
    {
      XMLReader paramtop(xml, path);
    
      int version;
      read(paramtop, "version", version);

      switch (version) 
      {
      case 1:
	read(paramtop, "t_source", param.t_source);
	read(paramtop, "t_sink", param.t_sink);
	read(paramtop, "mom_max", param.mom_max);
	read(paramtop, "boost_dir", param.boost_dir);
	read(paramtop, "gamma", param.gamma);
	read(paramtop, "num_vecs", param.num_vecs);
	read(paramtop, "decay_dir", param.decay_dir);
	read(paramtop, "mass_label", param.mass_label);
	read(paramtop, "restrict_plateau", param.restrict_plateau);
	break ;
      default :
	QDPIO::cerr << "Input parameter version " << version << " unsupported." << endl;
	QDP_abort(1);
      }


    }


    // Writer for input parameters
    void write(XMLWriter& xml, const string& path, const Params::Param_t& param)
    {
      push(xml, path);

      int version = 1;

      write(xml, "version", version);
      write(xml, "t_source", param.t_source);
      write(xml, "t_sink", param.t_sink);
      write(xml, "mom_max", param.mom_max);
      write(xml, "boost_dir", param.boost_dir);
      write(xml, "gamma", param.gamma);
      write(xml, "num_vecs", param.num_vecs);
      write(xml, "decay_dir", param.decay_dir);
      write(xml, "mass_label", param.mass_label);
      write(xml, "restrict_plateau", param.restrict_plateau);

      pop(xml);
    }

    //! Read named objects 
    void read(XMLReader& xml, const string& path, Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_id", input.gauge_id);
      read(inputtop, "sm_gauge_id", input.sm_gauge_id);
      read(inputtop, "source_prop_id", input.source_prop_id);
      read(inputtop, "sink_prop_id", input.sink_prop_id);
      read(inputtop, "genprop_op_file", input.genprop_op_file);
    }

    //! Write named objects
    void write(XMLWriter& xml, const string& path, const Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "gauge_id", input.gauge_id);
      write(xml, "sm_gauge_id", input.sm_gauge_id);
      write(xml, "source_prop_id", input.source_prop_id);
      write(xml, "sink_prop_id", input.sink_prop_id);
      write(xml, "genprop_op_file", input.genprop_op_file);

      pop(xml);
    }

    // Writer for input parameters
    void write(XMLWriter& xml, const string& path, const Params& param)
    {
      param.writeXML(xml, path);
    }
  }


  namespace InlineGenPropMatElemDAColorVecEnv  
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

    const std::string name = "GENPROP_MATELEM_DA_COLORVEC";

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
    }


    //----------------------------------------------------------------------------
    // Param stuff
    Params::Params()
    { 
      frequency = 0; 
      param.mom_max = 0;
    }

    Params::Params(XMLReader& xml_in, const std::string& path) 
    {
      try 
      {
	XMLReader paramtop(xml_in, path);

	if (paramtop.count("Frequency") == 1)
	  read(paramtop, "Frequency", frequency);
	else
	  frequency = 1;

	// Read program parameters
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


    void
    Params::writeXML(XMLWriter& xml_out, const std::string& path) const
    {
      push(xml_out, path);
    
      // Parameters for source construction
      write(xml_out, "Param", param);

      // Write out the output propagator/source configuration info
      write(xml_out, "NamedObject", named_obj);

      pop(xml_out);
    }


    //----------------------------------------------------------------------------
    //! Generalized propagator operator
    struct KeyGenPropElementalOperator_t
    {
      int                t_slice;       /*!< Propagator time slice */
      int                t_source;      /*!< Source time slice */
      int                t_sink;        /*!< Source time slice */
      int                spin_l;        /*!< Source spin index */
      int                spin_r;        /*!< Sink spin index */
      int                gamma;         /*!< The gamma matrix number - [0,Ns^2) */
      multi1d<int>       displacement;  /*!< Displacement dirs of right colorvector */
      multi1d<int>       mom;           /*!< D-1 momentum of this operator */
      std::string        mass_label;    /*!< A mass label */
    };

    //! Generalized propagator operator
    struct ValGenPropElementalOperator_t
    {
      multi2d<ComplexD>  op;              /*!< Colorvector source and sink with momentum projection */
    };


    //----------------------------------------------------------------------------
    //! Holds key and value as temporaries
    struct KeyValGenPropElementalOperator_t
    {
      SerialDBKey<KeyGenPropElementalOperator_t>  key;
      SerialDBData<ValGenPropElementalOperator_t> val;
    };


    //----------------------------------------------------------------------------
    //! KeyGenPropElementalOperator reader
    void read(BinaryReader& bin, KeyGenPropElementalOperator_t& param)
    {
      read(bin, param.t_slice);
      read(bin, param.t_source);
      read(bin, param.t_sink);
      read(bin, param.spin_l);
      read(bin, param.spin_r);
      read(bin, param.gamma);
      read(bin, param.displacement);
      read(bin, param.mom);
      read(bin, param.mass_label, 32);
    }

    //! GenPropElementalOperator write
    void write(BinaryWriter& bin, const KeyGenPropElementalOperator_t& param)
    {
      write(bin, param.t_slice);
      write(bin, param.t_source);
      write(bin, param.t_sink);
      write(bin, param.spin_l);
      write(bin, param.spin_r);
      write(bin, param.gamma);
      write(bin, param.displacement);
      write(bin, param.mom);
      write(bin, param.mass_label);
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


    //----------------------------------------------------------------------------
    //! Normalize just one displacement array
    multi1d<int> normDisp(const multi1d<int>& orig)
    {
      START_CODE();

      multi1d<int> disp;
      multi1d<int> empty; 
      multi1d<int> no_disp(1); no_disp[0] = 0;

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


    //-------------------------------------------------------------------------------
    // Function call
    void 
    InlineMeas::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
      // If xml file not empty, then use alternate
      if (params.xml_file != "")
      {
	string xml_file = makeXMLFileName(params.xml_file, update_no);

	push(xml_out, "GenPropMatElemDAColorVec");
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


    // Function call
    void 
    InlineMeas::func(unsigned long update_no,
		     XMLWriter& xml_out) 
    {
      START_CODE();

      StopWatch snoop;
      snoop.reset();
      snoop.start();

      StopWatch swiss;
			
      push(xml_out, "GenPropMatElemDAColorVec");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << name << ": Generalized propagator color-vector matrix element" << endl;

      // Test and grab a reference to the gauge field
      XMLBufferWriter gauge_xml, sm_gauge_xml;
      XMLReader source_prop_file_xml, source_prop_record_xml;
      XMLReader sink_prop_file_xml, sink_prop_record_xml;
      try
      {
	TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);
	TheNamedObjMap::Instance().get(params.named_obj.gauge_id).getRecordXML(gauge_xml);

	TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.sm_gauge_id);
	TheNamedObjMap::Instance().get(params.named_obj.sm_gauge_id).getRecordXML(sm_gauge_xml);

	*(TheNamedObjMap::Instance().getData< Handle<QDP::MapObject<KeyPropColorVec_t,LatticeFermion> > >(params.named_obj.source_prop_id));

	*(TheNamedObjMap::Instance().getData< Handle<QDP::MapObject<KeyPropColorVec_t,LatticeFermion> > >(params.named_obj.sink_prop_id));

	// Snarf the prop info. This is will throw if the prop_id is not there
	TheNamedObjMap::Instance().get(params.named_obj.source_prop_id).getFileXML(source_prop_file_xml);
	TheNamedObjMap::Instance().get(params.named_obj.source_prop_id).getRecordXML(source_prop_record_xml);

	// Snarf the prop info. This is will throw if the prop_id is not there
	TheNamedObjMap::Instance().get(params.named_obj.sink_prop_id).getFileXML(sink_prop_file_xml);
	TheNamedObjMap::Instance().get(params.named_obj.sink_prop_id).getRecordXML(sink_prop_record_xml);
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
      // Cast should be valid now
      const multi1d<LatticeColorMatrix>& u = 
	TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

      // Cast should be valid now
      const multi1d<LatticeColorMatrix>& u_smr = 
	TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.sm_gauge_id);

      QDP::MapObject<KeyPropColorVec_t,LatticeFermion>& source_ferm_map =
	*(TheNamedObjMap::Instance().getData< Handle<QDP::MapObject<KeyPropColorVec_t,LatticeFermion> > >(params.named_obj.source_prop_id));

      QDP::MapObject<KeyPropColorVec_t,LatticeFermion>& sink_ferm_map =
	*(TheNamedObjMap::Instance().getData< Handle<QDP::MapObject<KeyPropColorVec_t,LatticeFermion> > >(params.named_obj.sink_prop_id));

      push(xml_out, "Output_version");
      write(xml_out, "out_version", 2);
      pop(xml_out);

      proginfo(xml_out);    // Print out basic program info

      // Write out the input
      params.writeXML(xml_out, "Input");

      // Write out the source header
      write(xml_out, "Source_prop_file_info", source_prop_file_xml);
      write(xml_out, "Source_prop_record_info", source_prop_record_xml);
      write(xml_out, "Sink_prop_file_info", sink_prop_file_xml);
      write(xml_out, "Sink_prop_record_info", sink_prop_record_xml);

      // Write out the config info
      write(xml_out, "Config_info", gauge_xml);
      write(xml_out, "Smeared_Config_info", sm_gauge_xml);

      //First calculate some gauge invariant observables just for info.
      //This is really cheap.
      MesPlq(xml_out, "Observables", u);
      MesPlq(xml_out, "Smeared_Observables", u_smr);

      //
      // Initialize the slow Fourier transform phases
      //
      multi1d<int> origin_offset(Nd);
      origin_offset = 0;
      multi1d<int> mom_offset(Nd);
      mom_offset = 0;
      SftMom phases(0, origin_offset, mom_offset,false, params.param.decay_dir);

      //
      // DB storage
      //
      BinaryStoreDB< SerialDBKey<KeyGenPropElementalOperator_t>, SerialDBData<ValGenPropElementalOperator_t> > 
	qdp_db;

      // Open the file, and write the meta-data and the binary for this operator
      if (! qdp_db.fileExists(params.named_obj.genprop_op_file))
      {
	XMLBufferWriter file_xml;

	push(file_xml, "DBMetaData");
	write(file_xml, "id", string("genPropElemOp"));
	write(file_xml, "lattSize", QDP::Layout::lattSize());
	write(file_xml, "decay_dir", params.param.decay_dir);
	proginfo(file_xml);    // Print out basic program info
	write(file_xml, "Params", params.param);
	write(file_xml, "Op_Info", params.param.gamma);
	write(file_xml, "Source_prop_record_info", source_prop_record_xml);
	write(file_xml, "Sink_prop_record_info", sink_prop_record_xml);
	write(file_xml, "Config_info", gauge_xml);
	pop(file_xml);

	// NOTE: we always write a metadata, even if the file exists. 
	// Probably should not do this if a file already exists.
	// Yukky, but add on a bunch of padding. 
	// This is in case a latter task writes another metadata.
	std::string file_str(file_xml.str());
	qdp_db.setMaxUserInfoLen(file_str.size());

	qdp_db.open(params.named_obj.genprop_op_file, O_RDWR | O_CREAT, 0664);

	qdp_db.insertUserdata(file_str);
      }
      else
      {
	qdp_db.open(params.named_obj.genprop_op_file, O_RDWR, 0664);
      }


      //
      // Generalized propagatos
      //
      QDPIO::cout << "Building generalized propagators" << endl;

      push(xml_out, "ElementalOps");

      // Loop over all time slices for the source. This is the same 
      // as the subsets for  phases

      const bool restrict_plateau = params.param.restrict_plateau;
      const int num_vecs          = params.param.num_vecs;
      const int decay_dir         = params.param.decay_dir;
      const int t_source          = params.param.t_source;
      const int t_sink            = params.param.t_sink;
      const int z_dir             = params.param.boost_dir ;
      const int gamma             = params.param.gamma ;

      int Nt = phases.numSubsets();
      // Define the start and end region for the plateau
      int t_start;
      int t_end;

      if (restrict_plateau)
      {
	t_start = t_source;
	t_end   = t_sink;

	if (t_source > t_sink)
	{
	  t_end += phases.numSubsets();
	}
      }
      else
      {
	t_start = 0;
	t_end   = phases.numSubsets() - 1;
      }

      // Loop over momenta 
      for(int p=0; p < params.param.mom_max; ++p)
      {
	StopWatch watch;
	multi1d<int> mom(3); mom=0 ; mom[z_dir]=p;
	Double phi = (6.28318530717958647696*p)/toDouble(Layout::lattSize()[z_dir]);
	QDPIO::cout << "Boost Momentum: p = " << p << endl;

	swiss.reset();
	swiss.start();

	// Loop over spins
	for(int spin_r=0; spin_r < Ns; ++spin_r){
	  QDPIO::cout << "spin_r = " << spin_r << endl; 
	  for(int spin_l=0; spin_l < Ns; ++spin_l){
	    QDPIO::cout << "spin_l = " << spin_l << endl; 
	    // The keys for the spin and displacements 
	    // for this particular elemental operator
	    // No displacement for left colorvector, 
	    // only displace right colorvector
	    
	    // Invert the time - make it an independent key
	    multi2d<KeyValGenPropElementalOperator_t> buf(Layout::lattSize()[z_dir],phases.numSubsets());
	    for( int z(0);z<Layout::lattSize()[z_dir];z++){
	      multi1d<int> disp(Nd-1); disp=0 ; disp[z_dir]=z;
	      for(int tt=t_start; tt <= t_end; ++tt){
		// mod back into a normal interval
		int t = tt % phases.numSubsets(); 
		buf(z,t).key.key().t_slice       = t;
		buf(z,t).key.key().t_source      = t_source;
		buf(z,t).key.key().t_sink        = t_sink;
		buf(z,t).key.key().spin_r        = spin_r;
		buf(z,t).key.key().spin_l        = spin_l;
		buf(z,t).key.key().mass_label    = params.param.mass_label;
		buf(z,t).key.key().mom           = mom;
		buf(z,t).key.key().gamma         = gamma;
		buf(z,t).key.key().displacement  = disp; 
		buf(z,t).val.data().op.resize(num_vecs, num_vecs);
	      }
	    }
	    for(int j = 0; j < params.param.num_vecs; ++j){
	      KeyPropColorVec_t key_r;
	      key_r.t_source     = t_source;
	      key_r.colorvec_src = j;
	      key_r.spin_src     = spin_r;
	      
	      // Displace the right vector and multiply by the momentum phase
	      
	      LatticeFermion shift_ferm;
	      LatticeFermion src_ferm; source_ferm_map.get(key_r, src_ferm);
	      for( int z(0);z<Layout::lattSize()[z_dir];z++){
		multi1d<int> disp(Nd-1); disp=0 ; disp[z_dir]=z;
		if(z==0)
		  shift_ferm=src_ferm ;
		else{
		  LatticeFermion tt = shift(shift_ferm,BACKWARD,z_dir);
		  shift_ferm=u_smr[z_dir]*tt;
		}
		for(int i = 0; i < params.param.num_vecs; ++i){
		  KeyPropColorVec_t key_l;
		  key_l.t_source     = t_sink;
		  key_l.colorvec_src = i;
		  key_l.spin_src     = spin_l;
		  
		  watch.reset();
		  watch.start();
		  
		  // Contract over color indices
		  // Do the relevant quark contraction
		  LatticeFermion tt; sink_ferm_map.get(key_l, tt);
		  
		  // I need to rember to fold the Gamma_5 on 
		  // the source into the interpolating field
		  LatticeFermion snk=Gamma(Nd*Nd-1)*tt;
		  LatticeComplex lop = localInnerProduct(snk, shift_ferm);
		  LatticeReal foo = phi*Layout::latticeCoordinate(z_dir);
		  LatticeComplex phase = cmplx(cos(foo),sin(foo));
		  // Slow fourier-transform
		  multi1d<ComplexD> op_sum=sumMulti(phase*lop,phases.getSet());
		  watch.stop();
		  for(int tt=t_start; tt <= t_end; ++tt){
		    // mod back into a normal interval
		    int t = tt % phases.numSubsets(); 
		    buf(z,t).val.data().op(i,j) = op_sum[t];
		  }// loop over time
		}// loop over left (i)
	      }// loop over displacements 
	    }//loop over right (j)
	    QDPIO::cout << "insert: mom= " << mom << endl ;
	    for( int z(0);z<Layout::lattSize()[z_dir];z++){
	      for(int tt=t_start; tt <= t_end; ++tt){
		// mod back into a normal interval
		int t = tt % phases.numSubsets(); 
		qdp_db.insert(buf(z,t).key, buf(z,t).val);
	      } // loop over tt
	    }// loop over z
	  } // end for spin_l
	} // end for spin_r
	swiss.stop();
	QDPIO::cout << "GenProp mom= "<<p  
		    << "  time= "
		    << swiss.getTimeInSeconds() 
		    << " secs" << endl;
      } // loop over momenta
      
      
      pop(xml_out); // ElementalOps

      // Close the namelist output file XMLDAT
      pop(xml_out);     // GenPropMatElemDAColorVector

      snoop.stop();
      QDPIO::cout << name << ": total time = " 
		  << snoop.getTimeInSeconds() 
		  << " secs" << endl;

      QDPIO::cout << name << ": ran successfully" << endl;

      END_CODE();
    } // func
  } // namespace InlineGenPropMatElemDAColorVecEnv 

  /*! @} */  // end of group hadron

} // namespace Chroma
