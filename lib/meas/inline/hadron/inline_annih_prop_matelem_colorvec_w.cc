// $Id: inline_annih_prop_matelem_colorvec_w.cc,v 3.4 2009-09-14 20:50:14 edwards Exp $
/*! \file
 * \brief Compute the annihilation diagram propagator elements    M^-1 * multi1d<LatticeColorVector>
 *
 * Annihilation diagrams version of propagator calculation on a colorvector
 */

#include "fermact.h"
#include "meas/inline/hadron/inline_annih_prop_matelem_colorvec_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "meas/sources/zN_src.h"
#include "util/ferm/subset_vectors.h"
#include "qdp_map_obj_memory.h"
#include "util/ferm/key_val_db.h"
#include "util/ferm/key_prop_colorvec.h"
#include "util/ferm/key_prop_matelem.h"
#include "util/ferm/transf.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"
#include "meas/inline/make_xml_file.h"

#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{ 
  namespace InlineAnnihPropMatElemColorVecEnv 
  {
    //! Propagator input
    void read(XMLReader& xml, const string& path, Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_id", input.gauge_id);
      read(inputtop, "colorvec_id", input.colorvec_id);
      read(inputtop, "prop_op_file", input.prop_op_file);
    }

    //! Propagator output
    void write(XMLWriter& xml, const string& path, const Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "gauge_id", input.gauge_id);
      write(xml, "colorvec_id", input.colorvec_id);
      write(xml, "prop_op_file", input.prop_op_file);

      pop(xml);
    }


    //! Propagator input
    void read(XMLReader& xml, const string& path, Params::Param_t::Contract_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "num_vecs", input.num_vecs);
      read(inputtop, "t_sources_start", input.t_sources_start);
      read(inputtop, "dt", input.dt);
      read(inputtop, "decay_dir", input.decay_dir);
      read(inputtop, "N", input.N);
      read(inputtop, "ran_seed", input.ran_seed);
      read(inputtop, "mass_label", input.mass_label);
    }

    //! Propagator output
    void write(XMLWriter& xml, const string& path, const Params::Param_t::Contract_t& input)
    {
      push(xml, path);

      write(xml, "num_vecs", input.num_vecs);
      write(xml, "t_sources_start", input.t_sources_start);
      write(xml, "dt", input.dt);
      write(xml, "decay_dir", input.decay_dir);
      write(xml, "N", input.N);
      write(xml, "ran_seed", input.ran_seed);
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
  } // namespace InlinePropColorVecEnv 


  namespace InlineAnnihPropMatElemColorVecEnv 
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
      
    const std::string name = "ANNIH_PROP_MATELEM_COLORVEC";

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



    // Function call
    void 
    InlineMeas::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
      // If xml file not empty, then use alternate
      if (params.xml_file != "")
      {
	string xml_file = makeXMLFileName(params.xml_file, update_no);

	push(xml_out, "AnnihPropMatElemColorVec");
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

      push(xml_out, "AnnihPropMatElemColorVec");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << name << ": annihilation propagator matrix element calculation" << endl;

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
      XMLReader source_file_xml, source_record_xml;

      QDPIO::cout << "Snarf the source from a named buffer" << endl;
      try
      {
	TheNamedObjMap::Instance().getData< Handle< QDP::MapObject<int,EVPair<LatticeColorVector> > > >(params.named_obj.colorvec_id);

	// Snarf the source info. This is will throw if the colorvec_id is not there
	TheNamedObjMap::Instance().get(params.named_obj.colorvec_id).getFileXML(source_file_xml);
	TheNamedObjMap::Instance().get(params.named_obj.colorvec_id).getRecordXML(source_record_xml);

	// Write out the source header
	write(xml_out, "Source_file_info", source_file_xml);
	write(xml_out, "Source_record_info", source_record_xml);
      }    
      catch (std::bad_cast)
      {
	QDPIO::cerr << name << ": caught dynamic cast error" << endl;
	QDP_abort(1);
      }
      catch (const string& e) 
      {
	QDPIO::cerr << name << ": error extracting source_header: " << e << endl;
	QDP_abort(1);
      }

      // Cast should be valid now
      const QDP::MapObject<int,EVPair<LatticeColorVector> >& eigen_source = 
	*(TheNamedObjMap::Instance().getData< Handle< QDP::MapObject<int,EVPair<LatticeColorVector> > > >(params.named_obj.colorvec_id));

      QDPIO::cout << "Source successfully read and parsed" << endl;

      // Another sanity check
      if (params.param.contract.num_vecs > eigen_source.size())
      {
	QDPIO::cerr << __func__ << ": num_vecs= " << params.param.contract.num_vecs
		    << " is greater than the number of available colorvectors= "
		    << eigen_source.size() << endl;
	QDP_abort(1);
      }

      // Interval of sources should divide into lattice extent
      if ((QDP::Layout::lattSize()[params.param.contract.decay_dir] % params.param.contract.dt) != 0)
      {
	QDPIO::cerr << __func__ << ": dt= " << params.param.contract.dt
		    << " does not divide into lattice extent" << endl;
	QDP_abort(1);
      }


      //
      // DB storage
      //
      BinaryStoreDB< SerialDBKey<KeyPropElementalOperator_t>, SerialDBData<ValPropElementalOperator_t> > qdp_db;

      // Open the file, and write the meta-data and the binary for this operator
      if (! qdp_db.fileExists(params.named_obj.prop_op_file))
      {
	XMLBufferWriter file_xml;

	push(file_xml, "DBMetaData");
	write(file_xml, "id", string("propElemOp"));
	write(file_xml, "lattSize", QDP::Layout::lattSize());
	write(file_xml, "decay_dir", params.param.contract.decay_dir);
	proginfo(file_xml);    // Print out basic program info
	write(file_xml, "Params", params.param.contract);
	write(file_xml, "Config_info", gauge_xml);
	write(file_xml, "Weights", getEigenValues(eigen_source, params.param.contract.num_vecs));
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



      // Save current seed
      Seed ran_seed;
      QDP::RNG::savern(ran_seed);
	
      // Set the seed to desired value
      QDP::RNG::setrn(params.param.contract.ran_seed);

      // Total number of iterations
      int ncg_had = 0;

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
	const int num_vecs        = params.param.contract.num_vecs;
	const int decay_dir       = params.param.contract.decay_dir;
	const int dt              = params.param.contract.dt;
	const int Lt              = QDP::Layout::lattSize()[decay_dir];
	const int num_sources     = Lt / params.param.contract.dt;
	const multi1d<int>& t_sources_start = params.param.contract.t_sources_start;

	// Initialize the slow Fourier transform phases
	SftMom phases(0, true, decay_dir);

	for(int tt=0; tt < t_sources_start.size(); ++tt)
	{
	  int t_source_start = params.param.contract.t_sources_start[tt];
	  // Construct all the solution vectors for only the first time source
	  QDPIO::cout << "t_source_start = " << t_source_start << endl; 

	  //
	  // Have to hold temporarily the solution vectors
	  //
	  MapObjectMemory<KeyPropColorVec_t,LatticeFermion> map_obj;

	  // The random numbers used for the stochastic combination of sources is held here
	  MapObjectMemory<int,Complex> rng_map_obj;

	  // Build up list of time sources
	  multi1d<int> t_sources(num_sources);

	  for(int nn=0; nn < t_sources.size(); ++nn)
	  {
	    int t = (t_source_start + dt*nn + Lt) % Lt;
	    t_sources[nn] = t;

	    rng_map_obj.insert(t, zN_rng(params.param.contract.N));
	  }

	  for(int colorvec_source=0; colorvec_source < num_vecs; ++colorvec_source) {
	    
	    QDPIO::cout << "colorvec_source = " << colorvec_source << endl; 

	    // Accumulate the source from several time-slices
	    LatticeColorVector vec_srce = zero;
	    for(int nn=0; nn < t_sources.size(); ++nn) {
	      int t = t_sources[nn];

	      // Pull out a time-slice of the color vector source, give it a random weight

	      Complex weight; rng_map_obj.get(t, weight);
	      EVPair<LatticeColorVector> tmpvec; eigen_source.get(colorvec_source, tmpvec);
	      vec_srce[phases.getSet()[t]] = weight * tmpvec.eigenVector;
	    }
	
	    for(int spin_source=0; spin_source < Ns; ++spin_source) {

	      QDPIO::cout << "spin_source = " << spin_source << endl; 

	      // Insert a ColorVector into spin index spin_source
	      // This only overwrites sections, so need to initialize first
	      LatticeFermion chi = zero;
	      CvToFerm(vec_srce, chi, spin_source);

	      LatticeFermion quark_soln = zero;

	      // Do the propagator inversion
	      SystemSolverResults_t res = (*PP)(quark_soln, chi);
	      ncg_had = res.n_count;

	      KeyPropColorVec_t key;
	      key.t_source     = t_source_start;
	      key.colorvec_src = colorvec_source;
	      key.spin_src     = spin_source;
	    
	      map_obj.insert(key, quark_soln);
	    } // for spin_source
	  } // for colorvec_source


	  swatch.stop();
	  QDPIO::cout << "Propagators computed: time= " 
		      << swatch.getTimeInSeconds() 
		      << " secs" << endl;


	  //
	  // All the loops
	  //
	  // NOTE: I pull the spin source and sink loops to the outside intentionally.
	  // The idea is to create a colorvector index (2d) array. These are not
	  // too big, but are big enough to make the IO efficient, and the DB efficient
	  // on reading. For N=32 and Lt=128, the mats are 2MB.
	  //
	  swatch.reset();
	  swatch.start();
	  QDPIO::cout << "Extract matrix elements" << endl;

	  for(int spin_source=0; spin_source < Ns; ++spin_source) {

	    QDPIO::cout << "spin_source = " << spin_source << endl; 
	  
	    for(int spin_sink=0; spin_sink < Ns; ++spin_sink) {

	      QDPIO::cout << "spin_sink = " << spin_sink << endl; 

	      // Construct the keys and values. Have to hold the buffers here given that
	      // the source and sink spin indices are pulled out as well as the time
	      multi1d<KeyValPropElementalOperator_t> buf(t_sources.size());
	      for(int nn=0; nn < t_sources.size(); ++nn)
	      {
		int t = t_sources[nn];

		buf[nn].key.key().t_slice      = t;
		buf[nn].key.key().t_source     = t;
		buf[nn].key.key().spin_src     = spin_source;
		buf[nn].key.key().spin_snk     = spin_sink;
		buf[nn].key.key().mass_label   = params.param.contract.mass_label;
		buf[nn].val.data().mat.resize(num_vecs,num_vecs);
	      }

	      for(int colorvec_source=0; colorvec_source < num_vecs; ++colorvec_source)
	      {
		KeyPropColorVec_t key;
		key.t_source     = t_source_start;
		key.colorvec_src = colorvec_source;
		key.spin_src     = spin_source;

		LatticeColorVector vec_source;
		{
		  LatticeFermion tmp;
		  QDPIO::cout << "MAP_OBJ LOOKUP" << endl << flush;
		  map_obj.get(key,tmp);
		  vec_source=peekSpin(tmp, spin_sink);
		}

		for(int colorvec_sink=0; colorvec_sink < num_vecs; ++colorvec_sink)
		{
		  EVPair<LatticeColorVector> vec_sink; eigen_source.get(colorvec_sink, vec_sink);

		  multi1d<ComplexD> hsum(sumMulti(localInnerProduct(vec_sink.eigenVector, vec_source), phases.getSet()));
		
		  for(int nn=0; nn < t_sources.size(); ++nn)
		  {
		    int t = t_sources[nn];
		    Complex weight; rng_map_obj.get(t,weight);
		    buf[nn].val.data().mat(colorvec_sink,colorvec_source) = conj(weight) * hsum[t];
		  }

		} // for colorvec_sink
	      } // for colorvec_source
	      
	      QDPIO::cout << "insert: spin_source= " << spin_source << " spin_sink= " << spin_sink << endl; 
	      for(int nn=0; nn < t_sources.size(); ++nn)
	      {
		int t = t_sources[nn];
		qdp_db.insert(buf[nn].key, buf[nn].val);
	      }
	      

	    } // for spin_sink
	  } // for spin_source

	  swatch.stop();
	  QDPIO::cout << "Matrix elements computed: time= " 
		      << swatch.getTimeInSeconds() 
		      << " secs" << endl;
	} // tt

      }
      catch (const std::string& e) 
      {
	QDPIO::cout << name << ": caught exception around qprop: " << e << endl;
	QDP_abort(1);
      }

      push(xml_out,"Relaxation_Iterations");
      write(xml_out, "ncg_had", ncg_had);
      pop(xml_out);

      pop(xml_out);  // prop_colorvec

      // Reset the seed
      QDP::RNG::setrn(ran_seed);

      snoop.stop();
      QDPIO::cout << name << ": total time = "
		  << snoop.getTimeInSeconds() 
		  << " secs" << endl;

      QDPIO::cout << name << ": ran successfully" << endl;

      END_CODE();
    } 

  }

} // namespace Chroma
