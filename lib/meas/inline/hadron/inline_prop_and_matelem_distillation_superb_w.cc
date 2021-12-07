/*! \file
 * \brief Compute propagators from distillation
 *
 * Propagator calculation in distillation
 */

#include "qdp.h"
#include "fermact.h"
#include "meas/inline/hadron/inline_prop_and_matelem_distillation_superb_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "qdp_map_obj.h"
#include "qdp_map_obj_disk.h"
#include "qdp_map_obj_disk_multiple.h"
#include "qdp_map_obj_memory.h"
#include "qdp_disk_map_slice.h"
#include "util/ferm/subset_vectors.h"
#include "util/ferm/key_prop_distillation.h"
#include "util/ferm/key_prop_colorvec.h"
#include "util/ferm/key_prop_matelem.h"
#include "util/ferm/key_val_db.h"
#include "util/ferm/transf.h"
#include "util/ferm/spin_rep.h"
#include "util/ferm/diractodr.h"
#include "util/ferm/twoquark_contract_ops.h"
#include "util/ferm/superb_contractions.h"
#include "util/ft/sftmom.h"
#include "util/ft/time_slice_set.h"
#include "util/info/proginfo.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"
#include "meas/inline/make_xml_file.h"

#include "meas/inline/io/named_objmap.h"

#include "chroma_config.h"

#ifdef BUILD_SB

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
  }
  

  //----------------------------------------------------------------------------
  namespace InlinePropAndMatElemDistillationSuperbEnv 
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
      read(inputtop, "Nt_backward", input.Nt_backward);
      read(inputtop, "mass_label", input.mass_label);
      read(inputtop, "num_tries", input.num_tries);

      input.max_rhs = 8;
      if( inputtop.count("max_rhs") == 1 ) {
        read(inputtop, "max_rhs", input.max_rhs);
      }

      input.phase.resize(Nd - 1);
      for (int i = 0; i < Nd - 1; ++i)
	input.phase[i] = 0;
      if( inputtop.count("phase") == 1 ) {
        read(inputtop, "phase", input.phase);
      }

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
      write(xml, "num_tries", input.num_tries);
      write(xml, "max_rhs", input.max_rhs);
      write(xml, "phase", input.phase);

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
  } // namespace InlinePropDistillationSuperbEnv 


  //----------------------------------------------------------------------------
  namespace InlinePropAndMatElemDistillationSuperbEnv
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
      
    const std::string name = "PROP_AND_MATELEM_DISTILLATION_SUPERB";

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

      SB::ColorvecsStorage colorvecsSto = SB::openColorvecStorage(params.named_obj.colorvec_files);
     
      //
      // DB storage
      //
      BinaryStoreDB< SerialDBKey<KeyPropElementalOperator_t>, SerialDBData<ValPropElementalOperator_t> > qdp_db;

      if (!params.param.contract.zero_colorvecs)
	{
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
	}

      //
      // Parse the phase
      //
      if (params.param.contract.phase.size() != Nd - 1)
      {
	QDPIO::cerr << "phase tag should have " << Nd - 1 << " components" << std::endl;
	QDP_abort(1);
      }
      SB::Coor<Nd - 1> phase;
      for (int i = 0; i < Nd - 1; ++i)
      {
	phase[i] = params.param.contract.phase[i];
	if (std::fabs(phase[i] - params.param.contract.phase[i]) > 0)
	  std::runtime_error("phase should be integer");
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
	const int max_rhs             = params.param.contract.max_rhs;

	// Loop over each time source
	for (int tt = 0; tt < t_sources.size(); ++tt)
	{
	  int t_source = t_sources[tt]; // This is the actual time-slice.
	  QDPIO::cout << "t_source = " << t_source << std::endl;

	  // Compute the first tslice and the number of tslices involved in the contraction
	  int first_tslice = t_source - params.param.contract.Nt_backward;
	  int num_tslices = std::min(
	    params.param.contract.Nt_backward + std::max(1, params.param.contract.Nt_forward), Lt);

	  // Get `num_vecs` colorvecs, and `num_tslices` tslices starting from time-slice `first_tslice`
	  SB::Tensor<Nd + 3, SB::Complex> colorvec = SB::getColorvecs<SB::Complex>(
	    colorvecsSto, u, decay_dir, first_tslice, num_tslices, num_vecs, "cxyzXnt", phase);

	  // Get all eigenvectors for `t_source`
	  auto source_colorvec =
	    colorvec.kvslice_from_size({{'t', t_source - first_tslice}}, {{'t', 1}});

	  for (int spin_source = 0; spin_source < Ns; ++spin_source)
	  {
	    // Invert the source for `spin_source` spin and retrieve `num_tslices` tslices starting from tslice `first_tslice`
	    // NOTE: s is spin source, and S is spin sink
	    SB::Tensor<Nd + 5, SB::Complex> quark_solns = SB::doInversion<SB::Complex, SB::Complex>(
	      *PP, source_colorvec, t_source, first_tslice, num_tslices, {spin_source}, max_rhs,
	      "cxyzXnSst");

	    StopWatch snarss1;
	    snarss1.reset();
	    snarss1.start();

	    // Contract the distillation elements
	    // NOTE: N: is colorvec in sink, and n is colorvec in source
	    SB::Tensor<5, SB::Complex> elems("NnSst", {num_vecs, num_vecs, Ns, 1, num_tslices},
					     SB::OnHost, SB::OnMaster);
	    elems.contract(colorvec, {{'n', 'N'}}, SB::Conjugate, quark_solns, {},
			   SB::NotConjugate);

	    snarss1.stop();
	    QDPIO::cout << "Time to contract for one spin source : " << snarss1.getTimeInSeconds()
			<< " secs" << std::endl;

	    snarss1.reset();
	    snarss1.start();

	    ValPropElementalOperator_t val;
	    val.mat.resize(num_vecs, num_vecs);
	    val.mat = zero;
	    for (int i_tslice = 0; i_tslice < num_tslices; ++i_tslice)
	    {
	      for (int spin_sink = 0; spin_sink < Ns; ++spin_sink)
	      {
		KeyPropElementalOperator_t key;
		key.t_source = t_source;
		key.t_slice = SB::normalize_coor(i_tslice + first_tslice, Lt);
		key.spin_src = spin_source;
		key.spin_snk = spin_sink;
		key.mass_label = params.param.contract.mass_label;
		if (Layout::nodeNumber() == 0)
		{
		  for (int colorvec_sink = 0; colorvec_sink < num_vecs; ++colorvec_sink)
		  {
		    for (int colorvec_source = 0; colorvec_source < num_vecs; ++colorvec_source)
		    {
		      std::complex<REAL> e = elems.get(
			{colorvec_sink, colorvec_source, spin_sink, 0, i_tslice});
		      val.mat(colorvec_sink, colorvec_source).elem().elem().elem() =
			RComplex<REAL64>(e.real(), e.imag());
		    }
		  }
		}
		qdp_db.insert(key, val);
	      }
	    }

	    snarss1.stop();
	    QDPIO::cout << "Time to store the props : " << snarss1.getTimeInSeconds() << " secs"
			<< std::endl;
	  } // for spin_source
	}   // for tt

	swatch.stop();
	QDPIO::cout << "Propagators computed: time= " 
		    << swatch.getTimeInSeconds() 
		    << " secs" << std::endl;
      }
      catch (const std::exception& e) 
      {
	QDP_error_exit("%s: caught exception: %s\n", name, e.what());
      }

      // Close colorvecs storage
      SB::closeColorvecStorage(colorvecsSto);

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

#endif // BUILD_SB
