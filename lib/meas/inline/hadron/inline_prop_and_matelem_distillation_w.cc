/*! \file
 * \brief Compute propagators from distillation
 *
 * Propagator calculation in distillation
 */

#include "qdp.h"
#include "fermact.h"
#include "meas/inline/hadron/inline_prop_and_matelem_distillation_w.h"
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
#include "util/ferm/superb_contractions.h"
#include "util/ft/sftmom.h"
#include "util/ft/time_slice_set.h"
#include "util/info/proginfo.h"
#include "util/info/proginfo.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"
#include "meas/inline/make_xml_file.h"

#include "meas/inline/io/named_objmap.h"

#include "chroma_config.h"

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
  namespace InlinePropAndMatElemDistillationEnv 
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
  namespace InlinePropAndMatElemDistillationEnv
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
      
    const std::string name = "PROP_AND_MATELEM_DISTILLATION";

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
      QDPIO::cout << "Snarf the source from a std::map object disk file" << std::endl;

      SB::MODS_t eigen_source;
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
      // QDPIO::cout << "Initialize sub-lattice std::map" << std::endl;
      // SubEigenMap sub_eigen_map(eigen_source, decay_dir, params.param.contract.zero_colorvecs);
      // QDPIO::cout << "Finished initializing sub-lattice std::map" << std::endl;


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
	      if (!params.param.contract.zero_colorvecs)
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

	  // The tensor `active_colorvec` will have the eigenvectors for the active time-slides;
	  // the contraction requires `t` being the slowest index
	  int first_tslice_in_active_colorvec = t_source - params.param.contract.Nt_backward;
	  int num_tslices_in_active_colorvec =
	    std::min(params.param.contract.Nt_backward + params.param.contract.Nt_forward + 1, Lt);
	  const char order_in_active_colorvec[] = "cxyznt";
	  SB::Tensor<Nd + 2, SB::ComplexF> active_colorvec =
	    SB::getColorvecs(eigen_source, decay_dir, first_tslice_in_active_colorvec,
			     num_tslices_in_active_colorvec, num_vecs, order_in_active_colorvec);

	  // The tensor `tensor_quark_solns` will have the active time-slices of the solution vectors;
	  // the contraction requires `t` and `s` being the slowest indices
	  const char order_in_tensor_quark_solns[] = "cxyznst";
	  SB::Tensor<Nd + 3, SB::ComplexF> tensor_quark_solns(
	    order_in_tensor_quark_solns,
	    SB::latticeSize<Nd + 3>(order_in_tensor_quark_solns,
				    {{'t', num_tslices_in_active_colorvec}, {'n', num_vecs}}));

	  // Loop over each spin source
	  for (int spin_source = 0; spin_source < Ns; ++spin_source)
	  {
	    QDPIO::cout << "spin_source = " << spin_source << std::endl;

	    //
	    // The space distillation loop
	    //
	    for (int colorvec_src0 = 0, colorvec_src_step = std::min(max_rhs, num_vecs);
		 colorvec_src0 < num_vecs; colorvec_src0 += colorvec_src_step,
		     colorvec_src_step = std::min(colorvec_src_step, num_vecs - colorvec_src0))
	    {
	      std::vector<std::shared_ptr<LatticeFermion>> chis(colorvec_src_step),
		quark_solns(colorvec_src_step);
	      for (int col = 0; col < colorvec_src_step; col++)
		chis[col].reset(new LatticeFermion);
	      for (int col = 0; col < colorvec_src_step; col++)
		quark_solns[col].reset(new LatticeFermion);

	      StopWatch snarss1;
	      snarss1.reset();
	      snarss1.start();

	      for (int colorvec_src = colorvec_src0, col = 0; col < colorvec_src_step;
		   ++colorvec_src, ++col)
	      {
		QDPIO::cout << "Do spin_source= " << spin_source
			    << "  colorvec_src= " << colorvec_src << std::endl;

		// Put the colorvec sources for the t_source on chis for spin `spin_source`
		// chis[col][s=spin_source] = active_colorvec[t=t_source-first_tslice_in_active_colorvec, n=colorvec_src0]
		*chis[col] = zero;
		active_colorvec
		  .kvslice_from_size(
		    {{'t', t_source - first_tslice_in_active_colorvec}, {'n', colorvec_src0}},
		    {{'t', 1}, {'n', 1}})
		  .copyTo(SB::asTensorView(*chis[col]).kvslice_from_size({{'s', spin_source}}));

		*quark_solns[col] = zero;
	      } // colorvec_src

	      if (!params.param.contract.zero_colorvecs)
	      {
		// Do the propagator inversion
		// Check if bad things are happening
		bool badP = true;
		for (int nn = 1; nn <= params.param.contract.num_tries; ++nn)
		{
		  // Reset
		  for (int col = 0; col < colorvec_src_step; col++)
		    *quark_solns[col] = zero;
		  badP = false;

		  // Solve for the solution std::vector
		  std::vector<SystemSolverResults_t> res = PP->operator()(
		    quark_solns,
		    std::vector<std::shared_ptr<const LatticeFermion>>(chis.begin(), chis.end()));

		  for (int col = 0; col < colorvec_src_step; col++)
		  {
		    ncg_had += res[col].n_count;

		    // Some sanity checks
		    if (toDouble(res[col].resid) > 1.0e-3)
		    {
		      QDPIO::cerr << name << ": have a resid > 1.0e-3. That is unacceptable"
				  << std::endl;
		      QDP_abort(1);
		    }
		  }

		  // Check for finite values - neither NaN nor Inf
		  for (int col = 0; col < colorvec_src_step; col++)
		    if (isfinite(*quark_solns[col]))
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
		  QDPIO::cerr << name
			      << ": this is bad - did not get a finite solution std::vector "
				 "after num_tries= "
			      << params.param.contract.num_tries << std::endl;
		  QDP_abort(1);
		}
	      } // zero_colorvecs ??

	      snarss1.stop();
	      QDPIO::cout << "Time to compute prop for spin_source= " << spin_source
			  << "  colorvec_src= " << colorvec_src0 << " to "
			  << colorvec_src0 + colorvec_src_step - 1
			  << "  time = " << snarss1.getTimeInSeconds() << " secs" << std::endl;

	      for (int colorvec_src = colorvec_src0, col = 0; col < colorvec_src_step;
		   ++colorvec_src, ++col)
	      {
		// tensor_quark_solns[n=colorvec_src] = quark_solns[col][t=first_tslice_in_active_colorvec+(0:num_tslices_in_active_colorvec-1)]
		SB::asTensorView(*quark_solns[col])
		  .kvslice_from_size({{'t', first_tslice_in_active_colorvec}},
				     {{'t', num_tslices_in_active_colorvec}})
		  .copyTo(tensor_quark_solns.kvslice_from_size({{'n', colorvec_src}}));
	      }
	    }

	    // Contract the distillation elements
	    const char order_in_elems[] =
	      "nNst"; // N: colorvec in sink, n: colorvec in source, s: spin in sink
	    SB::Tensor<4, SB::ComplexF> elems(
	      order_in_elems, {num_vecs, num_vecs, Ns, num_tslices_in_active_colorvec}, SB::OnHost,
	      SB::OnMaster);
	    elems.contract(active_colorvec, {}, SB::Conjugate, tensor_quark_solns, {{'n', 'N'}},
			   SB::NotConjugate, {});

	    // Store them
	    if (!params.param.contract.zero_colorvecs)
	    {
	      QDPIO::cout << "Write all perambulator for spin_source= " << spin_source
			  << "  to disk" << std::endl;
	      ValPropElementalOperator_t val;
	      val.mat.resize(num_vecs, num_vecs);
	      val.mat = zero;
	      for (int t_slice = first_tslice_in_active_colorvec, i_tslice = 0;
		   i_tslice < num_tslices_in_active_colorvec;
		   ++i_tslice, t_slice = (t_slice + 1) % Lt)
	      {
		for (int spin_sink = 0; spin_sink < Ns; ++spin_sink)
		{
		  KeyPropElementalOperator_t key;
		  key.t_source = t_source;
		  key.t_slice = t_slice;
		  key.spin_src = spin_source;
		  key.spin_snk = spin_sink;
		  key.mass_label = params.param.contract.mass_label;
		  for (int colorvec_sink = 0; colorvec_sink < num_vecs; ++colorvec_sink)
		  {
		    for (int colorvec_source = 0; colorvec_source < num_vecs; ++colorvec_source)
		    {
		      std::complex<REAL> e =
			elems.get({colorvec_source, colorvec_sink, spin_sink, i_tslice});
		      val.mat(colorvec_sink, colorvec_source).elem().elem().elem() =
			RComplex<REAL64>(e.real(), e.imag());
		    }
		  }
		  qdp_db.insert(key, val);
		}
	      }
	    }
	  } // for spin_src
	}   // for tt

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
