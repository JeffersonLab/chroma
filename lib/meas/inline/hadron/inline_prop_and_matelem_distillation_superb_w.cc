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
#include "util/ferm/mgproton.h"
#include "util/ft/sftmom.h"
#include "util/ft/time_slice_set.h"
#include "util/info/proginfo.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"
#include "meas/inline/make_xml_file.h"

#include "meas/inline/io/named_objmap.h"

#include "chroma_config.h"

#include <set>

#ifdef BUILD_SB

namespace Chroma 
{ 

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
    void read(XMLReader& xml, const std::string& path, Params::Param_t::Phasing_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "sink", input.snk);
      if (input.snk.size() != Nd - 1)
      {
	throw std::runtime_error(std::string("phase sink should have ") + std::to_string(Nd - 1) +
				 " components");
      }
      read(inputtop, "source", input.src);
      if (input.src.size() != Nd - 1)
      {
	throw std::runtime_error(std::string("phase source should have ") + std::to_string(Nd - 1) +
				 " components");
      }
    }

    //! Propagator output
    void write(XMLWriter& xml, const std::string& path, const Params::Param_t::Phasing_t& input)
    {
      push(xml, path);

      write(xml, "sink", input.snk);
      write(xml, "source", input.src);

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

      input.max_rhs = 8;
      if( inputtop.count("max_rhs") == 1 ) {
        read(inputtop, "max_rhs", input.max_rhs);
      }

      if (inputtop.count("phases") == 1)
      {
	read(inputtop, "phases", input.phases);
      }
      else if (inputtop.count("phase") == 1)
      {
	read(inputtop, "phase", input.phase);
	if (input.phase.size() != Nd - 1)
	{
	  throw std::runtime_error(std::string("phase tag should have ") + std::to_string(Nd - 1) +
				   " components");
	}
      }
      else
      {
	input.phase.resize(Nd - 1);
	for (int i = 0; i < Nd - 1; ++i)
	  input.phase[i] = 0;
      }

      input.use_device_for_contractions = true;
      if( inputtop.count("use_device_for_contractions") == 1 ) {
        read(inputtop, "use_device_for_contractions", input.use_device_for_contractions);
      }

      input.use_superb_format = true;
      if( inputtop.count("use_superb_format") == 1 ) {
	read(inputtop, "use_superb_format", input.use_superb_format);
      }

      input.output_file_is_local = false;
      if( inputtop.count("output_file_is_local") == 1 ) {
        read(inputtop, "output_file_is_local", input.output_file_is_local);
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
      write(xml, "max_rhs", input.max_rhs);
      write(xml, "phase", input.phase);
      write(xml, "phases", input.phases);
      write(xml, "use_superb_format", input.use_superb_format);
      write(xml, "output_file_is_local", input.output_file_is_local);
      write(xml, "use_device_for_contractions", input.use_device_for_contractions);

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

      // Parse the phase
      std::vector<SB::Coor<3>> phasings;	     ///< list of all phasings
      std::map<int, std::set<int>> phasing_pairs;    ///< phasing src -> list of phasing sink
      if (params.param.contract.phase.size() > 0)
      {
	phasings.push_back(SB::toCoor(params.param.contract.phase));
	phasing_pairs[0] = std::set<int>{0};
      }
      else if (params.param.contract.phases.size() > 0)
      {
	if (!params.param.contract.use_superb_format)
	{
	  throw std::runtime_error(
	    "only the new format has support for `phases'; please set 'use_superb_format' to true");
	}
	for (const auto phasing_pair : params.param.contract.phases)
	{
	  auto src = SB::toCoor(phasing_pair.src);
	  auto snk_neg = SB::toCoor(phasing_pair.snk);
	  auto snk = SB::Coor<3>{-snk_neg[0], -snk_neg[1], -snk_neg[2]};
	  int src_idx =
	    std::find(phasings.begin(), phasings.end(), src) - phasings.begin();
	  if (src_idx == phasings.size())
	    phasings.push_back(src);
	  int snk_idx =
	    std::find(phasings.begin(), phasings.end(), snk) - phasings.begin();
	  if (snk_idx == phasings.size())
	    phasings.push_back(snk);
	  if (phasing_pairs.count(src_idx) == 0)
	  {
	    phasing_pairs[src_idx] = std::set<int>{snk_idx};
	  }
	  else
	  {
	    phasing_pairs[src_idx].insert(snk_idx);
	  }
	}
      }
      else
      {
	phasings.push_back(SB::Coor<Nd - 1>{{}});
	phasing_pairs[0] = std::set<int>{0};
      }

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

      //
      // Read in the source along with relevant information.
      // 

      SB::ColorvecsStorage colorvecsSto = SB::openColorvecStorage(params.named_obj.colorvec_files);
     
      //
      // DB storage
      //
      std::vector<BinaryStoreDB<SerialDBKey<KeyPropElementalOperator_t>,
				SerialDBData<ValPropElementalOperator_t>>>
	qdp_db{};
      SB::StorageTensor<8, SB::ComplexD> st;

      // Open the file, and write the meta-data and the binary for this operator
      if (!params.param.contract.use_superb_format)
      {
	qdp_db.resize(1);
	if (!qdp_db[0].fileExists(params.named_obj.prop_op_file))
	{
	  XMLBufferWriter file_xml;
	  push(file_xml, "DBMetaData");
	  write(file_xml, "id", std::string("propElemOp"));
	  write(file_xml, "lattSize", QDP::Layout::lattSize());
	  write(file_xml, "decay_dir", params.param.contract.decay_dir);
	  proginfo(file_xml); // Print out basic program info
	  write(file_xml, "Params", params.param);
	  write(file_xml, "Config_info", gauge_xml);
	  pop(file_xml);

	  std::string file_str(file_xml.str());
	  qdp_db[0].setMaxUserInfoLen(file_str.size());

	  qdp_db[0].open(params.named_obj.prop_op_file, O_RDWR | O_CREAT, 0664);

	  qdp_db[0].insertUserdata(file_str);
	}
	else
	{
	  qdp_db[0].open(params.named_obj.prop_op_file, O_RDWR, 0664);
	}
      }
      else
      {
	// Read order; letter meaning:
	//   n/N: source/sink eigenvector index
	//   s/q: source/sink spin index
	//   p/P: time slice source/sink
	//   h/H: phasing source/sink
	const char* order = "sqnNpPhH";
	XMLBufferWriter metadata_xml;
	push(metadata_xml, "DBMetaData");
	write(metadata_xml, "id", std::string("propElemOp"));
	write(metadata_xml, "lattSize", QDP::Layout::lattSize());
	write(metadata_xml, "decay_dir", params.param.contract.decay_dir);
	proginfo(metadata_xml); // Print out basic program info
	write(metadata_xml, "Config_info", gauge_xml);
	write(metadata_xml, "Params", params.param);
	write(metadata_xml, "mass_label", params.param.contract.mass_label);
	write(metadata_xml, "tensorOrder", order);
	std::vector<multi1d<int>> phasings_xml;
	for (const auto& it : phasings)
	  phasings_xml.push_back(SB::tomulti1d(it));
	write(metadata_xml, "phasings", phasings_xml);
	pop(metadata_xml);

	// NOTE: metadata_xml only has a valid value on Master node; so do a broadcast
	std::string metadata = SB::broadcast(metadata_xml.str());

	st = SB::StorageTensor<8, SB::ComplexD>(
	  params.named_obj.prop_op_file, metadata, order,
	  SB::kvcoors<8>(order,
			 {
			   {'s', Ns},
			   {'q', Ns},
			   {'n', params.param.contract.num_vecs},
			   {'N', params.param.contract.num_vecs},
			   {'p', Lt},
			   {'P', Lt},
			   {'h', (int)phasings.size()},
			   {'H', (int)phasings.size()},
			 }),
	  SB::Sparse, SB::checksum_type::BlockChecksum,
	  params.param.contract.output_file_is_local ? SB::LocalFSFile : SB::SharedFSFile);
      }

      QDPIO::cout << "Finished opening peram file" << std::endl;

      //
      // Try the factories
      //
      try
      {
	StopWatch swatch;
	swatch.reset();
	QDPIO::cout << "Try the various factories" << std::endl;

	//
	// Initialize fermion action and create the solver
	//
	SB::ChimeraSolver PP{params.param.prop.fermact, params.param.prop.invParam, u};

	swatch.start();

	//
	// Loop over the source color and spin, creating the source
	// and calling the relevant propagator routines.
	//
	const int num_vecs            = params.param.contract.num_vecs;
	const multi1d<int>& t_sources = params.param.contract.t_sources;
	const int max_rhs             = params.param.contract.max_rhs;

	// Set place for doing the contractions
	SB::DeviceHost dev =
	  params.param.contract.use_device_for_contractions ? SB::OnDefaultDevice : SB::OnHost;

	// Loop over each time source
	for (int tt = 0; tt < t_sources.size(); ++tt)
	{
	  int t_source = t_sources[tt]; // This is the actual time-slice.
	  QDPIO::cout << "t_source = " << t_source << std::endl;

	  // Compute the first tslice and the number of tslices involved in the contraction
	  int first_tslice = SB::normalize_coor(t_source - params.param.contract.Nt_backward, Lt);
	  int num_tslices = std::min(
	    params.param.contract.Nt_backward + std::max(1, params.param.contract.Nt_forward), Lt);

	  // Get `num_vecs` colorvecs, and `num_tslices` tslices starting from time-slice `first_tslice`
	  SB::Tensor<Nd + 3, SB::Complex> colorvec = SB::getColorvecs<SB::Complex>(
	    colorvecsSto, u, decay_dir, first_tslice, num_tslices, num_vecs, "cxyzXnt", SB::Coor<3>{{}}, dev);

	  for (const auto& it : phasing_pairs)
	  {
	    int phasing_src_idx = it.first;

	    // Get all eigenvectors for `t_source`
	    auto source_colorvec = SB::phaseColorvecs(
	      colorvec.kvslice_from_size({{'t', SB::normalize_coor(t_source - first_tslice, Lt)}},
					 {{'t', 1}}),
	      t_source, phasings[phasing_src_idx]);

	    for (int spin_source = 0; spin_source < Ns; ++spin_source)
	    {
	      // Invert the source for `spin_source` spin and retrieve `num_tslices` tslices starting from tslice `first_tslice`
	      // NOTE: s is spin source, and S is spin sink
	      SB::Tensor<Nd + 5, SB::Complex> quark_solns =
		SB::doInversion(PP, source_colorvec, t_source, first_tslice, num_tslices,
				{spin_source}, max_rhs, "cxyzXnSst");

	      StopWatch snarss1;
	      snarss1.reset();
	      snarss1.start();

	      for (int phasing_snk_idx : it.second)
	      {
		// Phase the sink
		auto colorvec_snk =
		  SB::phaseColorvecs(colorvec, first_tslice, phasings[phasing_snk_idx]);

		// Contract the distillation elements
		// NOTE: N: is colorvec in sink, and n is colorvec in source
		SB::Tensor<5, SB::Complex> elems(
		  "NnSst", {num_vecs, num_vecs, Ns, 1, num_tslices}, SB::OnHost,
		  !params.param.contract.use_superb_format ? SB::OnMaster : SB::OnEveryone);
		elems.contract(colorvec_snk, {{'n', 'N'}}, SB::Conjugate, quark_solns, {},
			       SB::NotConjugate);

		snarss1.stop();
		QDPIO::cout << "Time to contract for one spin source : "
			    << snarss1.getTimeInSeconds() << " secs" << std::endl;

		snarss1.reset();
		snarss1.start();

		if (!params.param.contract.use_superb_format)
		{
		  ValPropElementalOperator_t val;
		  val.mat.resize(num_vecs, num_vecs);
		  val.mat = zero;
		  auto local_elems = elems.getLocal();
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
		      if (local_elems)
		      {
			for (int colorvec_sink = 0; colorvec_sink < num_vecs; ++colorvec_sink)
			{
			  for (int colorvec_source = 0; colorvec_source < num_vecs;
			       ++colorvec_source)
			  {
			    std::complex<REAL> e = local_elems.get(
			      {colorvec_sink, colorvec_source, spin_sink, 0, i_tslice});
			    val.mat(colorvec_sink, colorvec_source).elem().elem().elem() =
			      RComplex<REAL64>(e.real(), e.imag());
			  }
			}
		      }
		      qdp_db[0].insert(key, val);
		    }
		  }
		}
		else
		{
		  st.kvslice_from_size({{'s', spin_source},
					{'p', t_source},
					{'P', first_tslice},
					{'h', phasing_src_idx},
					{'H', phasing_snk_idx}},
				       {{'s', 1}, {'p', 1}, {'P', num_tslices}, {'h', 1}, {'H', 1}})
		    .copyFrom(elems.rename_dims({{'S', 'q'}, {'t', 'P'}}));
		}
	      }

	      snarss1.stop();
	      QDPIO::cout << "Time to store the props : " << snarss1.getTimeInSeconds() << " secs"
			  << std::endl;
	    } // for spin_source
	  }   // for it
	}     // for tt

	swatch.stop();
	QDPIO::cout << "Propagators computed: time= " << swatch.getTimeInSeconds() << " secs"
		    << std::endl;
      } catch (const std::exception& e)
      {
	QDP_error_exit("%s: caught exception: %s\n", name.c_str(), e.what());
      }

      // Close colorvecs storage
      SB::closeColorvecStorage(colorvecsSto);

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
