/*! \file
 * \brief Inline measurement for compute correlation functions
 */

#include "chroma_config.h"

#ifdef BUILD_SB
// Activate the MPI support in Superbblas
#  define SUPERBBLAS_USE_MPI

// Activate redstar-datalib support for superbblas and with gpus
#  define USE_SUPERBBLAS
#  define USE_SUPERBBLAS_WITH_GPU_SUPPORT
#endif

#ifdef BUILD_REDSTAR_DATALIB
#include "algs/superb_contractions.h"
#include "hadron/process_corrs.h"
#include "io/adat_xmlio.h"
#include "io/baryon_superb.h"
#include "io/genprop4.h"
#include "io/meson_superb.h"
#include "io/prop_superb.h"

#include "handle.h"
#include "meas/glue/mesplq.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/hadron/inline_corr_superb_w.h"
#include "meas/inline/io/named_objmap.h"
#include "meas/inline/make_xml_file.h"
#include "meas/smear/link_smearing_aggregate.h"
#include "meas/smear/link_smearing_factory.h"
#include "util/ferm/mgproton.h"
#include "util/ferm/superb_contractions.h"
#include "util/info/proginfo.h"

namespace Chroma
{
  /*!
   * \ingroup hadron
   *
   * @{
   */
  namespace InlineCorrSuperbEnv
  {
    // Reader for input parameters
    void read(XMLReader& xml, const std::string& path,
	      InlineCorrSuperbEnv::Params::Param_t::FlavorToMass& param)
    {
      XMLReader paramtop(xml, path);

      read(paramtop, "flavor", param.flavor);
      read(paramtop, "mass", param.mass);
      const std::set<char> flavors{'c', 'e', 'l', 's', 'y', 'x'};
      if (flavors.count(param.flavor) == 0)
      {
	QDPIO::cerr << "invalid flavor: " << std::string{param.flavor} << std::endl;
	QDP_abort(1);
      }
    }

    // Writer for input parameters
    void write(XMLWriter& xml, const std::string& path,
	       const InlineCorrSuperbEnv::Params::Param_t::FlavorToMass& param)
    {
      push(xml, path);

      write(xml, "flavor", param.flavor);
      write(xml, "mass", param.mass);

      pop(xml);
    }

    // Reader for input parameters
    void read(XMLReader& xml, const std::string& path,
	      InlineCorrSuperbEnv::Params::Param_t::FlavorToProp& param)
    {
      XMLReader paramtop(xml, path);

      read(paramtop, "flavor", param.flavor);
      read(paramtop, "Propagator", param.prop);
      const std::set<char> flavors{'c', 'e', 'l', 's', 'y', 'x'};
      if (flavors.count(param.flavor) == 0)
      {
	QDPIO::cerr << "invalid flavor: " << std::string{param.flavor} << std::endl;
	QDP_abort(1);
      }
    }

    // Writer for input parameters
    void write(XMLWriter& xml, const std::string& path,
	       const InlineCorrSuperbEnv::Params::Param_t::FlavorToProp& param)
    {
      push(xml, path);

      write(xml, "flavor", param.flavor);
      write(xml, "Propagator", param.prop);

      pop(xml);
    }

    // Reader for input parameters
    void read(XMLReader& xml, const std::string& path,
	      InlineCorrSuperbEnv::Params::Param_t& param)
    {
      XMLReader paramtop(xml, path);

      param.num_vecs = 0;
      read(paramtop, "num_vecs", param.num_vecs);

      param.max_rhs = 0;
      if (paramtop.count("max_rhs") > 0)
      {
	read(paramtop, "max_rhs", param.max_rhs);
      }

      if (paramtop.count("flavor_to_mass") > 0)
      {
	read(paramtop, "flavor_to_mass", param.flavor_to_mass);
      }

      if (paramtop.count("flavor_to_prop") > 0)
      {
	read(paramtop, "flavor_to_prop", param.flavor_to_prop);
      }

      param.t_origin = 0;
      if (paramtop.count("t_origin") > 0)
      {
	read(paramtop, "t_origin", param.t_origin);
      }

      if (paramtop.count("mesons") > 0)
      {
	read(paramtop, "mesons", param.meson_files);
      }

      if (paramtop.count("baryons") > 0)
      {
	read(paramtop, "baryons", param.baryon_files);
      }

      if (paramtop.count("props") > 0)
      {
	read(paramtop, "props", param.prop_files);
      }

      if (paramtop.count("genprops") > 0)
      {
	read(paramtop, "genprops", param.genprop_files);
      }

      param.decay_dir = 3;

      param.Nt_forward = 0;
      if (paramtop.count("Nt_forward") > 0)
      {
	read(paramtop, "Nt_forward", param.Nt_forward);
      }

      if (paramtop.count("ensemble") > 0)
      {
	read(paramtop, "ensemble", param.ensemble);
      }

      param.testing = false;
      if (paramtop.count("testing") > 0)
      {
	read(paramtop, "testing", param.testing);
      }

      param.link_smearing = readXMLGroup(paramtop, "LinkSmearing", "LinkSmearingType");
    }

    // Writer for input parameters
    void write(XMLWriter& xml, const std::string& path,
	       const InlineCorrSuperbEnv::Params::Param_t& param)
    {
      push(xml, path);

      write(xml, "num_vecs", param.num_vecs);
      write(xml, "max_rhs", param.max_rhs);
      write(xml, "flavor_to_mass", param.flavor_to_mass);
      write(xml, "flavor_to_prop", param.flavor_to_prop);
      write(xml, "t_origin", param.t_origin);
      write(xml, "mesons", param.meson_files);
      write(xml, "baryons", param.baryon_files);
      write(xml, "props", param.prop_files);
      write(xml, "genprops", param.genprop_files);
      write(xml, "Nt_forward", param.Nt_forward);
      write(xml, "ensemble", param.ensemble);
      write(xml, "testing", param.testing);
      xml << param.link_smearing.xml;

      pop(xml);
    }

    //! Read named objects
    void read(XMLReader& xml, const std::string& path,
	      InlineCorrSuperbEnv::Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_id", input.gauge_id);
      read(inputtop, "colorvec_files", input.colorvec_files);
      read(inputtop, "corr_file", input.corr_file);
      read(inputtop, "corr_graph_file", input.corr_graph_file);
    }

    //! Write named objects
    void write(XMLWriter& xml, const std::string& path,
	       const InlineCorrSuperbEnv::Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "gauge_id", input.gauge_id);
      write(xml, "colorvec_files", input.colorvec_files);
      write(xml, "corr_file", input.corr_file);
      write(xml, "corr_graph_file", input.corr_graph_file);

      pop(xml);
    }

    // Writer for input parameters
    void write(XMLWriter& xml, const std::string& path,
	       const InlineCorrSuperbEnv::Params& param)
    {
      param.writeXML(xml, path);
    }
  }

  namespace InlineCorrSuperbEnv
  {
    // Anonymous namespace for registration
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, const std::string& path)
      {
	return new InlineMeas(Params(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "CORR_SUPERB";

    //! Register all the factories
    bool registerAll()
    {
      bool success = true;
      if (!registered)
      {
	success &= LinkSmearingEnv::registerAll();
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
    }

    //----------------------------------------------------------------------------
    // Param stuff
    Params::Params()
    {
      frequency = 0;
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
      } catch (const std::string& e)
      {
	QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << std::endl;
	QDP_abort(1);
      }
    }

    void Params::writeXML(XMLWriter& xml_out, const std::string& path) const
    {
      push(xml_out, path);

      // Parameters for source construction
      write(xml_out, "Param", param);

      // Write out the output propagator/source configuration info
      write(xml_out, "NamedObject", named_obj);

      pop(xml_out);
    }

    //-------------------------------------------------------------------------------
    // Function call
    void InlineMeas::operator()(unsigned long update_no, XMLWriter& xml_out)
    {
      // If xml file not empty, then use alternate
      if (params.xml_file != "")
      {
	std::string xml_file = makeXMLFileName(params.xml_file, update_no);

	push(xml_out, "CorrSuperbVec");
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

#  if defined(USE_SUPERBBLAS) && !defined(SUPERBNOVA_DEBUG)
    template <std::size_t N>
    SBN::superbblas_implementation::detail::Distribution
    get_sbn_distribution(const std::string& order, const SB::Distribution& dist,
			 const SB::detail::TensorPartition<N>& p)
    {
      const auto kind = p.isLocal
			  ? SBN::superbblas_implementation::detail::Distribution::Local
			  : SBN::superbblas_implementation::detail::Distribution::Distributed;
      const auto dim = SBN::Coor(p.dim.begin(), p.dim.end());
      std::vector<unsigned char> distributed_directions;
      if (kind != SBN::superbblas_implementation::detail::Distribution::Local)
      {
	if (dist != SB::OnMaster && !(dist.size() > 2 && dist.at(0) == '_' && dist.at(1) == '_'))
	{
	  for (char d : dist)
	  {
	    auto s = std::find(order.begin(), order.end(), d);
	    if (s != order.end())
	      distributed_directions.push_back(s - order.begin());
	  }
	}
      }
      SBN::superbblas_implementation::detail::Distribution::Fs fs(
	{(int)dim.size(), 2, (int)p.p.size()});
      for (int proc = 0; proc < (int)p.p.size(); ++proc)
      {
	for (int i = 0; i < 2; ++i)
	{
	  std::copy(p.p.at(proc).at(i).begin(), p.p.at(proc).at(i).end(), fs.begin({0, i, proc}));
	}
      }
      return SBN::superbblas_implementation::detail::Distribution(kind, dim, distributed_directions,
								  fs);
    }
#  endif // defined(USE_SUPERBBLAS) && !defined(SUPERBNOVA_DEBUG)

    /// Return a view of the tensor object
    /// \param t: tensor

    template <std::size_t N>
    SBN::Tensor toTensor(const SB::Tensor<N, SB::ComplexD>& t, bool allow_copy = true)
    {
      // Check that the input tensor is not fake complex or example gratia
      if (t.complexLabel != 0 || t.eg || t.dist == SB::Glocal)
	throw std::runtime_error("toTensor: unsupported tensor");

      const SBN::Coor from(t.from.begin(), t.from.end());
      const SBN::Coor size(t.size.begin(), t.size.end());
      const SBN::Coor dim(t.dim.begin(), t.dim.end());
      static_assert(std::is_same<SB::ComplexD, SBN::value_type>::value,
		    "superbnova has an invalid type");
#  if defined(USE_SUPERBBLAS) && !defined(SUPERBNOVA_DEBUG)
      auto op_ptr = t.data();
      const auto dev = t.getDev();
      const auto ctx =
	(dev == SB::OnHost ? SBN::superbblas_implementation::detail::getCpuContext()
			   : SBN::superbblas_implementation::detail::getGpuContext());
      auto alloc =
	std::make_shared<SBN::superbblas_implementation::detail::Allocation<SBN::value_type>>(
	  op_ptr, ctx);
      auto p = std::make_shared<SBN::superbblas_implementation::detail::Distribution>(
	get_sbn_distribution(t.order, t.dist, *t.p));
      return SBN::Tensor{from, size, dim, t.order, alloc, p, t.scalar, t.conjugate, 0, 0, {}};
#  else
      if (!allow_copy)
	throw std::runtime_error("unsupported");
      const auto trep = t.make_sure(SB::none, SB::OnHost, SB::OnEveryoneReplicated);
      auto p = (const std::complex<double>*)trep.data();
      auto alloc = std::make_shared<std::vector<std::complex<double>>>(p, p + SBN::volume(dim));
      return SBN::Tensor{from, size, dim, trep.order, alloc, trep.scalar, trep.conjugate, 0, 0, {}};
#  endif
    }

    /// Return the mesons without spins
    /// \param db: meson storage
    /// \param colorvecsSto: colorvec storage
    /// \param u: original gauge field
    /// \param u_smr: smeared original gauge field
    /// \param meson_keys: list of mesons keys
    /// \param perms: list of permutations, one for each meson key
    /// \param do_conj: list of whether to conjugate the meson, one for each meson key
    /// \param ev_from: first eigenvector to return
    /// \param ev_size: number of eigenvectors to return
    /// \param dist_labels: dimensions to be distributed, some of "vwi"
    /// \param alloc: allocation for the returned tensor

    inline SBN::Tensor
    get_meson_elementals(const ADATIO::StorageMeson& db, const SB::ColorvecsStorage& colorvecsSto,
			 const multi1d<LatticeColorMatrix>& u,
			 const multi1d<LatticeColorMatrix>& u_smr,
			 const std::vector<Hadron::KeyMesonElementalOperator_t>& meson_keys,
			 const std::vector<SBN::Coor>& perms, const std::vector<bool>& do_conj,
			 const SBN::Coor& ev_from, const SBN::Coor& ev_size,
			 const std::string& dist_labels, const SBN::Tensor& guide, bool testing)
    {
      // Check input
      if (meson_keys.size() != perms.size())
	throw std::runtime_error("invalid input");
      if (ev_from.size() != 2 || ev_size.size() != 2)
	throw std::runtime_error("invalid input");
      if (superbblas::getDebugLevel() > 0)
      {
	for (const auto& it : perms)
	{
	  if (it.size() != 2)
	    throw std::runtime_error("invalid input");
	  for (const auto& i : it)
	    if (i < 0 || i >= 2)
	      throw std::runtime_error("invalid input");
	}
      }

      const int num_vecs = std::max(ev_from.at(0) + ev_size.at(0), ev_from.at(1) + ev_size.at(1));

      // Create output tensor
      SBN::Tensor mesons = SBN::create_tensor_with_local_components(
	{ev_size.at(0), ev_size.at(1), -(int)meson_keys.size()}, "vwi", dist_labels,
	SBN::Options::Alloc::Device, 0, 0, SBN::Options::IsEg::False, guide);
      const auto first_local_meson = SBN::get_local_srange(mesons).at(0).at('i');

      // Try to get the mesons from the storage and annotate the missing keys
      std::vector<std::tuple<Hadron::KeyMesonElementalOperator_t, SBN::Coor, int>>
	local_missing_mesons;
      {
	local_missing_mesons.reserve(meson_keys.size());
	const auto& local_mesons =
	  SBN::slice_kv(SBN::get_local_tensor(mesons), {{'i', first_local_meson}},
			{{'i', (int)meson_keys.size()}});
	Hadron::ValMesonElementalOperator_t val;
	for (std::size_t i = 0; i < meson_keys.size(); ++i)
	{
	  const auto& key = meson_keys.at(i);
	  const auto& p = perms.at(i);
	  if (do_conj.at(i))
	    throw std::runtime_error("unsupported case");
	  if (db.get(key, val) != 0)
	  {
	    // We can't miss keys when testing
	    if (testing)
	    {
	      throw std::runtime_error(std::string("doing testing and missing meson: ") +
				       ::SB::getXML(key));
	    }

	    // Record missing meson
	    local_missing_mesons.push_back({key, p, first_local_meson + i});
	  }
	  else
	  {
	    if (val.op.size1() < num_vecs || val.op.size2() < num_vecs)
	      throw std::runtime_error("got a meson with insufficient number of vectors");
	    auto ti = SBN::toTensor(val.op, "vw", SBN::Options::Distribution::Local);
	    ti = Hadron::detail::apply_vertex_perm(ti, p, "vw");
	    ti = SBN::slice_kv(ti, //
			       {{'v', ev_from.at(0)}, {'w', ev_from.at(0)}},
			       {{'v', ev_size.at(1)}, {'w', ev_size.at(1)}});
	    SBN::copyTo(ti, SBN::slice_kv(local_mesons, {{'i', i}}, {{'i', 1}}));

	    // If testing, also record it as a missing meson
	    if (testing)
	      local_missing_mesons.push_back({key, p, first_local_meson + i});
	  }
	}
      }

      // If testing, save the output tensor and set it to zero
      SBN::Tensor true_mesons;
      if (testing)
      {
	true_mesons = mesons;
	mesons = SBN::like_this(mesons);
      }

      // Recompile for each missing time slice, the phases, momenta, and displacement to compute
      using displacement_t = std::vector<int>;
      using tslice_left_right_phases = std::tuple<int, SB::Coor<3>, SB::Coor<3>>;
      using momenta_displacements_indices = std::tuple<
	Hadron::detail::unordered_map<SB::Coor<3>, int>,
	Hadron::detail::unordered_map<displacement_t, int>,
	Hadron::detail::unordered_multimap<std::array<int, 2>, std::tuple<SBN::Coor, int>>>;
      Hadron::detail::unordered_map<tslice_left_right_phases, momenta_displacements_indices>
	from_tslice_left_right_phases_to_momenta_displacement;
      const auto get_index = [=](auto& map, const auto& value) {
	const auto s = map.size();
	if (map.count(value) == 0)
	  map[value] = s;
	return map.at(value);
      };
      const auto get_vector = [=](const auto& map) {
	std::vector<typename std::remove_reference<decltype(map)>::type::key_type> r(map.size());
	for (const auto& it : map)
	  r.at(it.second) = it.first;
	return r;
      };
      for (const auto& missing_mesons_in_some_process : SBN::gather(local_missing_mesons))
      {
	for (const auto& [meson_key, perm, index] : missing_mesons_in_some_process)
	{
	  const auto& k = tslice_left_right_phases{
	    meson_key.t_slice, ADATIO::detail::toCoor(meson_key.phasing_sink),
	    ADATIO::detail::toCoor(meson_key.phasing_source)};
	  auto& v = from_tslice_left_right_phases_to_momenta_displacement[k];
	  auto mom_index = get_index(std::get<0>(v), ADATIO::detail::toCoor(meson_key.mom));
	  auto disp_index = get_index(std::get<1>(v), meson_key.displacement);
	  std::get<2>(v).insert({{mom_index, disp_index}, {perm, index}});
	}
      }

      for (const auto& it : from_tslice_left_right_phases_to_momenta_displacement)
      {
	const auto& [t_source, left_phase, right_phase] = it.first;
	const auto& moms = get_vector(std::get<0>(it.second));
	const auto& disps = get_vector(std::get<1>(it.second));
	const auto& mom_disp_to_index = std::get<2>(it.second);

	// Get num_vecs colorvecs on time-slice t_source
	const int decay_dir = 3;
	SB::Tensor<Nd + 3, SB::Complex> source_colorvec = SB::getColorvecs<SB::Complex>(
	  colorvecsSto, u, decay_dir, t_source, 1, num_vecs, SB::none);

	// Callback
	const auto call = [&](SB::Tensor<4, SB::ComplexD> tensor, int disp, int first_tslice,
			      int first_mom) {
	  if (tensor.kvdim().at('m') != 1 || tensor.kvdim().at('t') != 1)
	    throw std::runtime_error("wtf");
	  auto range = mom_disp_to_index.equal_range({first_mom, disp});
	  for (auto it = range.first; it != range.second; ++it)
	  {
	    const auto& [p, index] = it->second;
	    auto ti = SBN::relabel(toTensor(tensor), {{'i', 'w'}, {'j', 'v'}});
	    ti = Hadron::detail::apply_vertex_perm(ti, p, "vw");
	    ti = SBN::slice_kv(ti, //
			       {{'v', ev_from.at(0)}, {'w', ev_from.at(1)}},
			       {{'v', ev_size.at(0)}, {'w', ev_size.at(1)}});
	    SBN::copyTo(ti, SBN::slice_kv(mesons, {{'i', index}}, {{'i', 1}}));
	  }
	};

	// Do the contractions
	const bool use_derivP = true;
	SB::doMomDisp_contractions<Nd + 3, SB::ComplexD>(
	  u_smr, source_colorvec.make_sure<SB::ComplexD>(), left_phase, right_phase, moms, t_source,
	  disps, use_derivP, call, SB::none, SB::OnDefaultDevice, SB::OnEveryone,
	  0 /* max_tslices_in_contraction==0 means to do all */, 1 /* max_moms_in_contraction */);
      }

      // If testing, make sure that recomputed mesons are similar to the ones got from storage
      if (testing)
      {
	for (int i = 0, n = SBN::get_kv_size(mesons).at('i'); i < n; ++i)
	{
	  const auto e_i = SBN::slice_kv(mesons, {{'i', i}}, {{'i', 1}});
	  const auto true_e_i = SBN::slice_kv(true_mesons, {{'i', i}}, {{'i', 1}});
	  auto diff = SBN::clone(e_i);
	  SBN::addTo(SBN::scale(true_e_i, -1), diff);
	  if (SBN::frob(true_e_i) * 1e-5 < SBN::frob(diff))
	    throw std::runtime_error("mesons not passing test");
	}
      }

      return mesons;
    }

    inline std::tuple<Hadron::KeyBaryonElementalOperator_t, int>
    get_perm_baryon_key(const Hadron::KeyBaryonElementalOperator_t& key, const SB::Coor<3>& perm)
    {
      // Get the new key
      auto perm_key = key;
      const auto& disps = std::array<std::vector<int>, 3>{key.left, key.middle, key.right};
      auto perm_disps =
	std::array<std::vector<int>, 3>{std::vector<int>{}, std::vector<int>{}, std::vector<int>{}};
      for (std::size_t i = 0; i < 3; ++i)
	perm_disps.at(perm.at(i)) = disps.at(i);
      perm_key.left = std::move(perm_disps.at(0));
      perm_key.middle = std::move(perm_disps.at(1));
      perm_key.right = std::move(perm_disps.at(2));

      // List all permutation of three elements; and the sign is the evenness of the number of single exchanges
      const std::vector<std::pair<std::array<int, 3>, int>> perm_and_sign_list{
	{SB::Coor<3>{0, 1, 2}, 1},  //
	{SB::Coor<3>{0, 2, 1}, -1}, //
	{SB::Coor<3>{1, 0, 2}, -1}, //
	{SB::Coor<3>{1, 2, 0}, 1},  //
	{SB::Coor<3>{2, 0, 1}, 1},  //
	{SB::Coor<3>{2, 1, 0}, -1}  //
      };

      // Find scalar
      int scalar = 0;
      for (const auto& it : perm_and_sign_list)
      {
	if (it.first == perm)
	{
	  scalar = it.second;
	  break;
	}
      }
      if (scalar == 0)
	throw std::runtime_error("wtf");

      return {perm_key, scalar};
    }

    /// Return the baryons without spins
    /// \param db: baryon storage
    /// \param colorvecsSto: colorvec storage
    /// \param u: original gauge field
    /// \param u_smr: smeared original gauge field
    /// \param baryon_keys: list of baryons keys
    /// \param perms: list of permutations, one for each baryon key
    /// \param do_conj: list of whether to conjugate the baryon, one for each baryon key
    /// \param ev_from: first eigenvector to return for "vwx"
    /// \param ev_size: number of eigenvectors to return for "vwx"
    /// \param dist_labels: dimensions to be distributed, some of "vwxi"
    /// \param alloc: allocation for the returned tensor

    inline SBN::Tensor
    get_baryon_elementals(const ADATIO::StorageBaryon& db, const SB::ColorvecsStorage& colorvecsSto,
			  const multi1d<LatticeColorMatrix>& u,
			  const multi1d<LatticeColorMatrix>& u_smr,
			  const std::vector<Hadron::KeyBaryonElementalOperator_t>& baryon_keys,
			  const std::vector<SBN::Coor>& perms, const std::vector<bool>& do_conj,
			  const SBN::Coor& ev_from, const SBN::Coor& ev_size,
			  const std::string& dist_labels, const SBN::Tensor& guide, bool testing)
    {
      if (baryon_keys.size() != perms.size() || baryon_keys.size() != do_conj.size())
	throw std::runtime_error("invalid input");
      if (ev_from.size() != 3 || ev_size.size() != 3)
	throw std::runtime_error("invalid input");
      if (SBN::detail::get_debug_level() > 0)
      {
	for (const auto& it : perms)
	{
	  if (it.size() != 3)
	    throw std::runtime_error("invalid input");
	  for (const auto& i : it)
	    if (i < 0 || i >= 3)
	      throw std::runtime_error("invalid input");
	}
      }

      const auto toCoor3 = [=](const SBN::Coor& v) {
	return SB::Coor<3>{v.at(0), v.at(1), v.at(2)};
      };

      const int num_vecs =
	std::max(std::max(ev_from.at(0) + ev_size.at(0), ev_from.at(1) + ev_size.at(1)),
		 ev_from.at(2) + ev_size.at(2));

      // Create output tensor
      SBN::Tensor baryons = SBN::create_tensor_with_local_components(
	SBN::concat(ev_size, {-(int)baryon_keys.size()}), "vwxi", dist_labels,
	SBN::Options::Alloc::Device, 0, 0, SBN::Options::IsEg::False, guide);

      // Try to get the mesons from the storage and annotate the missing keys
      std::vector<std::tuple<Hadron::KeyBaryonElementalOperator_t, int, bool, int>>
	local_missing_baryons;
      {
	local_missing_baryons.reserve(baryon_keys.size());
	const auto first_local_baryon = SBN::get_local_srange(baryons).at(0).at('i');
	const auto& local_baryons =
	  SBN::slice_kv(SBN::get_local_tensor(baryons), {{'i', first_local_baryon}},
			{{'i', (int)baryon_keys.size()}});
	ADATIO::SerialDBData<Hadron::ValBaryonElementalOperator_t> val;
	for (std::size_t i = 0; i < baryon_keys.size(); ++i)
	{
	  const auto& baryon_i = SBN::slice_kv(local_baryons, {{'i', i}}, {{'i', 1}});
	  if (!has_local_support(baryon_i))
	    continue;
	  const auto& key = baryon_keys.at(i);
	  if (db.get(key, val) != 0)
	  {
	    // We can't miss keys when testing
	    if (testing)
	    {
	      throw std::runtime_error(std::string("doing testing and missing baryon: ") +
				       ::SB::getXML(key));
	    }

	    // Record missing baryon
	    const auto& [norm_key, scalar] = get_perm_baryon_key(key, toCoor3(perms.at(i)));
	    local_missing_baryons.push_back(
	      {norm_key, scalar, do_conj.at(i), first_local_baryon + i});
	  }
	  else
	  {
	    if (val.data().op.size1() < num_vecs || val.data().op.size2() < num_vecs ||
		val.data().op.size3() < num_vecs)
	      throw std::runtime_error("got a baryon with insufficient number of vectors");
	    auto ti = SBN::toTensor(val.data().op, "vwx", SBN::Options::Distribution::Local);
	    ti = Hadron::detail::apply_vertex_perm(ti, perms.at(i), "vwx");
	    ti = SBN::slice_kv(ti, SBN::get_scoor("vwx", ev_from), SBN::get_scoor("vwx", ev_size));
	    SBN::copyTo(do_conj.at(i) ? SBN::conj(ti) : ti, baryon_i);

	    // If testing, also record it as a missing baryon
	    if (testing)
	    {
	      const auto& [norm_key, scalar] = get_perm_baryon_key(key, toCoor3(perms.at(i)));
	      local_missing_baryons.push_back(
		{norm_key, scalar, do_conj.at(i), first_local_baryon + i});
	    }
	  }
	}
      }

      // If testing, save the output tensor and set it to zero
      SBN::Tensor true_baryons;
      if (testing)
      {
	true_baryons = baryons;
	baryons = SBN::like_this(baryons);
      }

      // Recompile for each missing time slice, the phases, momenta, and displacement to compute
      using displacement_t = std::array<std::vector<int>, 3>;
      using tslice_phase = std::tuple<int, SB::Coor<3>>;
      using momenta_displacements_indices = std::tuple<
	Hadron::detail::unordered_map<SB::Coor<3>, int>,
	Hadron::detail::unordered_map<displacement_t, int>,
	Hadron::detail::unordered_multimap<std::array<int, 2>, std::tuple<int, bool, int>>>;
      Hadron::detail::unordered_map<tslice_phase, momenta_displacements_indices>
	from_tslice_and_phase_to_momenta_displacement;
      const auto get_index = [=](auto& map, const auto& value) {
	const auto s = map.size();
	if (map.count(value) == 0)
	  map[value] = s;
	return map.at(value);
      };
      const auto get_vector = [=](const auto& map) {
	std::vector<typename std::remove_reference<decltype(map)>::type::key_type> r(map.size());
	for (const auto& it : map)
	  r.at(it.second) = it.first;
	return r;
      };
      for (const auto& missing_baryons_in_some_process : SBN::gather(local_missing_baryons))
      {
	for (const auto& [baryon_key, scalar, do_conj, index] : missing_baryons_in_some_process)
	{
	  const auto& k =
	    tslice_phase{baryon_key.t_slice, ADATIO::detail::toCoor(baryon_key.phasing)};
	  auto& v = from_tslice_and_phase_to_momenta_displacement[k];
	  auto mom_index = get_index(std::get<0>(v), ADATIO::detail::toCoor(baryon_key.mom));
	  auto disp_index = get_index(
	    std::get<1>(v), displacement_t{baryon_key.left, baryon_key.middle, baryon_key.right});
	  std::get<2>(v).insert({{mom_index, disp_index}, {scalar, do_conj, index}});
	}
      }

      for (const auto& it : from_tslice_and_phase_to_momenta_displacement)
      {
	const int t_slice = std::get<0>(it.first);
	const auto& phase = std::get<1>(it.first);
	const auto& moms = get_vector(std::get<0>(it.second));
	const auto& disps = get_vector(std::get<1>(it.second));
	const auto& mom_disp_to_index = std::get<2>(it.second);

	// Get num_vecs colorvecs on time-slice t_slice
	const int decay_dir = 3;
	SB::Tensor<Nd + 3, SB::Complex> source_colorvec =
	  SB::getColorvecs<SB::Complex>(colorvecsSto, u, decay_dir, t_slice, 1, num_vecs, SB::none);
	source_colorvec = SB::phaseColorvecs(source_colorvec, t_slice, phase);

	// Callback
	const auto call = SB::ColorContractionFn<SB::Complex>(
	  [&](SB::Tensor<5, SB::ComplexD> tensor, int disp, int first_tslice, int first_mom) {
	    auto range = mom_disp_to_index.equal_range({first_mom, disp});
	    for (auto it = range.first; it != range.second; ++it)
	    {
	      const auto& [scalar, do_conj, index] = it->second;
	      auto ti = SBN::relabel(toTensor(tensor), {{'i', 'v'}, {'j', 'w'}, {'k', 'x'}});
	      ti =
		SBN::slice_kv(ti, //
			      {{'v', ev_from.at(0)}, {'w', ev_from.at(1)}, {'x', ev_from.at(2)}},
			      {{'v', ev_size.at(0)}, {'w', ev_size.at(1)}, {'x', ev_size.at(2)}});
	      if (do_conj)
		ti = SBN::conj(ti);
	      SBN::copyTo(SBN::scale(ti, (double)scalar),
			  SBN::slice_kv(baryons, {{'i', index}}, {{'i', 1}}));
	    }
	  });

	// Do the contractions
	const bool use_derivP = true;
	const int max_moms_in_contraction = 1;
	const int max_vecs = 4;
	SB::doMomDisp_colorContractions(u_smr, source_colorvec.make_sure<SB::ComplexD>(), moms,
					t_slice, disps, use_derivP, call,
					0 /* it means to do all */, max_moms_in_contraction,
					max_vecs, SB::none, SB::OnDefaultDevice, SB::OnEveryone);
      }

      // If testing, make sure that recomputed baryons are similar to the ones got from storage
      if (testing)
      {
	for (int i = 0, n = SBN::get_kv_size(baryons).at('i'); i < n; ++i)
	{
	  const auto e_i = SBN::slice_kv(baryons, {{'i', i}}, {{'i', 1}});
	  const auto true_e_i = SBN::slice_kv(true_baryons, {{'i', i}}, {{'i', 1}});
	  auto diff = SBN::clone(e_i);
	  SBN::addTo(SBN::scale(true_e_i, -1), diff);
	  if (SBN::frob(true_e_i) * 1e-5 < SBN::frob(diff))
	    throw std::runtime_error("baryons not passing test");
	}
      }

      return baryons;
    }

    /// Return the props
    /// \param db: prop storage
    /// \param colorvecsSto: colorvec storage
    /// \param u: original gauge field
    /// \param get_prop: get solver from mass label
    /// \param prop_keys: list of props keys
    /// \param max_rhs: maximum RHS to solve at once
    /// \param perms: list of permutations, one for each prop key
    /// \param do_conj: list of whether to conjugate the prop, one for each prop key
    /// \param ev_from: first eigenvector to return for "vwx"
    /// \param ev_size: number of eigenvectors to return for "vwx"
    /// \param dist_labels: dimensions to be distributed, some of "vwxi"
    /// \param alloc: allocation for the returned tensor

    inline SBN::Tensor
    get_prop_elementals(const ADATIO::StorageProp4& db, const SB::ColorvecsStorage& colorvecsSto,
			const multi1d<LatticeColorMatrix>& u,
			const std::function<SB::ChimeraSolver(std::string)>& get_prop, int max_rhs,
			const std::vector<Hadron::KeyProp4ElementalOperator_t>& prop_keys,
			const std::vector<SBN::Coor>& perms, const std::vector<bool>& do_conj,
			const SBN::Coor& ev_from, const SBN::Coor& ev_size,
			const std::string& dist_labels, const SBN::Tensor& guide, bool testing)
    {
      if (prop_keys.size() != perms.size() || prop_keys.size() != do_conj.size())
	throw std::runtime_error("invalid input");
      if (ev_from.size() != 4 || ev_size.size() != 4)
	throw std::runtime_error("invalid input");
      if (SBN::detail::get_debug_level() > 0)
      {
	for (const auto& it : perms)
	{
	  if (it != SBN::Coor{0, 1, 2, 3} && it != SBN::Coor{1, 0, 3, 2})
	    throw std::runtime_error("invalid input");
	}
      }

      const int num_vecs = std::max(ev_from.at(0) + ev_size.at(0), ev_from.at(1) + ev_size.at(1));
      if (max_rhs == 0)
	max_rhs = num_vecs;

      // This object is in the DR basis; create matrices to convert it to DP
      const auto& dr_left_global =
	SBN::slice_kv(Hadron::detail::adjForSpins(
			Hadron::detail::getDiracToDRMat(SBN::Options::Distribution::Replicated)), //
		      {{'s', ev_from.at(2)}}, {{'s', ev_size.at(2)}});
      const auto& dr_left = SBN::get_local_tensor(dr_left_global);
      const auto& dr_right = SBN::slice_kv(Hadron::detail::getDiracToDRMat(),
					   {{'S', ev_from.at(3)}}, {{'S', ev_size.at(3)}});

      // Create matrices to convert it to DP and with pre/post applying \gamma_5
      const auto dr_g5_left_global = Hadron::detail::contractSpins(
	dr_left_global, Hadron::detail::chromaGamma5(SBN::Options::Distribution::Replicated));
      const auto& dr_g5_left = SBN::get_local_tensor(dr_g5_left_global);
      const auto dr_g5_right =
	Hadron::detail::contractSpins(Hadron::detail::chromaGamma5(), dr_right);

      // Create chroma versions of dr_right and dr_g5_right
      auto dr_right_chroma = SB::Tensor<2, SB::Complex>("Ss", {{Ns, ev_size.at(3)}}, SB::OnHost,
							SB::OnEveryoneReplicated);
      SBN::copyTo(SBN::relabel(dr_right, {{'S', 's'}, {'s', 'S'}}),
		  SBN::get_local_tensor(toTensor(dr_right_chroma, false /* don't copy */)));
      auto dr_g5_right_chroma = SB::Tensor<2, SB::Complex>("Ss", {{Ns, ev_size.at(3)}}, SB::OnHost,
							   SB::OnEveryoneReplicated);
      SBN::copyTo(SBN::relabel(dr_g5_right, {{'S', 's'}, {'s', 'S'}}),
		  SBN::get_local_tensor(toTensor(dr_g5_right_chroma, false /* don't copy */)));

      // Create output tensor
      SBN::Tensor props = SBN::create_tensor_with_local_components(
	SBN::concat(ev_size, {-(int)prop_keys.size()}), "vwrsi", dist_labels,
	SBN::Options::Alloc::Host, 0, 0, SBN::Options::IsEg::False, guide);

      // Try to get the mesons from the storage and annotate the missing keys
      std::vector<std::tuple<Hadron::KeyProp4ElementalOperator_t, bool, int>> local_missing_props;
      {
	local_missing_props.reserve(prop_keys.size());
	const auto first_local_prop = SBN::get_local_srange(props).at(0).at('i');
	const auto& local_props = //
	  SBN::relabel(		  //
	    SBN::slice_kv(	  //
	      SBN::get_local_tensor(props), {{'i', first_local_prop}},
	      {{'i', (int)prop_keys.size()}}),
	    {{'r', 's'}, {'s', 'S'}});
	Hadron::ValProp4ElementalOperator_t val;
	const auto record_missing_key = [&](Hadron::KeyProp4ElementalOperator_t key, bool this_conj,
					    bool is_swap, int local_index) {
	  // Don't record requiring doing sink-source swapping, just whether to apply \gamma_5
	  if (is_swap)
	  {
	    std::swap(key.t_slice, key.t_source);
	    key.phasing_sink = -key.phasing_sink;
	    std::swap(key.phasing_source, key.phasing_sink);
	    key.phasing_sink = -key.phasing_sink;
	    is_swap = !is_swap;
	    this_conj = !this_conj;
	  }

	  // Record missing prop
	  local_missing_props.push_back({key, this_conj, first_local_prop + local_index});
	};
	for (std::size_t i = 0; i < prop_keys.size(); ++i)
	{
	  const auto& prop_i = SBN::slice_kv(local_props, {{'i', i}}, {{'i', 1}});
	  if (!has_local_support(prop_i))
	    continue;
	  auto key = prop_keys.at(i);
	  bool is_swap = (perms.at(i) != SBN::Coor{0, 1, 2, 3});
	  bool this_conj = do_conj.at(i);
	  bool gotit = false;
	  for (int attempt = 0; attempt < 2; ++attempt)
	  {
	    if (db.get(key, val) == 0)
	    {
	      gotit = true;
	      break;
	    }
	    if (attempt == 0)
	    {
	      // If it fails, try reversing sink and source
	      std::swap(key.t_slice, key.t_source);
	      key.phasing_sink = -key.phasing_sink;
	      std::swap(key.phasing_source, key.phasing_sink);
	      key.phasing_sink = -key.phasing_sink;
	      is_swap = !is_swap;
	      this_conj = !this_conj;
	    }
	    else
	    {
	      // We can't miss keys when testing
	      if (testing)
	      {
		throw std::runtime_error(std::string("doing testing and missing prop: ") +
					 ::SB::getXML(key));
	      }

	      // Record missing prop
	      record_missing_key(key, this_conj, is_swap, i);
	      break;
	    }
	  }
	  if (gotit)
	  {
	    if (val.op.size3() < num_vecs || val.op.size4() < num_vecs)
	    {
	      throw std::runtime_error("got a propagator with insufficient number of vectors");
	    }
	    // If sink and source are swapped, then conjugate and apply \gamma_5 left and right:
	    // V_t0' D^{-1} V_t1 =
	    //          [(V_t0' D^{-1} V_t1)']' =
	    //          [V_t1' D^{-\dagger} V_t0]' =
	    //          [\g_5 V_t1' D^{-1} V_t0 \g_5]' =
	    //          \g_5 [V_t1' D^{-1} V_t0]' \g_5
	    const auto& ti0 =
	      SBN::toTensor(val.op, is_swap ? "wvSs" : "vwsS", SBN::Options::Distribution::Local);
	    auto ti = SBN::slice_kv(is_swap ? SBN::conj(ti0) : ti0, //
				    {{'v', ev_from.at(0)}, {'w', ev_from.at(1)}},
				    {{'v', ev_size.at(0)}, {'w', ev_size.at(1)}});

	    // If the propagator is conjugated, then applied \gamma_5 left and right
	    // NOTE: \g_5 (V' D^{-1} V) \g_5 = V' \g_5 D^{-1} \g_5 V = V' D^{-\dagger} V
	    // NOTE: colorvec function gamma5Herm also adjoint the matrix; we do that just above
	    if (this_conj)
	    {
	      ti = Hadron::detail::contractSpins(Hadron::detail::contractSpins(dr_g5_left, ti),
						 dr_g5_right);
	    }
	    else
	    {
	      ti =
		Hadron::detail::contractSpins(Hadron::detail::contractSpins(dr_left, ti), dr_right);
	    }

	    SBN::copyTo(ti, prop_i);

	    // If testing, also record it as a missing meson
	    if (testing)
	      record_missing_key(key, this_conj, is_swap, i);
	  }
	}
      }

      // If testing, save the output tensor and set it to zero
      SBN::Tensor true_props;
      if (testing)
      {
	true_props = props;
	props = SBN::like_this(props);
      }

      // Recompile for each missing time slice, the phases, momenta, and displacement to compute
      using mass_tsource_source_sink_phases_conj =
	std::tuple<std::string, int, SB::Coor<3>, SB::Coor<3>, bool>;
      using tsink_vs_indices_t = Hadron::detail::unordered_multimap<int, int>;
      Hadron::detail::unordered_map<mass_tsource_source_sink_phases_conj, tsink_vs_indices_t>
	from_mass_tsource_source_sink_phases_conj_to_tsink_and_indices;
      const auto get_keys = [=](const auto& map) {
	using T = typename std::remove_reference<decltype(map)>::type::key_type;
	std::set<T> r;
	for (const auto& it : map)
	  r.insert(it.first);
	return std::vector<T>(r.begin(), r.end());
      };
      for (const auto& missing_props_in_some_process : SBN::gather(local_missing_props))
      {
	for (const auto& [prop_key, do_conj, index] : missing_props_in_some_process)
	{
	  const auto& k = mass_tsource_source_sink_phases_conj{
	    prop_key.mass_label, prop_key.t_source, ADATIO::detail::toCoor(prop_key.phasing_source),
	    ADATIO::detail::toCoor(prop_key.phasing_sink), do_conj};
	  from_mass_tsource_source_sink_phases_conj_to_tsink_and_indices[k].insert(
	    {prop_key.t_slice, index});
	}
      }

      for (const auto& it : from_mass_tsource_source_sink_phases_conj_to_tsink_and_indices)
      {
	const auto& [mass_label, t_source, source_phase, sink_phase, do_conj_] = it.first;
	const auto& do_conj = do_conj_;
	const auto& from_tsink_to_indices = it.second;
	const auto& t_sinks = get_keys(from_tsink_to_indices);

	// Get num_vecs colorvecs on time-slice t_source
	const int decay_dir = 3;
	SB::Tensor<Nd + 3, SB::Complex> source_colorvec = SB::getColorvecs<SB::Complex>(
	  colorvecsSto, u, decay_dir, t_source, 1, num_vecs, SB::none);
	source_colorvec =
	  source_colorvec.kvslice_from_size({{'n', ev_from.at(1)}}, {{'n', ev_size.at(1)}});
	source_colorvec = SB::phaseColorvecs(source_colorvec, t_source, source_phase);

	// Get num_vecs colorvecs on time-slice t_sink
	SB::Tensor<Nd + 3, SB::Complex> sinks_colorvecs =
	  source_colorvec.like_this(SB::none, {{'n', ev_size.at(0)}, {'t', t_sinks.size()}});
	for (int t_sink_index = 0; t_sink_index < t_sinks.size(); ++t_sink_index)
	{
	  SB::Tensor<Nd + 3, SB::Complex> sink_colorvec = SB::getColorvecs<SB::Complex>(
	    colorvecsSto, u, decay_dir, t_sinks.at(t_sink_index), 1, num_vecs, SB::none);
	  sink_colorvec =
	    sink_colorvec.kvslice_from_size({{'n', ev_from.at(0)}}, {{'n', ev_size.at(0)}});
	  SB::phaseColorvecs(sink_colorvec, t_sinks.at(t_sink_index), sink_phase)
	    .copyTo(sinks_colorvecs.kvslice_from_size({{'t', t_sink_index}}, {{'t', 1}}));
	}
	sinks_colorvecs = sinks_colorvecs.rename_dims({{'n', 'N'}, {'t', 'T'}});

	// Callback
	const auto call = [&](SB::Tensor<Nd + 5, SB::Complex> tensor, int sink_spin, int first_n) {
	  for (int t_sink_index = 0; t_sink_index < t_sinks.size(); ++t_sink_index)
	  {
	    const auto& r =
	      SB::contract<6>(
		sinks_colorvecs.kvslice_from_size({{'T', t_sink_index}}, {{'T', 1}}).conj(),
		tensor.kvslice_from_size({{'t', t_sinks.at(t_sink_index)}}, {{'t', 1}}), "cXxyz")
		.rename_dims({{'s', 'S'}, {'S', 's'}});
	    const auto& ti = Hadron::detail::contractSpins(
	      !do_conj ? dr_left_global : dr_g5_left_global, toTensor(r));
	    auto range = from_tsink_to_indices.equal_range(t_sinks.at(t_sink_index));
	    for (auto it = range.first; it != range.second; ++it)
	    {
	      const auto& index = it->second;
	      auto tii = SBN::relabel(ti, {{'n', 'w'}, {'N', 'v'}, {'s', 'r'}, {'S', 's'}});
	      SBN::copyTo(
		tii, SBN::slice_kv(props, {{'s', sink_spin}, {'w', first_n}, {'i', index}},
				   {{'s', 1}, {'w', SBN::get_kv_size(tii).at('w')}, {'i', 1}}));
	    }
	  }
	};

	// Do the inversions
	doInversion<SB::Complex>(get_prop(mass_label), source_colorvec, t_source,
				 !do_conj ? dr_right_chroma : dr_g5_right_chroma, max_rhs, call);
      }

      // If testing, make sure that recomputed mesons are similar to the ones got from storage
      if (testing)
      {
	for (int i = 0, n = SBN::get_kv_size(props).at('i'); i < n; ++i)
	{
	  const auto e_i = SBN::slice_kv(props, {{'i', i}}, {{'i', 1}});
	  const auto true_e_i = SBN::slice_kv(true_props, {{'i', i}}, {{'i', 1}});
	  auto diff = SBN::clone(e_i);
	  SBN::addTo(SBN::scale(true_e_i, -1), diff);
	  if (SBN::frob(true_e_i) * 1e-5 < SBN::frob(diff))
	    throw std::runtime_error("props not passing test");
	}
      }

      return props;
    }

    /// Return the genprops
    /// \param db: genprop storage
    /// \param colorvecsSto: colorvec storage
    /// \param u: original gauge field
    /// \param get_genprop: get solver from mass label
    /// \param genprop_keys: list of genprops keys
    /// \param max_rhs: maximum RHS to solve at once
    /// \param perms: list of permutations, one for each genprop key
    /// \param do_conj: list of whether to conjugate the genprop, one for each genprop key
    /// \param ev_from: first eigenvector to return for "vwx"
    /// \param ev_size: number of eigenvectors to return for "vwx"
    /// \param dist_labels: dimensions to be distributed, some of "vwxi"
    /// \param alloc: allocation for the returned tensor

    inline SBN::Tensor
    get_genprop_elementals(const ADATIO::StorageGenprop4& db,
			   const SB::ColorvecsStorage& colorvecsSto,
			   const multi1d<LatticeColorMatrix>& u,
			   const std::function<SB::ChimeraSolver(std::string)>& get_prop,
			   int max_rhs, bool zero_values_for_outside_t_slices,
			   const std::vector<Hadron::KeyGenProp4ElementalOperator_t>& genprop_keys,
			   const std::vector<SBN::Coor>& perms, const std::vector<bool>& do_conj,
			   const SBN::Coor& ev_from, const SBN::Coor& ev_size,
			   const std::string& dist_labels, const SBN::Tensor& guide, bool testing)
    {
      if (!zero_values_for_outside_t_slices)
	throw std::runtime_error("unsupported case: zero_values_for_outside_t_slices is false!");
      if (genprop_keys.size() != perms.size() || genprop_keys.size() != do_conj.size())
	throw std::runtime_error("invalid input");
      if (ev_from.size() != 4 || ev_size.size() != 4)
	throw std::runtime_error("invalid input");
      if (SBN::detail::get_debug_level() > 0)
      {
	for (const auto& it : perms)
	{
	  if (it != SBN::Coor{0, 1, 2, 3} && it != SBN::Coor{1, 0, 3, 2})
	    throw std::runtime_error("invalid input");
	}
      }

      const int num_vecs = std::max(ev_from.at(0) + ev_size.at(0), ev_from.at(1) + ev_size.at(1));
      if (max_rhs == 0)
	max_rhs = num_vecs;

      // This object is in the DR basis; create matrices to convert it to DP
      const auto& dr_left_global =
	SBN::slice_kv(Hadron::detail::adjForSpins(
			Hadron::detail::getDiracToDRMat(SBN::Options::Distribution::Replicated)), //
		      {{'s', ev_from.at(2)}}, {{'s', ev_size.at(2)}});
      const auto& dr_left = SBN::get_local_tensor(dr_left_global);
      const auto& dr_right = SBN::slice_kv(Hadron::detail::getDiracToDRMat(),
					   {{'S', ev_from.at(3)}}, {{'S', ev_size.at(3)}});

      // Create matrices to convert it to DP and with pre/post applying \gamma_5
      const auto dr_g5_left_global = Hadron::detail::contractSpins(
	dr_left_global, Hadron::detail::chromaGamma5(SBN::Options::Distribution::Replicated));
      const auto& dr_g5_left = SBN::get_local_tensor(dr_g5_left_global);

      // Create chroma versions of dr_right and dr_g5_left
      auto dr_right_chroma = SB::Tensor<2, SB::Complex>("Ss", {{Ns, ev_size.at(3)}}, SB::OnHost,
							SB::OnEveryoneReplicated);
      SBN::copyTo(SBN::relabel(dr_right, {{'S', 's'}, {'s', 'S'}}),
		  SBN::get_local_tensor(toTensor(dr_right_chroma, false /* don't copy */)));
      auto dr_g5_left_conj_chroma = SB::Tensor<2, SB::Complex>(
	"Ss", {{Ns, ev_size.at(2)}}, SB::OnHost, SB::OnEveryoneReplicated);
      SBN::copyTo(SBN::relabel(SBN::conj(dr_g5_left), {{'s', 's'}, {'S', 'S'}}),
		  SBN::get_local_tensor(toTensor(dr_g5_left_conj_chroma, false /* don't copy */)));

      const int decay_dir = 3;
      const int Lt = Layout::lattSize()[decay_dir];

      // Create output tensor
      SBN::Tensor genprops = SBN::create_tensor_with_local_components(
	SBN::concat(ev_size, {-(int)genprop_keys.size()}), "vwrsi", dist_labels,
	SBN::Options::Alloc::Host, 0, 0, SBN::Options::IsEg::False, guide);

      // Set all genprops to zero as the default value for the keys with tslice outside of t_source and t_sink
      SBN::set_zero(genprops);

      // Try to get the mesons from the storage and annotate the missing keys
      std::vector<std::tuple<Hadron::KeyGenProp4ElementalOperator_t, int>> local_missing_genprops;
      {
	local_missing_genprops.reserve(genprop_keys.size());
	const auto first_local_genprop = SBN::get_local_srange(genprops).at(0).at('i');
	const auto& local_genprops = //
	  SBN::relabel(		     //
	    SBN::slice_kv(	     //
	      SBN::get_local_tensor(genprops), {{'i', first_local_genprop}},
	      {{'i', (int)genprop_keys.size()}}),
	    {{'r', 's'}, {'s', 'S'}});
	auto val = SBN::create_tensor({ev_size.at(0), ev_size.at(1), ENSEM::Ns, ENSEM::Ns}, "vwsS",
				      SBN::Options::Distribution::Local, SBN::Options::Alloc::Host);
	auto toCoor4 = [](const SBN::Coor& coor) {
	  if (coor.size() != 4)
	    throw std::runtime_error("wtf");
	  std::array<int, 4> r;
	  std::copy_n(coor.begin(), 4, r.begin());
	  return r;
	};
	for (std::size_t i = 0; i < genprop_keys.size(); ++i)
	{
	  const auto& genprop_i = SBN::slice_kv(local_genprops, {{'i', i}}, {{'i', 1}});
	  if (!has_local_support(genprop_i))
	    continue;

	  // Check that tslice is inside t_source and t_sink
	  const auto& key = genprop_keys.at(i);
	  if (SB::normalize_coor(key.t_sink - key.t_source, Lt) <
	      SB::normalize_coor(key.t_slice - key.t_source, Lt))
	    continue;

	  bool is_swap = (perms.at(i) != SBN::Coor{0, 1, 2, 3});
	  bool this_conj = do_conj.at(i);
	  if (is_swap || this_conj)
	    throw std::runtime_error(
	      "get_genprop_elementals: unsupported permuting or conjugated genprops");
	  const auto& val_from =
	    toCoor4(SBN::get_coor(val.order, {{'v', ev_from.at(0)}, {'w', ev_from.at(1)}}, 0));
	  const auto& val_size = toCoor4(val.size);
	  const auto& val_order =
	    SBN::relabel(SBN::remap{{'v', 'N'}, {'w', 'n'}, {'s', 'q'}, {'S', 's'}}, val.order);
	  if (db.get(key, val_from, toCoor4(val.size), val_order, data(val)) != 0)
	  {
	    // We can't miss keys when testing
	    if (testing)
	    {
	      throw std::runtime_error(std::string("doing testing and missing genprops: ") +
				       ::SB::getXML(key));
	    }

	    local_missing_genprops.push_back({key, first_local_genprop + i});
	  }
	  else
	  {
	    const auto& ti =
	      Hadron::detail::contractSpins(Hadron::detail::contractSpins(dr_left, val), dr_right);

	    SBN::copyTo(ti, genprop_i);

	    // If testing, also record it as a missing genprop
	    if (testing)
	      local_missing_genprops.push_back({key, first_local_genprop + i});
	  }
	}
      }

      // If testing, save the output tensor and set it to zero
      SBN::Tensor true_genprops;
      if (testing)
      {
	true_genprops = genprops;
	genprops = SBN::like_this(genprops);
	SBN::set_zero(genprops);
      }

      // Recompile for each missing time slice, the phases, momenta, and displacement to compute
      using displacement_t = std::vector<int>;
      using mass_tsource_sink_phase_source_sink_t =
	std::tuple<std::string, int, int, SB::Coor<3>, SB::Coor<3>>;
      using mom_disp_gamma_tslice_vs_indices_t =
	Hadron::detail::unordered_multimap<SB::Coor<4>, int>;
      using moms_disps_and_gamma_vs_indices_t =
	std::tuple<Hadron::detail::unordered_map<SB::Coor<3>, int>,    // moms
		   Hadron::detail::unordered_map<displacement_t, int>, // disps
		   Hadron::detail::unordered_map<int, int>,	       // gammas
		   mom_disp_gamma_tslice_vs_indices_t // mom, disp, gamma, tslice to index
		   >;
      Hadron::detail::unordered_map<mass_tsource_sink_phase_source_sink_t,
				    moms_disps_and_gamma_vs_indices_t>
	from_mass_tsource_sink_phase_source_sink_to_moms_disps_gammas_and_tslides_vs_indices;
      const auto get_index = [=](auto& map, const auto& value) {
	const auto s = map.size();
	if (map.count(value) == 0)
	  map[value] = s;
	return map.at(value);
      };
      const auto get_vector = [=](const auto& map) {
	std::vector<typename std::remove_reference<decltype(map)>::type::key_type> r(map.size());
	for (const auto& it : map)
	  r.at(it.second) = it.first;
	return r;
      };
      for (const auto& missing_genprops_in_some_process : SBN::gather(local_missing_genprops))
      {
	for (const auto& [genprop_key, index] : missing_genprops_in_some_process)
	{
	  const auto& k = mass_tsource_sink_phase_source_sink_t{
	    genprop_key.mass, genprop_key.t_source, genprop_key.t_sink,
	    ADATIO::detail::toCoor(genprop_key.phasing_source),
	    ADATIO::detail::toCoor(genprop_key.phasing_sink)};
	  auto& v =
	    from_mass_tsource_sink_phase_source_sink_to_moms_disps_gammas_and_tslides_vs_indices[k];
	  const auto mom_index = get_index(std::get<0>(v), ADATIO::detail::toCoor(genprop_key.mom));
	  const auto disp_index = get_index(std::get<1>(v), genprop_key.displacement);
	  const auto gamma_index = get_index(std::get<2>(v), genprop_key.g);
	  std::get<3>(v).insert({{mom_index, disp_index, gamma_index, genprop_key.t_slice}, index});
	}
      }

      for (const auto& it :
	   from_mass_tsource_sink_phase_source_sink_to_moms_disps_gammas_and_tslides_vs_indices)
      {
	const auto& [mass_label, t_source_, t_sink, source_phase, sink_phase] = it.first;
	const auto& t_source = t_source_;
	const auto& moms = get_vector(std::get<0>(it.second));
	const auto& disps = get_vector(std::get<1>(it.second));
	const auto& gammas = get_vector(std::get<2>(it.second));
	const auto& mom_disp_gamma_tslice_to_index = std::get<3>(it.second);

	const int num_tslices = SB::normalize_coor(t_sink - t_source, Lt) + 1;

	const auto& PP = get_prop(mass_label);

	// Invert a source and get the time slices
	const auto& get_inv_tslice =
	  [&](int t_slice, const SB::Coor<3>& phase, const SB::Tensor<2, SB::Complex>& spins,
	      int ev_from, int ev_size) {
	    // Get num_vecs colorvecs on time-slice t_slice
	    const int decay_dir = 3;
	    SB::Tensor<Nd + 3, SB::Complex> source_colorvec = SB::getColorvecs<SB::Complex>(
	      colorvecsSto, u, decay_dir, t_slice, 1, ev_from + ev_size, SB::none);
	    source_colorvec = source_colorvec.kvslice_from_size({{'n', ev_from}}, {{'n', ev_size}});
	    source_colorvec = SB::phaseColorvecs(source_colorvec, t_slice, phase);

	    const auto order_out = "cSxyztXns";
	    SB::Tensor<Nd + 5, SB::Complex> r(
	      order_out,
	      SB::latticeSize<Nd + 5>(
		order_out,
		{{'t', num_tslices}, {'S', Ns}, {'s', spins.kvdim().at('s')}, {'n', ev_size}}),
	      SB::OnDefaultDevice);
	    const auto call = [&](SB::Tensor<Nd + 5, SB::Complex> tensor, int sink_spin,
				  int first_n) {
	      tensor.kvslice_from_size({{'t', t_source}}, {{'t', num_tslices}})
		.copyTo(r.kvslice_from_size({{'s', sink_spin}, {'n', first_n}}, {{'s', 1}}));
	    };
	    doInversion<SB::Complex>(PP, source_colorvec, t_slice, spins, max_rhs, call);
	    return r;
	  };

	// Get num_vecs colorvecs on time-slice t_source
	auto inv_src =
	  get_inv_tslice(t_source, source_phase, dr_right_chroma, ev_from.at(1), ev_size.at(1));

	// Get num_vecs colorvecs on time-slice t_sink
	auto inv_snk =
	  get_inv_tslice(t_sink, sink_phase, dr_g5_left_conj_chroma, ev_from.at(0), ev_size.at(0))
	    .rename_dims({{'n', 'N'}, {'s', 'q'}, {'S', 'Q'}});

	// Get the gamma matrices, premultiplied by g5
	const int g5 = Ns * Ns - 1;
	std::vector<SB::Tensor<2, SB::Complex>> gamma_mats;
	{
	  for (const int g : gammas)
	  {
	    SpinMatrix gmat = Gamma(g5) * (Gamma(g) * SB::SpinMatrixIdentity());
	    gamma_mats.push_back(SB::asTensorView(gmat).cloneOn<SB::Complex>(SB::OnDefaultDevice));
	  }
	}

	auto call = [&](SB::Tensor<7, SB::Complex> r_chroma_, int disp_index, int tfrom,
			int mfrom) {
	  const auto& r_chroma = r_chroma_.toComplex().template cast<SB::ComplexD>();
	  const int tsize = r_chroma.kvdim().at('t');
	  const int msize = r_chroma.kvdim().at('m');
	  const auto& r =
	    SBN::relabel(toTensor(r_chroma), {{'N', 'v'}, {'n', 'w'}, {'q', 'r'}, {'s', 's'}});
	  for (int t = 0; t < tsize; ++t)
	  {
	    for (int g = 0; g < gammas.size(); ++g)
	    {
	      for (int m = 0; m < msize; ++m)
	      {
		const auto& ti =
		  SBN::slice_kv(r, {{'g', g}, {'m', m}, {'t', t}}, {{'g', 1}, {'m', 1}, {'t', 1}});
		const int this_t = SB::normalize_coor(tfrom + t, Lt);
		const auto& k = SB::Coor<4>{mfrom + m, disp_index, g, this_t};
		auto range = mom_disp_gamma_tslice_to_index.equal_range(k);
		for (auto it = range.first; it != range.second; ++it)
		{
		  SBN::copyTo(ti, SBN::slice_kv(genprops, {{'i', it->second}}, {{'i', 1}}));
		}
	      }
	    }
	  }
	};

	const auto use_derivP = false;
	const int max_moms_in_contraction = 1;
	const int max_tslices_in_contraction = 1;
	SB::doMomGammaDisp_contractions<7, Nd + 5, Nd + 5, SB::Complex>(
	  u, inv_snk, inv_src, t_source, 0, num_tslices, moms, gamma_mats, disps, use_derivP, call,
	  "qgmNnst", max_tslices_in_contraction, max_moms_in_contraction, inv_src.getDev());
      }

      // If testing, make sure that recomputed genprops are similar to the ones got from storage
      if (testing)
      {
	for (int i = 0, n = SBN::get_kv_size(genprops).at('i'); i < n; ++i)
	{
	  const auto e_i = SBN::slice_kv(genprops, {{'i', i}}, {{'i', 1}});
	  const auto true_e_i = SBN::slice_kv(true_genprops, {{'i', i}}, {{'i', 1}});
	  auto diff = SBN::clone(e_i);
	  SBN::addTo(SBN::scale(true_e_i, -1), diff);
	  if (SBN::frob(true_e_i) * 1e-5 < SBN::frob(diff))
	    throw std::runtime_error("genprops not passing test");
	}
      }

      return genprops;
    }

    // Function call
    void InlineMeas::func(unsigned long update_no, XMLWriter& xml_out)
    {
      START_CODE();
      if (Nc != 3)
      { /* Code is specific to Ns=4 and Nc=3. */
	QDPIO::cerr << " code only works for Nc=3 and Ns=4\n";
	QDP_abort(111);
      }
#  if QDP_NC == 3

      StopWatch snoop;
      snoop.reset();
      snoop.start();

      // Test and grab a reference to the gauge field
      XMLBufferWriter gauge_xml;
      try
      {
	TheNamedObjMap::Instance().getData<multi1d<LatticeColorMatrix>>(params.named_obj.gauge_id);
	TheNamedObjMap::Instance().get(params.named_obj.gauge_id).getRecordXML(gauge_xml);
      } catch (std::bad_cast)
      {
	QDPIO::cerr << name << ": caught dynamic cast error" << std::endl;
	QDP_abort(1);
      } catch (const std::string& e)
      {
	QDPIO::cerr << name << ": std::map call failed: " << e << std::endl;
	QDP_abort(1);
      }

      // Cast should be valid now
      const multi1d<LatticeColorMatrix>& u =
	TheNamedObjMap::Instance().getData<multi1d<LatticeColorMatrix>>(params.named_obj.gauge_id);

      //
      // Read in the source along with relevant information.
      //

      SB::ColorvecsStorage colorvecsSto = SB::openColorvecStorage(params.named_obj.colorvec_files);

      push(xml_out, "CorrSuperb");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << name << ": compute correlation functions" << std::endl;

      proginfo(xml_out); // Print out basic program info

      // Write out the input
      params.writeXML(xml_out, "Input");

      // Write out the config info
      write(xml_out, "Config_info", gauge_xml);

      push(xml_out, "Output_version");
      write(xml_out, "out_version", 1);
      pop(xml_out);

      //First calculate some gauge invariant observables just for info.
      //This is really cheap.
      MesPlq(xml_out, "Observables", u);

      //
      // Smear the gauge field if needed
      //
      multi1d<LatticeColorMatrix> u_smr = u;
      try
      {
	std::istringstream xml_l(params.param.link_smearing.xml);
	XMLReader linktop(xml_l);
	QDPIO::cout << "Link smearing type = " << params.param.link_smearing.id << std::endl;

	Handle<LinkSmearing> linkSmearing(TheLinkSmearingFactory::Instance().createObject(
	  params.param.link_smearing.id, linktop, params.param.link_smearing.path));

	(*linkSmearing)(u_smr);
      } catch (const std::string& e)
      {
	QDPIO::cerr << name << ": Caught Exception link smearing: " << e << std::endl;
	QDP_abort(1);
      }

      MesPlq(xml_out, "Smeared_Observables", u_smr);

      const auto nev = params.param.num_vecs;

      std::map<char, std::string> flavor_to_mass;
      for (const auto& it : params.param.flavor_to_mass)
	flavor_to_mass[it.flavor] = it.mass;

      std::map<std::string, ChromaProp_t> mass_to_prop_options;
      for (const auto& it : params.param.flavor_to_prop)
      {
	if (flavor_to_mass.count(it.flavor) == 0)
	  flavor_to_mass[it.flavor] = std::string{'_', it.flavor};
	mass_to_prop_options[flavor_to_mass.at(it.flavor)] = it.prop;
      }

      std::map<std::string, SB::ChimeraSolver> mass_to_prop;
      const auto& get_prop = [&](const std::string& mass) {
	if (mass_to_prop.count(mass) == 0)
	{
	  if (mass_to_prop_options.count(mass) == 0)
	  {
	    QDPIO::cout << "Unspecified mass or flavor: " << mass << std::endl;
	    QDP_abort(1);
	  }
	  QDPIO::cout << "Initializing propagator for mass " << mass << std::endl;
	  const auto& prop_options = mass_to_prop_options.at(mass);
	  mass_to_prop.insert(
	    {mass, SB::ChimeraSolver{prop_options.fermact, prop_options.invParam, u}});
	}
	return mass_to_prop.at(mass);
      };

      QDPIO::cout << "Opening corr graph " << params.named_obj.corr_graph_file << std::endl;
      Hadron::CorrGraphMap_t corr_graph;
      ADATIO::BinaryFileReader bin(params.named_obj.corr_graph_file);
      read(bin, corr_graph);
      bin.close();

      // Check that the corr graph lattice size coincides with chroma's
      if (SB::tovector(QDP::Layout::lattSize()) != corr_graph.layout.latt_size)
      {
	QDPIO::cerr << "The lattice size of the corr graph file does not coincide" << std::endl;
	QDP_abort(1);
      }

      QDPIO::cout << "Opening mesons" << std::endl;
      Hadron::StorageMeson storage_meson;
      storage_meson.open(params.param.meson_files, nev);
      QDPIO::cout << "Opening baryons" << std::endl;
      Hadron::StorageBaryon storage_baryon;
      storage_baryon.open(params.param.baryon_files, nev);
      QDPIO::cout << "Opening props" << std::endl;
      Hadron::StorageProp4 storage_prop;
      storage_prop.open(params.param.prop_files, nev);
      QDPIO::cout << "Opening genprops" << std::endl;
      Hadron::StorageGenprop4 storage_genprop;
      const bool zero_values_for_outside_t_slices = true;
      storage_genprop.open(params.param.genprop_files, nev, zero_values_for_outside_t_slices);

      const auto meson_callback =
	[&](const std::vector<Hadron::KeyMesonElementalOperator_t>& meson_keys,
	    const std::vector<SBN::Coor>& perms, const std::vector<bool>& do_conj,
	    const SBN::Coor& ev_from, const SBN::Coor& ev_size, const std::string& dist_labels,
	    const SBN::Tensor& guide) {
	  return get_meson_elementals(storage_meson, colorvecsSto, u, u_smr, meson_keys, perms,
				      do_conj, ev_from, ev_size, dist_labels, guide,
				      params.param.testing);
	};
      const auto baryon_callback =
	[&](const std::vector<Hadron::KeyBaryonElementalOperator_t>& baryon_keys,
	    const std::vector<SBN::Coor>& perms, const std::vector<bool>& do_conj,
	    const SBN::Coor& ev_from, const SBN::Coor& ev_size, const std::string& dist_labels,
	    const SBN::Tensor& perm) {
	  return get_baryon_elementals(storage_baryon, colorvecsSto, u, u_smr, baryon_keys, perms,
				       do_conj, ev_from, ev_size, dist_labels, perm,
				       params.param.testing);
	};
      const auto prop_callback =
	[&](const std::vector<Hadron::KeyProp4ElementalOperator_t>& prop_keys,
	    const std::vector<SBN::Coor>& perms, const std::vector<bool>& do_conj,
	    const SBN::Coor& ev_from, const SBN::Coor& ev_size, const std::string& dist_labels,
	    const SBN::Tensor& perm) {
	  return get_prop_elementals(storage_prop, colorvecsSto, u, get_prop, params.param.max_rhs,
				     prop_keys, perms, do_conj, ev_from, ev_size, dist_labels, perm,
				     params.param.testing);
	};
      const auto genprop_callback =
	[&](const std::vector<Hadron::KeyGenProp4ElementalOperator_t>& genprop_keys,
	    const std::vector<SBN::Coor>& perms, const std::vector<bool>& do_conj,
	    const SBN::Coor& ev_from, const SBN::Coor& ev_size, const std::string& dist_labels,
	    const SBN::Tensor& perm) {
	  return get_genprop_elementals(storage_genprop, colorvecsSto, u, get_prop,
					params.param.max_rhs, zero_values_for_outside_t_slices,
					genprop_keys, perms, do_conj, ev_from, ev_size, dist_labels,
					perm, params.param.testing);
	};

      const bool zeroUnsmearedGraphsP = true;
#    if defined(QDP_IS_QDPJIT) && defined(SUPERBBLAS_USE_GPU)
      SBN::get_default_gpu_device() = SB::detail::get_default_gpu_device();
#    endif
      const auto& corr = Hadron::evaluate_graphs_with_superb(
	corr_graph, zeroUnsmearedGraphsP, prop_callback, baryon_callback, meson_callback,
	genprop_callback, flavor_to_mass, nev, params.param.t_origin, params.param.Nt_forward);

      QDPIO::cout << "Storing the correlation functions" << std::endl;
      if (Layout::nodeNumber() == 0)
      {
	const int decay_dir = 3;
	Hadron::writeCorrMap(corr, params.param.ensemble, corr_graph.layout.latt_size, decay_dir,
			     params.param.t_origin, params.named_obj.corr_file);
      }

      // Close colorvecs storage
      SB::closeColorvecStorage(colorvecsSto);

      // Close the namelist output file XMLDAT
      pop(xml_out); // CorrSuperb

      snoop.stop();
      QDPIO::cout << name << ": total time = " << snoop.getTimeInSeconds() << " secs" << std::endl;

      QDPIO::cout << name << ": ran successfully" << std::endl;
#  endif

      END_CODE();
    } // func
  }   // namespace InlineCorrSuperbEnv

  /*! @} */ // end of group hadron

} // namespace Chroma

#endif // BUILD_REDSTAR_DATALIB
