#include "meas/inline/hadron/baryon_derivative_ylm.h"
#include <fstream>
#include <tuple>
/*! \file
 * \brief Inline measurement of baryon operators via colorvector matrix elements
 */

#include "meas/inline/hadron/inline_baryon_matelem_colorvec_ylm_superb_w.h"
#include "handle.h"
#include "meas/glue/mesplq.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/make_xml_file.h"
#include "meas/smear/disp_colvec_map.h"
#include "meas/smear/link_smearing_aggregate.h"
#include "meas/smear/link_smearing_factory.h"
#include "util/ferm/key_val_db.h"
#include "util/ferm/subset_vectors.h"
#include "util/ferm/superb_contractions.h"
#include "util/info/proginfo.h"

#include "meas/inline/io/named_objmap.h"

#include <set>

#define COLORVEC_MATELEM_TYPE_ZERO 0
#define COLORVEC_MATELEM_TYPE_ONE 1
#define COLORVEC_MATELEM_TYPE_MONE -1
#define COLORVEC_MATELEM_TYPE_YLM 12

#ifdef BUILD_SB
namespace Chroma
{
  /*!
   * \ingroup hadron
   *
   * @{
   */
  namespace InlineBaryonMatElemColorVecYlmSuperbEnv
  {
    //! Phase-mom combo input
    void read(XMLReader& xml, const std::string& path,
	      phase_moms_combo_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "phase", input.phase);
      read(inputtop, "mom_list", input.mom_list);
    }

    //! Phase-mom combo output
    void write(XMLWriter& xml, const std::string& path,
	       const phase_moms_combo_t& input)
    {
      push(xml, path);

      write(xml, "phase", input.phase);
      write(xml, "mom_list", input.mom_list);

      pop(xml);
    }


    // Reader for input parameters
    void read(XMLReader& xml, const std::string& path,
	      InlineBaryonMatElemColorVecYlmSuperbEnv::Params::Param_t& param)
    {
      XMLReader paramtop(xml, path);

      if (paramtop.count("use_derivP") || paramtop.count("displacement_list"))
        throw std::runtime_error("Baryon YLM requires ylm_list and always uses derivatives");
      if (paramtop.count("drop_negative_m")) {
        bool drop = false; read(paramtop, "drop_negative_m", drop);
        if (drop) throw std::runtime_error("Baryon YLM requires all magnetic components");
      }
      param.displacement_length = 1;
      if (paramtop.count("displacement_length"))
        read(paramtop, "displacement_length", param.displacement_length);
      if (param.displacement_length != 1)
        throw std::runtime_error("Baryon YLM requires displacement_length=1");

      param.mom2_min = 0;
      if (paramtop.count("mom2_min") > 0)
      {
	read(paramtop, "mom2_min", param.mom2_min);
      }

      param.mom2_max = 0;
      if (paramtop.count("mom2_max") > 0)
      {
	read(paramtop, "mom2_max", param.mom2_max);
      }

      if (paramtop.count("mom_list") > 0)
      {
	read(paramtop, "mom_list", param.mom_list);
      }

      read(paramtop, "ylm_list", param.ylm_list);
      read(paramtop, "num_vecs", param.num_vecs);
      read(paramtop, "decay_dir", param.decay_dir);

      if (Nd != 4 || param.decay_dir != 3 || param.num_vecs <= 0)
        throw std::runtime_error("Baryon YLM requires Nd=4, decay_dir=3, num_vecs>0");
      param.t_source = 0;
      if (paramtop.count("t_source") > 0)
      {
	read(paramtop, "t_source", param.t_source);
      }

      param.Nt_forward = Layout::lattSize()[param.decay_dir];
      if (paramtop.count("Nt_forward") > 0)
      {
	read(paramtop, "Nt_forward", param.Nt_forward);
      }

      if (paramtop.count("t_slices") > 0)
      {
	read(paramtop, "t_slices", param.t_slices);
      }

      if (paramtop.count("phases") == 1)
      {
	read(paramtop, "phases", param.phases);
	for (const auto& phase : param.phases)
	{
	  if (phase.size() != Nd - 1)
	  {
	    throw std::runtime_error(std::string("each element of `phases' should have ") +
				     std::to_string(Nd - 1) + " components");
	  }
	}
      }
      else if (paramtop.count("phase") == 1)
      {
	read(paramtop, "phase", param.phase);
	if (param.phase.size() != Nd - 1)
	{
	  throw std::runtime_error(std::string("phase should have ") + std::to_string(Nd - 1) +
				   " components");
	}
      }
      else if (paramtop.count("combos") == 0)
      {
	param.phase.resize(Nd - 1);
	for (int i = 0; i < Nd - 1; ++i)
	  param.phase[i] = 0;
      }

      if (paramtop.count("combos") == 1)
      {
	read(paramtop, "combos", param.alt_phase_moms_combos);
      }

      param.max_tslices_in_contraction = 0;
      if (paramtop.count("max_tslices_in_contraction") == 1)
      {
	read(paramtop, "max_tslices_in_contraction", param.max_tslices_in_contraction);
      }

      param.max_moms_in_contraction = 0;
      if (paramtop.count("max_moms_in_contraction") == 1)
      {
	read(paramtop, "max_moms_in_contraction", param.max_moms_in_contraction);
      }

      param.max_vecs = 0;
      if (paramtop.count("max_vecs") == 1)
      {
	read(paramtop, "max_vecs", param.max_vecs);
      }

      param.use_superb_format = true;
      if( paramtop.count("use_superb_format") == 1 ) {
        read(paramtop, "use_superb_format", param.use_superb_format);
      }

      param.output_file_is_local = false;
      if( paramtop.count("output_file_is_local") == 1 ) {
        read(paramtop, "output_file_is_local", param.output_file_is_local);
      }

      if (param.max_vecs < 0 || param.max_tslices_in_contraction < 0 || param.max_moms_in_contraction < 0)
        throw std::runtime_error("Baryon YLM batch sizes must be nonnegative");
      param.link_smearing = readXMLGroup(paramtop, "LinkSmearing", "LinkSmearingType");
    }

    // Writer for input parameters
    void write(XMLWriter& xml, const std::string& path,
	       const InlineBaryonMatElemColorVecYlmSuperbEnv::Params::Param_t& param)
    {
      push(xml, path);

      write(xml, "mom2_min", param.mom2_min);
      write(xml, "mom2_max", param.mom2_max);
      write(xml, "mom_list", param.mom_list);
      write(xml, "ylm_list", param.ylm_list);
      write(xml, "derivative_basis", std::string("baryon_cg_covariant_derivatives_v1"));
      write(xml, "ylm_encoding_version", 1);
      write(xml, "angular_labels", std::string("physical_integer_L_M"));
      write(xml, "type_of_data", 12);
      write(xml, "drop_negative_m", false);
      write(xml, "displacement_length", param.displacement_length);
      write(xml, "num_vecs", param.num_vecs);
      write(xml, "decay_dir", param.decay_dir);
      write(xml, "t_source", param.t_source);
      write(xml, "Nt_forward", param.Nt_forward);
      write(xml, "t_slices", param.t_slices);
      write(xml, "phase", param.phase);
      write(xml, "phases", param.phases);
      write(xml, "combos", param.alt_phase_moms_combos);
      write(xml, "max_tslices_in_contraction", param.max_tslices_in_contraction);
      write(xml, "max_moms_in_contraction", param.max_moms_in_contraction);
      write(xml, "max_vecs", param.max_vecs);
      write(xml, "use_superb_format", param.use_superb_format);
      write(xml, "output_file_is_local", param.output_file_is_local);
      xml << param.link_smearing.xml;

      pop(xml);
    }

    //! Read named objects
    void read(XMLReader& xml, const std::string& path,
	      InlineBaryonMatElemColorVecYlmSuperbEnv::Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_id", input.gauge_id);
      read(inputtop, "colorvec_files", input.colorvec_files);
      read(inputtop, "baryon_op_file", input.baryon_op_file);
    }

    //! Write named objects
    void write(XMLWriter& xml, const std::string& path,
	       const InlineBaryonMatElemColorVecYlmSuperbEnv::Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "gauge_id", input.gauge_id);
      write(xml, "colorvec_files", input.colorvec_files);
      write(xml, "baryon_op_file", input.baryon_op_file);

      pop(xml);
    }

    // Writer for input parameters
    void write(XMLWriter& xml, const std::string& path,
	       const InlineBaryonMatElemColorVecYlmSuperbEnv::Params& param)
    {
      param.writeXML(xml, path);
    }
  }

  namespace InlineBaryonMatElemColorVecYlmSuperbEnv
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

    const std::string name = "BARYON_MATELEM_COLORVEC_YLM_SUPERB";

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

    //! Anonymous namespace
    /*! Diagnostic stuff */
    namespace
    {
      StandardOutputStream& operator<<(StandardOutputStream& os, const multi1d<int>& d)
      {
	if (d.size() > 0)
	{
	  os << d[0];

	  for (int i = 1; i < d.size(); ++i)
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
      param.mom2_min = 0;
      param.mom2_max = 0;
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

    //----------------------------------------------------------------------------
    //! Baryon operator
    struct KeyBaryonElementalOperator_t {
      int t_slice;	   /*!< Baryon operator time slice */
      multi1d<int> left;   /*!< Displacement dirs of left colorstd::vector */
      multi1d<int> middle; /*!< Displacement dirs of middle colorstd::vector */
      multi1d<int> right;  /*!< Displacement dirs of right colorstd::vector */
      multi1d<int> mom;	   /*!< D-1 momentum of this operator */
    };

    //! Baryon operator
    /*!< Momentum projected operator */
    struct ValBaryonElementalOperator_t : public SB::Tensor<3, SB::ComplexD> {
      int type_of_data; /*!< Flag indicating type of data (maybe trivial) */
      ValBaryonElementalOperator_t(int n = 0, int type_of_data = COLORVEC_MATELEM_TYPE_YLM)
	: SB::Tensor<3, SB::ComplexD>("kji", {n, n, n}, SB::OnHost, SB::Local),
	  type_of_data(type_of_data)
      {
      }
    };

    //----------------------------------------------------------------------------
    //! BaryonElementalOperator reader
    void read(BinaryReader& bin, KeyBaryonElementalOperator_t& param)
    {
      read(bin, param.t_slice);
      read(bin, param.left);
      read(bin, param.middle);
      read(bin, param.right);
      read(bin, param.mom);
    }

    //! BaryonElementalOperator write
    void write(BinaryWriter& bin, const KeyBaryonElementalOperator_t& param)
    {
      write(bin, param.t_slice);
      write(bin, param.left);
      write(bin, param.middle);
      write(bin, param.right);
      write(bin, param.mom);
    }

    //! BaryonElementalOperator reader
    void read(XMLReader& xml, const std::string& path, KeyBaryonElementalOperator_t& param)
    {
      XMLReader paramtop(xml, path);

      read(paramtop, "t_slice", param.t_slice);
      read(paramtop, "left", param.left);
      read(paramtop, "middle", param.middle);
      read(paramtop, "right", param.right);
      read(paramtop, "mom", param.mom);
    }

    //! BaryonElementalOperator writer
    void write(XMLWriter& xml, const std::string& path, const KeyBaryonElementalOperator_t& param)
    {
      push(xml, path);

      write(xml, "t_slice", param.t_slice);
      write(xml, "left", param.left);
      write(xml, "middle", param.middle);
      write(xml, "right", param.right);
      write(xml, "mom", param.mom);

      pop(xml);
    }

    //----------------------------------------------------------------------------
    //! BaryonElementalOperator reader
    void read(BinaryReader& bin, ValBaryonElementalOperator_t& param)
    {
      int type_of_data;
      read(bin, type_of_data);

      int n;
      read(bin, n); // the size is always written, even if 0
      param = ValBaryonElementalOperator_t(n, type_of_data);
      SB::Tensor<3, SB::ComplexD>& t = param;
      read(bin, t);
    }

    //! BaryonElementalOperator write
    void write(BinaryWriter& bin, const ValBaryonElementalOperator_t& param)
    {
      int type_of_data = param.type_of_data;
      write(bin, type_of_data);

      auto kvdim = param.kvdim();
      assert(kvdim['i'] == kvdim['j'] && kvdim['j'] == kvdim['k']);
      int n = kvdim['i']; // all sizes the same
      write(bin, n);
      SB::Tensor<3, SB::ComplexD> t = param.reorder("kji");
      write(bin, t);
    }

    //----------------------------------------------------------------------------
    // Function call
    void InlineMeas::operator()(unsigned long update_no, XMLWriter& xml_out)
    {
      // If xml file not empty, then use alternate
      if (params.xml_file != "")
      {
	std::string xml_file = makeXMLFileName(params.xml_file, update_no);

	push(xml_out, "BaryonMatElemColorYlmSuperbVec");
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

      push(xml_out, "BaryonMatElemColorVecSuperb");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << name << ": Baryon color-std::vector matrix element" << std::endl;

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
      // If a list of momenta has been specified only need phases corresponding to these
      //
      SB::CoorMoms mom_list;
      if (params.param.mom_list.size() == 0)
      {
	mom_list = SB::getMomenta(params.param.mom2_min, params.param.mom2_max);
      }
      else if (params.param.mom_list.size() > 0)
      {
	mom_list = SB::getMomenta(params.param.mom_list);
      }

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

      //
      // Parse the phase
      //
      std::vector<SB::Coor<Nd - 1>> phasings;
      if (params.param.phase.size() > 0)
      {
	phasings.push_back(SB::toCoor(params.param.phase));
      }
      else if (params.param.phases.size() > 0)
      {
	std::set<SB::Coor<Nd - 1>> phasings_set;
	for (const auto& phase : params.param.phases)
	  phasings_set.insert(SB::toCoor(phase));
	if (phasings_set.size() > 1 && !params.param.use_superb_format)
	  throw std::runtime_error("Unsupported computing several phasings with the old database "
				   "format; set `use_superb_format' to true");
	phasings = std::vector<SB::Coor<Nd - 1>>(phasings_set.begin(), phasings_set.end());
      }
      else if (params.param.alt_phase_moms_combos.size() == 0)
      {
	phasings.push_back(SB::Coor<Nd - 1>{{}});
      }

      // Create the phase-momenta combos
      std::map<SB::Coor<Nd - 1>, std::set<SB::Coor<Nd - 1>>>
	phase_moms_combos; ///< map from phase to moms
      if (phasings.size() > 0 && mom_list.size() > 0)
      {
	for (const auto& phase : phasings)
	  phase_moms_combos[phase].insert(mom_list.begin(), mom_list.end());
      }
      for (const auto& it : params.param.alt_phase_moms_combos)
      {
	const auto& moms = SB::getMomenta(it.mom_list);
	phase_moms_combos[SB::toCoor(it.phase)].insert(moms.begin(), moms.end());
      }

      // Create a list of unique phases and momenta (used by superb format storage)
      std::vector<SB::Coor<Nd - 1>> unique_phases;
      std::vector<SB::Coor<Nd - 1>> unique_moms;
      {
	std::set<SB::Coor<Nd - 1>> unique_moms_s;
	for (const auto& it : phase_moms_combos)
	{
	  unique_phases.push_back(it.first);
	  unique_moms_s.insert(it.second.begin(), it.second.end());
	}
	unique_moms = std::vector<SB::Coor<Nd - 1>>(unique_moms_s.begin(), unique_moms_s.end());
      }

      std::vector<BaryonDerivativeYlm::Key> requests;
      for (int i=0; i<params.param.ylm_list.size(); ++i) {
        BaryonDerivativeYlm::Key key;
        for (int j=0; j<params.param.ylm_list[i].size(); ++j) key.push_back(params.param.ylm_list[i][j]);
        requests.push_back(key);
      }
      const auto components = BaryonDerivativeYlm::expand(requests);
      std::vector<BaryonDerivativeYlm::Path> cartesian_paths;
      std::map<BaryonDerivativeYlm::Path,int> path_index;
      std::vector<std::vector<std::pair<int,std::complex<double>>>> consumers;
      multi1d<multi1d<int>> resolved(components.size());
      for (int c=0; c<int(components.size()); ++c) {
        resolved[c] = SB::tomulti1d(components[c].key);
        for (const auto& term : components[c].paths) {
          auto ins = path_index.emplace(term.first, int(cartesian_paths.size()));
          if (ins.second) {cartesian_paths.push_back(term.first); consumers.emplace_back();}
          consumers[ins.first->second].push_back({c,term.second});
        }
      }
      write(xml_out, "ResolvedYlmList", resolved);
      if (unique_moms.empty() || unique_phases.empty())
        throw std::runtime_error("Baryon YLM requires nonempty momenta and phasings");
      if (!params.param.use_superb_format) {
        std::string exists;
        if (Layout::nodeNumber()==0) {
          std::ifstream f(params.named_obj.baryon_op_file.c_str(), std::ios::binary);
          if (f.good()) exists="exists";
        }
        if (!SB::broadcast(exists).empty())
          throw std::runtime_error("Baryon YLM FileDB output exists; choose a new filename");
      }


      // Compute the interval of t points to compute
      int tfrom = 0; // First t-slice to compute
      int tsize = 0; // Number of t-slices to compute
      std::set<int> t_slices_to_write{};
      const int Nt = Layout::lattSize()[params.param.decay_dir];
      if (params.param.t_slices.size() == 0)
      {
	tfrom = params.param.t_source;
	tsize = params.param.Nt_forward;
        if (tfrom < 0 || tfrom >= Nt || tsize <= 0 || tsize > Nt)
          throw std::runtime_error("Baryon YLM invalid time interval");
	params.param.t_slices.resize(tsize);
	for (int i = 0; i < tsize; ++i)
	{
	  params.param.t_slices[i] = (tfrom + i) % Nt;
	  t_slices_to_write.insert((tfrom + i) % Nt);
	}
      }
      else
      {
	for (int i = 0; i < params.param.t_slices.size(); ++i)
	{
	  // Check the values
	  int t0 = params.param.t_slices[i];
	  if (t0 < 0 || t0 >= Nt)
	    throw std::runtime_error("Invalid source at tag `t_slices'");

	  SB::union_interval(tfrom, tsize, t0, 1, Nt, tfrom, tsize);
	  t_slices_to_write.insert(t0);
	}
      }

      //
      // DB storage
      // NOTE: Only the master node opens the storage and writes on it
      //
      std::vector<LocalBinaryStoreDB<LocalSerialDBKey<KeyBaryonElementalOperator_t>,
				     LocalSerialDBData<ValBaryonElementalOperator_t>>>
	qdp_db{};

      // This function opens the output file
      // NOTE: Only called by the master node
      std::function<void()> open_db = [&]() {
	// If the qdp_db is already opened, do nothing
	if (qdp_db.size() > 0)
	  return;
	qdp_db.resize(1);

	// Open the file, and write the meta-data and the binary for this operator
	if (!qdp_db[0].fileExists(params.named_obj.baryon_op_file))
	{
	  XMLBufferWriter file_xml;

	  push(file_xml, "DBMetaData");
	  write(file_xml, "id", std::string("baryonElemOpYlm"));
	  write(file_xml, "lattSize", QDP::Layout::lattSize());
	  write(file_xml, "decay_dir", params.param.decay_dir);
	  proginfo(file_xml); // Print out basic program info
	  write(file_xml, "Params", params.param);
	  write(file_xml, "ylm_components", resolved);
	  write(file_xml, "Config_info", gauge_xml);

	  // Some tasks read the eigenvalues from metadata but they not used; so we are going to give fake values
	  multi1d<multi1d<double>> evals(params.param.num_vecs);
	  const int Nt = Layout::lattSize()[params.param.decay_dir];
	  for (int i = 0; i < params.param.num_vecs; ++i)
	  {
	    evals[i].resize(Nt);
	    for (int t = 0; t < Nt; ++t)
	      evals[i][t] = 0;
	  }
	  write(file_xml, "Weights", evals);

	  pop(file_xml);

	  std::string file_str(file_xml.str());
	  qdp_db[0].setMaxUserInfoLen(file_str.size());

	  qdp_db[0].open(params.named_obj.baryon_op_file, O_RDWR | O_CREAT, 0664);

	  qdp_db[0].insertUserdata(file_str);
	}
	else
	{
	  qdp_db[0].open(params.named_obj.baryon_op_file, O_RDWR, 0664);
	}
      };

      if (!params.param.use_superb_format && phase_moms_combos.size() > 1)
      {
	throw std::runtime_error("unsupported to write multiple phasings on the same file with the "
				 "old file format: set the xml tag use_superb_format to true");
      }

      /// Superb storage; dimension labels ijktdm:
      /// i,j,k: eigenvector indices
      /// t: time slice
      /// d: explicit coupled derivative descriptor in ylm_components
      /// m: momentum
      /// h: phasing

      SB::StorageTensor<7, SB::ComplexD> st;
      if (params.param.use_superb_format)
      {
	const char* order = "ijktdmh";
	XMLBufferWriter metadata_xml;
	push(metadata_xml, "DBMetaData");
	write(metadata_xml, "id", std::string("baryonElemOpYlmSuperb"));
	write(metadata_xml, "lattSize", QDP::Layout::lattSize());
	write(metadata_xml, "decay_dir", params.param.decay_dir);
	proginfo(metadata_xml); // Print out basic program info
	write(metadata_xml, "Config_info", gauge_xml);
	write(metadata_xml, "Params", params.param);
	write(metadata_xml, "tensorOrder", order);
	write(metadata_xml, "ylm_components", resolved);
	std::vector<multi1d<int>>  moms;
	for (const auto& mom : unique_moms)
	  moms.push_back(SB::tomulti1d(mom));
	write(metadata_xml, "moms", moms);
	std::vector<multi1d<int>> phasings;
	for (const auto& it : unique_phases)
	  phasings.push_back(SB::tomulti1d(it));
	write(metadata_xml, "phasings", phasings);

	// Some tasks read the eigenvalues from metadata but they not used; so we are going to give fake values
	multi1d<multi1d<double>> evals(params.param.num_vecs);
	const int Nt = Layout::lattSize()[params.param.decay_dir];
	for (int i = 0; i < params.param.num_vecs; ++i)
	{
	  evals[i].resize(Nt);
	  for (int t = 0; t < Nt; ++t)
	    evals[i][t] = 0;
	}
	write(metadata_xml, "Weights", evals);

	pop(metadata_xml);

	// NOTE: metadata_xml only has a valid value on Master node; so do a broadcast
	std::string metadata = SB::broadcast(metadata_xml.str());

	st = SB::StorageTensor<7, SB::ComplexD>(
	  params.named_obj.baryon_op_file, metadata, order,
	  SB::kvcoors<7>(order, {{'i', params.param.num_vecs},
				 {'j', params.param.num_vecs},
				 {'k', params.param.num_vecs},
				 {'t', Nt},
				 {'d', components.size()},
				 {'m', moms.size()},
				 {'h', phasings.size()}}),
	  SB::Sparse, SB::checksum_type::BlockChecksum,
	  params.param.output_file_is_local ? SB::LocalFSFile : SB::SharedFSFile);
	int num_moms = 0;
	for (const auto& it : phase_moms_combos)
	  num_moms += it.second.size();
	st.preallocate(params.param.num_vecs * params.param.num_vecs * params.param.num_vecs *
		       t_slices_to_write.size() * components.size() * num_moms *
		       sizeof(SB::ComplexD) /
		       (params.param.output_file_is_local ? Layout::numNodes() : 1));
      }


      //
      // Baryon operators
      //
      // Derivative/color tensors only. Redstar supplies spin, flavor,
      // permutation coupling and creation/annihilation conventions.
      //
      QDPIO::cout << "Building baryon operators" << std::endl;

      push(xml_out, "ElementalOps");

      // Build the operator
      StopWatch swiss;
      swiss.reset();
      swiss.start();

      double time_storing = 0; // total time in writing elementals

      // NOTE: st needs MPI synchronization when closing, so capture exception and abort in that case
      //       to avoid hangs
      try
      {
	int max_tslices_in_contraction = params.param.max_tslices_in_contraction == 0
					   ? tsize
					   : params.param.max_tslices_in_contraction;

	// Iterate over time-slices
	for (int tfrom0 = 0, this_tsize = std::min(tsize, max_tslices_in_contraction);
	     tfrom0 < tsize; tfrom0 += this_tsize,
		 this_tsize = std::min(max_tslices_in_contraction, tsize - tfrom0))
	{
	  int this_tfrom = (tfrom + tfrom0) % Nt;

	  // Get num_vecs colorvecs on time-slice t_source
	  SB::Tensor<Nd + 3, SB::Complex> colorvec =
	    SB::getColorvecs<SB::Complex>(colorvecsSto, u, params.param.decay_dir, this_tfrom,
					  this_tsize, params.param.num_vecs, "cxyzXnt");

	  for (int phase_idx = 0; phase_idx < unique_phases.size(); ++phase_idx)
	  {
	    SB::Tensor<Nd + 3, SB::Complex> source_colorvec =
	      SB::phaseColorvecs(colorvec, this_tfrom, unique_phases.at(phase_idx));
	    // Call for storing the baryons
	    const auto& moms_set = phase_moms_combos.at(unique_phases.at(phase_idx));
	    const auto& moms = SB::CoorMoms(moms_set.begin(), moms_set.end());
	    auto store_component = [&](SB::Tensor<5, SB::ComplexD> tensor, int disp, int first_tslice, int first_mom) {
		StopWatch tstoring;
		tstoring.reset();
		tstoring.start();

		if (params.param.use_superb_format)
		{
		  for (int m = 0, numm = tensor.kvdim()['m']; m < numm; ++m)
		  {
		    int mom_idx =
		      std::find(unique_moms.begin(), unique_moms.end(), moms.at(first_mom + m)) -
		      unique_moms.begin();
		    if (mom_idx == unique_moms.size())
		      throw std::runtime_error("this shouldn't happen");
		    for (int t = 0, numt = tensor.kvdim()['t']; t < numt; ++t)
		    {
		      if (t_slices_to_write.count((first_tslice + t) % Nt) == 0)
			continue;

		      st.kvslice_from_size({{'t', (first_tslice + t) % Nt},
					    {'d', disp},
					    {'m', mom_idx},
					    {'h', phase_idx}},
					   {{'t', 1}, {'d', 1}, {'m', 1}, {'h', 1}})
			.copyFrom(
			  tensor.kvslice_from_size({{'t', t}, {'m', m}}, {{'t', 1}, {'m', 1}}));
		    }
		  }
		}
		else
		{
		  // Only the master node writes the elementals and we assume that tensor is only supported on master
		  assert(tensor.dist == SB::OnMaster);
		  tensor = tensor.getLocal();
		  if (tensor) // if the local tensor isn't empty, ie this node holds the tensor
		  {
		    // Open the database
		    open_db();

		    KeyBaryonElementalOperator_t key;
		    ValBaryonElementalOperator_t val(params.param.num_vecs);

		    for (int t = 0, numt = tensor.kvdim()['t']; t < numt; ++t)
		    {
		      if (t_slices_to_write.count((first_tslice + t) % Nt) == 0)
			continue;
		      for (int m = 0, numm = tensor.kvdim()['m']; m < numm; ++m)
		      {
			key.t_slice = (first_tslice + t) % Nt;
			key.left.resize(0);
			key.middle.resize(0);
			key.right = SB::tomulti1d(components[disp].key);
			key.mom = SB::tomulti1d(moms.at(first_mom + m));
			tensor.kvslice_from_size({{'t', t}, {'m', m}}, {{'t', 1}, {'m', 1}})
			  .copyTo(val);
			qdp_db[0].insert(key, val);
		      }
		    }
		  }
		}

		tstoring.stop();
		time_storing += tstoring.getTimeInSeconds();
              };
            struct Pending { SB::Tensor<5, SB::ComplexD> value; int count=0; };
            std::map<std::tuple<int,int,int>,Pending> pending;
            SB::ColorContractionFn<SB::Complex> call(
              [&](SB::Tensor<5, SB::Complex> tensor, int path, int first_tslice, int first_mom) {
                auto term = tensor.template like_this<5, SB::ComplexD>();
                tensor.copyTo(term);
                for (const auto& target : consumers.at(path)) {
                  auto key=std::make_tuple(first_tslice,first_mom,target.first);
                  auto& sum=pending[key];
                  if (sum.count==0) {
                    sum.value=tensor.template like_this<5, SB::ComplexD>();
                    sum.value.set_zero();
                  }
                  term.scale(target.second).addTo(sum.value);
                  if (++sum.count==int(components[target.first].paths.size())) {
                    store_component(sum.value,target.first,first_tslice,first_mom);
                    pending.erase(key);
                  }
                }
              });

	    // Do the color-contraction
	    SB::doMomDisp_colorContractions(
	      u_smr, source_colorvec, moms, this_tfrom,
	      cartesian_paths, true, call, 0 /* it means to do all */,
	      params.param.max_moms_in_contraction, params.param.max_vecs, SB::none,
	      SB::OnDefaultDevice,
	      params.param.use_superb_format ? SB::none
					     : SB::Maybe<SB::Distribution>(SB::OnMaster));
            if (!pending.empty()) throw std::runtime_error("Incomplete baryon YLM accumulation");
	  }
	}
      } catch (const std::exception& e)
      {
	std::cerr << "caught error: " << e.what() << std::endl;
	QDP_abort(1);
      }

      // Close db
      for (auto& db : qdp_db)
	db.close();

      // Close colorvecs storage
      SB::closeColorvecStorage(colorvecsSto);

      swiss.stop();

      QDPIO::cout << "All baryon operators computed in time= "
		  << swiss.getTimeInSeconds() - time_storing << " secs, writing time "
		  << time_storing << " secs " << std::endl;
      pop(xml_out); // ElementalOps

      // Close the namelist output file XMLDAT
      pop(xml_out); // BaryonMatElemColorVector

      snoop.stop();
      QDPIO::cout << name << ": total time = " << snoop.getTimeInSeconds() << " secs" << std::endl;

      QDPIO::cout << name << ": ran successfully" << std::endl;

#  endif

      END_CODE();
    } // func
  }   // namespace InlineBaryonMatElemColorVecYlmSuperbEnv

  /*! @} */ // end of group hadron

} // namespace Chroma

#endif // BUILD_SB
