/*! \file
 * \brief Inline measurement of baryon operators via colorstd::vector matrix elements
 */

#include "meas/inline/hadron/inline_baryon_matelem_colorvec_superb_w.h"
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
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"

#include "meas/inline/io/named_objmap.h"

#define COLORVEC_MATELEM_TYPE_ZERO 0
#define COLORVEC_MATELEM_TYPE_ONE 1
#define COLORVEC_MATELEM_TYPE_MONE -1
#define COLORVEC_MATELEM_TYPE_GENERIC 10

#ifdef BUILD_SB
namespace Chroma
{
  /*!
   * \ingroup hadron
   *
   * @{
   */
  namespace InlineBaryonMatElemColorVecSuperbEnv
  {
    // Reader for input parameters
    void read(XMLReader& xml, const std::string& path,
	      InlineBaryonMatElemColorVecSuperbEnv::Params::Param_t::Displacement_t& param)
    {
      XMLReader paramtop(xml, path);

      read(paramtop, "left", param.left);
      read(paramtop, "middle", param.middle);
      read(paramtop, "right", param.right);
    }

    // Writer for input parameters
    void write(XMLWriter& xml, const std::string& path,
	       const InlineBaryonMatElemColorVecSuperbEnv::Params::Param_t::Displacement_t& param)
    {
      push(xml, path);

      write(xml, "left", param.left);
      write(xml, "middle", param.middle);
      write(xml, "right", param.right);

      pop(xml);
    }

    // Reader for input parameters
    void read(XMLReader& xml, const std::string& path,
	      InlineBaryonMatElemColorVecSuperbEnv::Params::Param_t& param)
    {
      XMLReader paramtop(xml, path);

      int version;
      read(paramtop, "version", version);

      switch (version)
      {
      case 1: param.use_derivP = false; break;

      case 2: read(paramtop, "use_derivP", param.use_derivP); break;

      default:
	QDPIO::cerr << "Input parameter version " << version << " unsupported." << std::endl;
	QDP_abort(1);
      }

      read(paramtop, "mom2_max", param.mom2_max);
      read(paramtop, "displacement_length", param.displacement_length);
      read(paramtop, "displacement_list", param.displacement_list);
      read(paramtop, "num_vecs", param.num_vecs);
      read(paramtop, "decay_dir", param.decay_dir);
      if (paramtop.count("t_source") > 0)
      {
	read(paramtop, "t_source", param.t_source);
      }
      else
      {
	param.t_source = 0;
      }

      if (paramtop.count("Nt_forward") > 0)
      {
	read(paramtop, "Nt_forward", param.Nt_forward);
      }
      else
      {
	param.Nt_forward = Layout::lattSize()[param.decay_dir];
      }
      if (paramtop.count("t_slices") > 0)
      {
	read(paramtop, "t_slices", param.t_slices);
      }

      if (paramtop.count("phase") == 1)
      {
	read(paramtop, "phase", param.phase);
      }
      else
      {
	param.phase.resize(Nd - 1);
	for (int i = 0; i < Nd - 1; ++i)
	  param.phase[i] = 0;
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

      param.link_smearing = readXMLGroup(paramtop, "LinkSmearing", "LinkSmearingType");
    }

    // Writer for input parameters
    void write(XMLWriter& xml, const std::string& path,
	       const InlineBaryonMatElemColorVecSuperbEnv::Params::Param_t& param)
    {
      push(xml, path);

      int version = 2;

      write(xml, "version", version);
      write(xml, "mom2_max", param.mom2_max);
      write(xml, "use_derivP", param.use_derivP);
      write(xml, "displacement_length", param.displacement_length);
      write(xml, "displacement_list", param.displacement_list);
      write(xml, "num_vecs", param.num_vecs);
      write(xml, "decay_dir", param.decay_dir);
      write(xml, "t_source", param.t_source);
      write(xml, "Nt_forward", param.Nt_forward);
      write(xml, "t_slices", param.t_slices);
      write(xml, "phase", param.phase);
      write(xml, "max_tslices_in_contraction", param.max_tslices_in_contraction);
      write(xml, "max_moms_in_contraction", param.max_moms_in_contraction);
      xml << param.link_smearing.xml;

      pop(xml);
    }

    //! Read named objects
    void read(XMLReader& xml, const std::string& path,
	      InlineBaryonMatElemColorVecSuperbEnv::Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_id", input.gauge_id);
      read(inputtop, "colorvec_files", input.colorvec_files);
      read(inputtop, "baryon_op_file", input.baryon_op_file);
    }

    //! Write named objects
    void write(XMLWriter& xml, const std::string& path,
	       const InlineBaryonMatElemColorVecSuperbEnv::Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "gauge_id", input.gauge_id);
      write(xml, "colorvec_files", input.colorvec_files);
      write(xml, "baryon_op_file", input.baryon_op_file);

      pop(xml);
    }

    // Writer for input parameters
    void write(XMLWriter& xml, const std::string& path,
	       const InlineBaryonMatElemColorVecSuperbEnv::Params& param)
    {
      param.writeXML(xml, path);
    }
  }

  namespace InlineBaryonMatElemColorVecSuperbEnv
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

    const std::string name = "BARYON_MATELEM_COLORVEC_SUPERB";

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

      StandardOutputStream& operator<<(StandardOutputStream& os,
				       const Params::Param_t::Displacement_t& d)
      {
	os << "left= " << d.left << " middle= " << d.middle << " right= " << d.right;

	return os;
      }
    }

    //----------------------------------------------------------------------------
    // Param stuff
    Params::Params()
    {
      frequency = 0;
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
      ValBaryonElementalOperator_t(int n = 0, int type_of_data = COLORVEC_MATELEM_TYPE_GENERIC)
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
    //! Normalize just one displacement array and return std::vector
    multi1d<int> tomulti1d(const std::vector<int>& orig)
    {
      multi1d<int> r(orig.size());
      for (int i = 0; i < orig.size(); ++i)
	r[i] = orig[i];
      return r;
    }

    //----------------------------------------------------------------------------
    //! Normalize just one displacement array and return std::vector
    std::vector<int> normDisp(const multi1d<int>& orig)
    {
      std::vector<int> r; // path to return

      // NO: a no-displacement is recorded as a zero-length array
      // Convert a length one array with no displacement into a no-displacement array
      if (orig.size() == 1 && orig[0] == 0)
	return r;

      r.resize(orig.size());
      for (int i = 0; i < orig.size(); ++i)
	r[i] = orig[i];

      return r;
    }

    //----------------------------------------------------------------------------
    //! Make sure displacements are something sensible and return std structures
    std::vector<std::array<std::vector<int>, 3>>
    normalizeDisplacements(const multi1d<Params::Param_t::Displacement_t>& displacement_list)
    {
      std::vector<std::array<std::vector<int>, 3>> r(displacement_list.size());
      for (int i = 0; i < displacement_list.size(); ++i)
      {
	r[i][0] = normDisp(displacement_list[i].left);
	r[i][1] = normDisp(displacement_list[i].middle);
	r[i][2] = normDisp(displacement_list[i].right);
      }
      return r;
    }

    //-------------------------------------------------------------------------------
    // Function call
    void InlineMeas::operator()(unsigned long update_no, XMLWriter& xml_out)
    {
      // If xml file not empty, then use alternate
      if (params.xml_file != "")
      {
	std::string xml_file = makeXMLFileName(params.xml_file, update_no);

	push(xml_out, "BaryonMatElemColorSuperbVec");
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
      QDPIO::cout << "Snarf the source from a std::map object disk file" << std::endl;

      SB::MODS_t eigen_source;
      eigen_source.setDebug(0);

      std::string eigen_meta_data; // holds the eigenvalues

      try
      {
	// Open
	QDPIO::cout << "Open file= " << params.named_obj.colorvec_files[0] << std::endl;
	eigen_source.open(params.named_obj.colorvec_files);

	// Snarf the source info.
	QDPIO::cout << "Get user data" << std::endl;
	eigen_source.getUserdata(eigen_meta_data);
      } catch (std::bad_cast)
      {
	QDPIO::cerr << name << ": caught dynamic cast error" << std::endl;
	QDP_abort(1);
      } catch (const std::string& e)
      {
	QDPIO::cerr << name << ": error extracting source_header: " << e << std::endl;
	QDP_abort(1);
      } catch (const char* e)
      {
	QDPIO::cerr << name << ": Caught some char* exception:" << std::endl;
	QDPIO::cerr << e << std::endl;
	QDPIO::cerr << "Rethrowing" << std::endl;
	throw;
      }

      QDPIO::cout << "Source successfully read and parsed" << std::endl;

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
      // Initialize the slow Fourier transform phases
      //
      SftMom phases(params.param.mom2_max, false, params.param.decay_dir);

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
      if (params.param.phase.size() != Nd - 1)
      {
	QDPIO::cerr << "phase tag should have " << Nd - 1 << " components" << std::endl;
	QDP_abort(1);
      }
      SB::Coor<Nd - 1> phase;
      for (int i = 0; i < Nd - 1; ++i)
      {
	phase[i] = params.param.phase[i];
	if (std::fabs(phase[i] - params.param.phase[i]) > 0)
	  std::runtime_error("phase should be integer");
      }

      //
      // DB storage
      // NOTE: Only the master node opens the storage and writes on it
      //
      std::vector<BinaryStoreDB<LocalSerialDBKey<KeyBaryonElementalOperator_t>,
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
	  write(file_xml, "id", std::string("baryonElemOp"));
	  write(file_xml, "lattSize", QDP::Layout::lattSize());
	  write(file_xml, "decay_dir", params.param.decay_dir);
	  proginfo(file_xml); // Print out basic program info
	  write(file_xml, "Params", params.param);
	  //write(file_xml, "Op_Info", displacement_list);
	  write(file_xml, "Config_info", gauge_xml);
	  //write(file_xml, "Weights", getEigenValues(eigen_source, params.param.num_vecs));
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

      // Compute the interval of t points to compute
      int tfrom = 0; // First t-slice to compute
      int tsize = 0; // Number of t-slices to compute
      std::set<int> t_slices_to_write{};
      const int Nt = Layout::lattSize()[params.param.decay_dir];
      if (params.param.t_slices.size() == 0)
      {
	tfrom = params.param.t_source;
	tsize = params.param.Nt_forward;
	params.param.t_slices.resize(tsize);
	for (int i = 0; i < tsize; ++i)
	  t_slices_to_write.insert((tfrom + i) % Nt);
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

      // The function SB::doMomDisp_colorContractions constrains that the number of t points to compute
      // should be zero, one, or even
      if (tsize > 1 && tsize % 2 == 1)
	tsize++;

      //
      // Baryon operators
      //
      // The creation and annihilation operators are the same without the
      // spin matrices.
      //
      QDPIO::cout << "Building baryon operators" << std::endl;

      push(xml_out, "ElementalOps");

      // Build the operator
      StopWatch swiss;
      swiss.reset();
      swiss.start();

      // Make sure displacements are something sensible and transform to std objects
      std::vector<std::array<std::vector<int>, 3>> displacement_list =
	normalizeDisplacements(params.param.displacement_list);

      // Get num_vecs colorvecs on time-slice t_source
      SB::Tensor<Nd + 3, SB::Complex> source_colorvec = SB::getColorvecs<SB::Complex>(
	eigen_source, params.param.decay_dir, tfrom, tsize, params.param.num_vecs, SB::none, phase);

      // Call to that stores the baryons
      double time_storing = 0; // total time in writing elementals
      SB::ColorContractionFn<SB::Complex> call(
	[&](SB::Tensor<5, SB::Complex> tensor, int disp, int first_tslice, int first_mom) {
	  // Only the master node writes the elementals and we assume that tensor is only supported on master
	  assert(tensor.dist == SB::OnMaster);
	  tensor = tensor.getLocal();
	  if (tensor) // if the local tensor isn't empty, ie this node holds the tensor
	  {
	    // Open the database
	    open_db();

	    StopWatch tstoring;
	    tstoring.reset();
	    tstoring.start();

	    KeyBaryonElementalOperator_t key;
	    ValBaryonElementalOperator_t val(params.param.num_vecs);

	    // The keys for the displacements for this particular elemental operator
	    // Invert the time - make it an independent key
	    for (int t = 0, numt = tensor.kvdim()['t']; t < numt; ++t)
	    {
	      if (t_slices_to_write.count(first_tslice + t) == 0)
		continue;
	      for (int m = 0, numm = tensor.kvdim()['m']; m < numm; ++m)
	      {
		key.t_slice = first_tslice + t;
		key.left = tomulti1d(displacement_list[disp][0]);
		key.middle = tomulti1d(displacement_list[disp][1]);
		key.right = tomulti1d(displacement_list[disp][2]);
		key.mom = phases.numToMom(first_mom + m);
		tensor.kvslice_from_size({{'t', t}, {'m', m}}, {{'t', 1}, {'m', 1}}).copyTo(val);
		qdp_db[0].insert(key, val);
	      }
	    }

	    time_storing += tstoring.getTimeInSeconds();
	  }
	});

      // Do the color-contraction
      auto moms = SB::getMoms(params.param.decay_dir, phases, SB::none, SB::none, tfrom, tsize);
      SB::doMomDisp_colorContractions(
	u_smr, source_colorvec, moms, tfrom, displacement_list, params.param.use_derivP, call,
	params.param.max_tslices_in_contraction, params.param.max_moms_in_contraction, SB::none,
	SB::OnDefaultDevice, SB::OnMaster);

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
  }   // namespace InlineBaryonMatElemColorVecSuperbEnv

  /*! @} */ // end of group hadron

} // namespace Chroma

#endif // BUILD_SB
