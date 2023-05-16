/*! \file
 * \brief Inline measurement of meson operators via colorstd::vector matrix elements
 */

#include "meas/inline/hadron/inline_meson_matelem_colorvec_superb_w.h"
#include "handle.h"
#include "meas/glue/mesplq.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/make_xml_file.h"
#include "meas/smear/displace.h"
#include "meas/smear/link_smearing_aggregate.h"
#include "meas/smear/link_smearing_factory.h"
#include "util/ferm/key_val_db.h"
#include "util/ferm/subset_vectors.h"
#include "util/ferm/superb_contractions.h"
#include "util/info/proginfo.h"

#include "meas/inline/io/named_objmap.h"

#define COLORVEC_MATELEM_TYPE_ZERO 0
#define COLORVEC_MATELEM_TYPE_ONE 1
#define COLORVEC_MATELEM_TYPE_MONE -1
#define COLORVEC_MATELEM_TYPE_GENERIC 10
#define COLORVEC_MATELEM_TYPE_DERIV 11

#ifdef BUILD_SB
namespace Chroma
{
  /*!
   * \ingroup hadron
   *
   * @{
   */
  namespace InlineMesonMatElemColorVecSuperbEnv
  {
    // Reader for input parameters
    void read(XMLReader& xml, const std::string& path,
	      InlineMesonMatElemColorVecSuperbEnv::Params::Param_t& param)
    {
      XMLReader paramtop(xml, path);

      int version;
      read(paramtop, "version", version);

      switch (version)
      {
      case 1:
	/**************************************************************************/
	param.mom2_min = 0;
	param.mom_list.resize(0);
	break;

      case 2:
	/**************************************************************************/
	read(paramtop, "mom2_min", param.mom2_min);
	param.mom_list.resize(0);
	break;

      case 3:
      case 4: // For compatibility with harom's task
	/**************************************************************************/
	read(paramtop, "mom2_min", param.mom2_min);
	read(paramtop, "mom_list", param.mom_list);
	break;

      default:
	/**************************************************************************/

	QDPIO::cerr << "Input parameter version " << version << " unsupported." << std::endl;
	QDP_abort(1);
      }

      param.use_derivP = true;
      if (paramtop.count("use_derivP") > 0)
	read(paramtop, "use_derivP", param.use_derivP);
      read(paramtop, "mom2_max", param.mom2_max);
      read(paramtop, "displacement_length", param.displacement_length);
      read(paramtop, "displacement_list", param.displacement_list);
      read(paramtop, "num_vecs", param.num_vecs);
      read(paramtop, "decay_dir", param.decay_dir);

      param.t_source = 0;
      if (paramtop.count("t_source") > 0)
	read(paramtop, "t_source", param.t_source);

      param.Nt_forward = Layout::lattSize()[param.decay_dir];
      if (paramtop.count("Nt_forward") > 0)
	read(paramtop, "Nt_forward", param.Nt_forward);

      param.max_tslices_in_contraction = Layout::logicalSize()[param.decay_dir];
      if (paramtop.count("max_tslices_in_contraction") == 1)
      {
	read(paramtop, "max_tslices_in_contraction", param.max_tslices_in_contraction);
      }

      param.max_moms_in_contraction = 1;
      if (paramtop.count("max_moms_in_contraction") == 1)
      {
	read(paramtop, "max_moms_in_contraction", param.max_moms_in_contraction);
      }

      if (paramtop.count("phase") == 1)
      {
	read(paramtop, "phase", param.quarkPhase);
	if (paramtop.count("quarkPhase") == 1 || paramtop.count("quarkPhase") == 1)
	{
	  QDPIO::cerr << "Error: please don't give the tag `phase' and either `quarkPhase' or "
			 "`aQuarkPhase'"
		      << std::endl;
	  QDP_abort(1);
	}
      }
      else if (paramtop.count("quarkPhase") == 1)
      {
	read(paramtop, "quarkPhase", param.quarkPhase);
      }
      else if (paramtop.count("aQuarkPhase") == 1)
      {
	QDPIO::cerr << "Label `aQuarkPhase' without the label `quarkPhase'" << std::endl;
	QDP_abort(1);
      }
      else
      {
	param.quarkPhase.resize(Nd - 1);
      }

      if (paramtop.count("aQuarkPhase") == 1)
      {
	read(paramtop, "aQuarkPhase", param.aQuarkPhase);
      }
      else
      {
	for (float i : param.quarkPhase)
	  param.aQuarkPhase.push_back(-i);
      }

      param.link_smearing = readXMLGroup(paramtop, "LinkSmearing", "LinkSmearingType");
    }

    // Writer for input parameters
    void write(XMLWriter& xml, const std::string& path,
	       const InlineMesonMatElemColorVecSuperbEnv::Params::Param_t& param)
    {
      push(xml, path);

      int version = 3;

      write(xml, "version", version);
      write(xml, "use_derivP", param.use_derivP);
      write(xml, "mom2_min", param.mom2_min);
      write(xml, "mom2_max", param.mom2_max);
      write(xml, "mom_list", param.mom_list);
      write(xml, "displacement_length", param.displacement_length);
      write(xml, "displacement_list", param.displacement_list);
      write(xml, "num_vecs", param.num_vecs);
      write(xml, "decay_dir", param.decay_dir);
      write(xml, "t_source", param.t_source);
      write(xml, "Nt_forward", param.Nt_forward);
      write(xml, "max_tslices_in_contraction", param.max_tslices_in_contraction);
      write(xml, "max_moms_in_contraction", param.max_moms_in_contraction);
      write(xml, "quarkPhase", SB::tomulti1d(param.quarkPhase));
      write(xml, "aQuarkPhase", SB::tomulti1d(param.aQuarkPhase));
      xml << param.link_smearing.xml;

      pop(xml);
    }

    //! Read named objects
    void read(XMLReader& xml, const std::string& path,
	      InlineMesonMatElemColorVecSuperbEnv::Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_id", input.gauge_id);
      read(inputtop, "colorvec_files", input.colorvec_files);
      read(inputtop, "meson_op_file", input.meson_op_file);
    }

    //! Write named objects
    void write(XMLWriter& xml, const std::string& path,
	       const InlineMesonMatElemColorVecSuperbEnv::Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "gauge_id", input.gauge_id);
      write(xml, "colorvec_files", input.colorvec_files);
      write(xml, "meson_op_file", input.meson_op_file);

      pop(xml);
    }

    // Writer for input parameters
    void write(XMLWriter& xml, const std::string& path,
	       const InlineMesonMatElemColorVecSuperbEnv::Params& param)
    {
      param.writeXML(xml, path);
    }
  }

  namespace InlineMesonMatElemColorVecSuperbEnv
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

    const std::string name = "MESON_MATELEM_COLORVEC_SUPERB";

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
      param.mom_list.resize(0);
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
    //! Meson operator
    struct KeyMesonElementalOperator_t {
      int t_slice;		 /*!< Meson operator time slice */
      multi1d<int> displacement; /*!< Displacement dirs of right colorstd::vector */
      multi1d<int> mom;		 /*!< D-1 momentum of this operator */
    };

    //! Meson operator, colorstd::vector source and sink with momentum projection
    struct ValMesonElementalOperator_t : public SB::Tensor<2, SB::ComplexD> {
      int type_of_data; /*!< Flag indicating type of data (maybe trivial) */
      ValMesonElementalOperator_t(int n = 0, int type_of_data = COLORVEC_MATELEM_TYPE_GENERIC)
	: SB::Tensor<2, SB::ComplexD>("ij", {n, n}, SB::OnHost, SB::Local),
	  type_of_data(type_of_data)
      {
      }
    };

    //----------------------------------------------------------------------------
    //! KeyMesonElementalOperator reader
    void read(BinaryReader& bin, KeyMesonElementalOperator_t& param)
    {
      read(bin, param.t_slice);
      read(bin, param.displacement);
      read(bin, param.mom);
    }

    //! MesonElementalOperator write
    void write(BinaryWriter& bin, const KeyMesonElementalOperator_t& param)
    {
      write(bin, param.t_slice);
      write(bin, param.displacement);
      write(bin, param.mom);
    }

    //! MesonElementalOperator reader
    void read(XMLReader& xml, const std::string& path, KeyMesonElementalOperator_t& param)
    {
      XMLReader paramtop(xml, path);

      read(paramtop, "t_slice", param.t_slice);
      read(paramtop, "displacement", param.displacement);
      read(paramtop, "mom", param.mom);
    }

    //! MesonElementalOperator writer
    void write(XMLWriter& xml, const std::string& path, const KeyMesonElementalOperator_t& param)
    {
      push(xml, path);

      write(xml, "t_slice", param.t_slice);
      write(xml, "displacement", param.displacement);
      write(xml, "mom", param.mom);

      pop(xml);
    }

    //----------------------------------------------------------------------------
    //! MesonElementalOperator reader
    void read(BinaryReader& bin, ValMesonElementalOperator_t& param)
    {
      int type_of_data;
      read(bin, type_of_data);

      int n;
      read(bin, n); // the size is always written, even if 0
      param = ValMesonElementalOperator_t(n, type_of_data);
      SB::Tensor<2, SB::ComplexD>& t = param;
      read(bin, t);
    }

    //! MesonElementalOperator write
    void write(BinaryWriter& bin, const ValMesonElementalOperator_t& param)
    {
      int type_of_data = param.type_of_data;
      write(bin, param.type_of_data);

      auto kvdim = param.kvdim();
      assert(kvdim['i'] == kvdim['j']);
      int n = kvdim['i']; // all sizes the same
      write(bin, n);
      SB::Tensor<2, SB::ComplexD> t = param.reorder("ij");
      write(bin, t);
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
    std::vector<std::vector<int>>
    normalizeDisplacements(const multi1d<multi1d<int>>& displacement_list)
    {
      std::vector<std::vector<int>> r(displacement_list.size());
      for (int i = 0; i < displacement_list.size(); ++i)
	r[i] = normDisp(displacement_list[i]);
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

	push(xml_out, "MesonMatElemColorVec");
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

      StopWatch swiss;
      swiss.reset();
      swiss.start();

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

      push(xml_out, "MesonMatElemColorVec");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << name << ": Meson color-std::vector matrix element" << std::endl;

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
      else
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

      // Record the smeared observables
      MesPlq(xml_out, "Smeared_Observables", u_smr);

      //
      // Parse the phase
      //
      if (params.param.quarkPhase.size() != Nd - 1 || params.param.aQuarkPhase.size() != Nd - 1)
      {
	QDPIO::cerr << "`phase', `quarkPhase', and `aQuarkPhase' tags should have " << Nd - 1
		    << " components" << std::endl;
	QDP_abort(1);
      }
      SB::Coor<Nd - 1> leftphase, rightphase;
      for (int i = 0; i < Nd - 1; ++i)
      {
	if (std::fabs((int)params.param.quarkPhase[i] - params.param.quarkPhase[i]) > 0 ||
	    std::fabs((int)params.param.aQuarkPhase[i] - params.param.aQuarkPhase[i]) > 0)
	  std::runtime_error("phase', `quarkPhase', and `aQuarkPhase' should be integer");
	rightphase[i] = params.param.quarkPhase[i];
	leftphase[i] = params.param.aQuarkPhase[i];
      }

      //
      // DB storage
      // NOTE: Only the master node opens the storage and writes on it
      //
      std::vector<LocalBinaryStoreDB<LocalSerialDBKey<KeyMesonElementalOperator_t>,
				     LocalSerialDBData<ValMesonElementalOperator_t>>>
	qdp_db;

      // This function opens the output file
      // NOTE: Only called by the master node
      std::function<void()> open_db = [&]() {
	// If the qdp_db is already opened, do nothing
	if (qdp_db.size() > 0)
	  return;
	qdp_db.resize(1);

	// Open the file, and write the meta-data and the binary for this operator
	if (!qdp_db[0].fileExists(params.named_obj.meson_op_file))
	{
	  XMLBufferWriter file_xml;

	  push(file_xml, "DBMetaData");
	  write(file_xml, "id", std::string("mesonElemOp"));
	  write(file_xml, "lattSize", QDP::Layout::lattSize());
	  //	write(file_xml, "blockSize", params.param.block_size);
	  write(file_xml, "decay_dir", params.param.decay_dir);
	  proginfo(file_xml); // Print out basic program info
	  write(file_xml, "Params", params.param);
	  write(file_xml, "Op_Info", params.param.displacement_list);
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

	  qdp_db[0].open(params.named_obj.meson_op_file, O_RDWR | O_CREAT, 0664);

	  qdp_db[0].insertUserdata(file_str);
	}
	else
	{
	  qdp_db[0].open(params.named_obj.meson_op_file, O_RDWR, 0664);
	}
      };

      // Compute the interval of t points to compute
      int tfrom = params.param.t_source;   // First t-slice to compute
      int tsize = params.param.Nt_forward; // Number of t-slices to compute
      const int Nt = Layout::lattSize()[params.param.decay_dir];

      // Make sure displacements are something sensible and transform to std objects
      std::vector<std::vector<int>> displacement_list =
	normalizeDisplacements(params.param.displacement_list);

      //
      // Meson operators
      //
      // The creation and annihilation operators are the same without the
      // spin matrices.
      //
      QDPIO::cout << "Building meson operators" << std::endl;

      push(xml_out, "ElementalOps");

      double time_storing = 0; // total time in writing elementals

      // Iterate over time-slices
      for (int tfrom0 = 0, this_tsize = std::min(tsize, params.param.max_tslices_in_contraction);
	   tfrom0 < tsize; tfrom0 += this_tsize,
	       this_tsize = std::min(params.param.max_tslices_in_contraction, tsize - tfrom0))
      {
	int this_tfrom = (tfrom + tfrom0) % Nt;

	// Get num_vecs colorvecs on time-slice t_source
	SB::Tensor<Nd + 3, SB::Complex> source_colorvec =
	  SB::getColorvecs<SB::Complex>(colorvecsSto, u, params.param.decay_dir, this_tfrom,
					this_tsize, params.param.num_vecs, SB::none);

	// Call for storing the baryons
	SB::ContractionFn<SB::Complex> call(
	  [&](SB::Tensor<4, SB::Complex> tensor, int disp, int first_tslice, int first_mom) {
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

	      KeyMesonElementalOperator_t key;
	      ValMesonElementalOperator_t val(
		params.param.num_vecs, params.param.use_derivP ? COLORVEC_MATELEM_TYPE_DERIV
							       : COLORVEC_MATELEM_TYPE_GENERIC);

	      // The keys for the displacements for this particular elemental operator
	      // Invert the time - make it an independent key
	      for (int t = 0, numt = tensor.kvdim()['t']; t < numt; ++t)
	      {
		for (int m = 0, numm = tensor.kvdim()['m']; m < numm; ++m)
		{
		  key.t_slice = first_tslice + t;
		  key.mom = SB::tomulti1d(mom_list[first_mom + m]);
		  key.displacement =
		    SB::tomulti1d(displacement_list[disp]); // only right colorstd::vector
		  tensor.kvslice_from_size({{'t', t}, {'m', m}}, {{'t', 1}, {'m', 1}}).copyTo(val);
		  qdp_db[0].insert(key, val);
		}
	      }

	      tstoring.stop();
	      time_storing += tstoring.getTimeInSeconds();
	    }
	  });

	// Do the contractions
	SB::doMomDisp_contractions(u_smr, source_colorvec, leftphase, rightphase, mom_list,
				   this_tfrom, displacement_list, params.param.use_derivP, call,
				   SB::none, SB::OnDefaultDevice, SB::OnMaster,
				   0 /* max_tslices_in_contraction==0 means to do all */,
				   params.param.max_moms_in_contraction);
      }

      // Close db
      for (auto& db : qdp_db)
	db.close();

      // Close colorvecs storage
      SB::closeColorvecStorage(colorvecsSto);

      swiss.stop();

      QDPIO::cout << "All meson operators computed in time= "
		  << swiss.getTimeInSeconds() - time_storing << " secs, writing time "
		  << time_storing << " secs " << std::endl;

      pop(xml_out); // ElementalOps

      // Close the namelist output file XMLDAT
      pop(xml_out); // MesonMatElemColorVector

      QDPIO::cout << name << ": total time = " << swiss.getTimeInSeconds() << " secs" << std::endl;

      QDPIO::cout << name << ": ran successfully" << std::endl;

      END_CODE();
    } // func
  }   // namespace InlineMesonMatElemColorVecSuperbEnv

  /*! @} */ // end of group hadron

} // namespace Chroma

#endif // BUILD_SB
