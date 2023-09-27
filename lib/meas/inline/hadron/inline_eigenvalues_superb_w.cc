/*! \file
 * \brief Compute the smallest eigenvalues of D^\dagger*D
 *
 * Compute an approximation of the smallest eigenvalues/vectors of D^\dagger * D with
 * an iterative eigensolver on \gamma_5 D^{-1}, where D^{-1} is approximate with
 * a linear solver
 */

#include "qdp.h"
#include "fermact.h"
#include "meas/inline/hadron/inline_eigenvalues_superb_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "qdp_map_obj.h"
#include "qdp_map_obj_disk.h"
#include "qdp_map_obj_disk_multiple.h"
#include "qdp_map_obj_memory.h"
#include "qdp_disk_map_slice.h"
#include "util/ferm/subset_vectors.h"
#include "util/ferm/key_val_db.h"
#include "util/ferm/transf.h"
#include "util/ferm/spin_rep.h"
#include "util/ferm/diractodr.h"
#include "util/ferm/twoquark_contract_ops.h"
#include "util/ferm/superb_contractions.h"
#include "util/ferm/mgproton.h"
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
  namespace InlineEigenvaluesSuperbEnv 
  {
    //! Propagator input
    void read(XMLReader& xml, const std::string& path, Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_id", input.gauge_id);
      read(inputtop, "eigs_file", input.eigs_file);
    }

    //! Propagator output
    void write(XMLWriter& xml, const std::string& path, const Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "gauge_id", input.gauge_id);
      write(xml, "eigs_file", input.eigs_file);

      pop(xml);
    }


    //! Propagator input
    void read(XMLReader& xml, const std::string& path, Params::Param_t::Contract_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "num_vecs", input.num_vecs);
      read(inputtop, "tolerance", input.tol);
      read(inputtop, "mass_label", input.mass_label);

      input.max_rhs = 8;
      if( inputtop.count("max_rhs") == 1 ) {
        read(inputtop, "max_rhs", input.max_rhs);
      }
    }

    //! Propagator output
    void write(XMLWriter& xml, const std::string& path, const Params::Param_t::Contract_t& input)
    {
      push(xml, path);

      write(xml, "num_vecs", input.num_vecs);
      write(xml, "tolerance", input.tol);
      write(xml, "mass_label", input.mass_label);
      write(xml, "max_rhs", input.max_rhs);

      pop(xml);
    }


    //! Propagator input
    void read(XMLReader& xml, const std::string& path, Params::Param_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "Propagator", input.prop);

      input.eigensolver = "<eigensolver></eigensolver>";
      if (inputtop.count("eigensolver") == 1)
      {
	XMLReader xml_tmp(inputtop, "eigensolver");
	std::ostringstream os;
	xml_tmp.print(os);
	input.eigensolver = os.str();
      }

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

    //----------------------------------------------------------------------------
    //! Unsmeared meson operator
    struct KeyEigenpair_t
    {
      int idx;			       /*!< eigenpair index */
      std::string mass_label;	       /*!< Some kind of mass label */
    };

    //! Eigenvalue and eigenvector
    struct ValEigenpair_t : public SB::Tensor<6, SB::ComplexD>
    {
      double value;
      ValEigenpair_t()
	: SB::Tensor<6, SB::ComplexD>("csxyzt",
				      SB::latticeSize<6>("csxyzt", {{'x', Layout::lattSize()[0]}}),
				      SB::OnHost, SB::Local)
      {
      }
    };


    //----------------------------------------------------------------------------
    //! KeyEigenpair_t reader
    void read(BinaryReader& bin, KeyEigenpair_t& param)
    {
      read(bin, param.idx);
      readDesc(bin, param.mass_label);
    }

    //! KeyEigenpair_t write
    void write(BinaryWriter& bin, const KeyEigenpair_t& param)
    {
      write(bin, param.idx);
      writeDesc(bin, param.mass_label);
    }

    //----------------------------------------------------------------------------
    //! ValEigenpair_t reader
    void read(BinaryReader& bin, ValEigenpair_t& param)
    {
	double value;
	read(bin, value);
	param = ValEigenpair_t();
	SB::Tensor<6, SB::ComplexD> &t = param;
	read(bin, t);
    }

    //! ValEigenpair_t write
    void write(BinaryWriter& bin, const ValEigenpair_t& param)
    {
	write(bin, param.value);
	SB::Tensor<6, SB::ComplexD> t = param.reorder("csxyzt");
	write(bin, t);
    }


  } // namespace InlinePropDistillationSuperbEnv 


  //----------------------------------------------------------------------------
  namespace InlineEigenvaluesSuperbEnv
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
      
    const std::string name = "EIGENVALUES_SUPERB";

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

	push(xml_out, "Eigenvalues");
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

      push(xml_out, "Eigenvpairs");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << name << ": eigenpairs calculation" << std::endl;

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
      const int decay_dir = 3;
      const int Lt        = Layout::lattSize()[decay_dir];

      //
      // DB storage
      //
      std::vector<
	LocalBinaryStoreDB<LocalSerialDBKey<KeyEigenpair_t>, LocalSerialDBData<ValEigenpair_t>>>
	qdp_db;

      // Open the file, and write the meta-data and the binary for this operator
      auto open_db = [&]() {
	if (qdp_db.size() > 0)
	  return;
	XMLBufferWriter file_xml;
	push(file_xml, "DBMetaData");
	write(file_xml, "id", std::string("eigenpairsOp"));
	write(file_xml, "lattSize", QDP::Layout::lattSize());
	proginfo(file_xml); // Print out basic program info
	write(file_xml, "Params", params.param);
	write(file_xml, "Config_info", gauge_xml);
	pop(file_xml);

	std::string file_str(file_xml.str());
	qdp_db.resize(1);
	qdp_db[0].setMaxUserInfoLen(file_str.size());

	qdp_db[0].open(params.named_obj.eigs_file, O_RDWR | O_CREAT, 0664);

	qdp_db[0].insertUserdata(file_str);
      };

      //
      // Try the factories
      //
      try
      {
	StopWatch swatch;
	swatch.reset();
	QDPIO::cout << "Try the various factories" << std::endl;

	// Initialize fermion action and create the solver
	SB::ChimeraSolver PP{params.param.prop.fermact, params.param.prop.invParam, u};

	// Prepare eigensolver
	std::shared_ptr<SB::Options> ops =
	  SB::getOptionsFromXML(SB::broadcast(params.param.eigensolver));
	auto eigensolver = SB::detail::getInexactEigensolverGD(
	  SB::getOperator(PP, params.param.contract.max_rhs), ops->getValue("eigensolver"));

	// Run
	auto values_vectors = eigensolver(params.param.contract.num_vecs, params.param.contract.tol);
	auto values = std::get<0>(values_vectors);
	auto vectors = std::get<1>(values_vectors);
	assert(values.size() == (std::size_t)vectors.kvdim().at('n'));

	swatch.stop();
	QDPIO::cout << "Eigenpairs computed: time= " << swatch.getTimeInSeconds() << " secs"
		    << std::endl;

	// Store the eigenpairs
	swatch.reset();
	LocalSerialDBKey<KeyEigenpair_t> key;
	key.key().mass_label = params.param.contract.mass_label;
	for (int idx = 0; idx < values.size(); ++idx)
	{
	  key.key().idx = idx;
	  auto v =
	    SB::detail::toNaturalOrdering(vectors.kvslice_from_size({{'n', idx}}, {{'n', 1}}))
	      .make_sure(SB::none, SB::OnHost, SB::OnMaster)
	      .getLocal();
	  if (v)
	  {
	    open_db();
	    LocalSerialDBData<ValEigenpair_t> val;
	    val.data().value = values[idx];
	    v.copyTo(val.data());
	    qdp_db[0].insert(key, val);
	  }
	}

	swatch.stop();
	QDPIO::cout << "Eigenpairs stored: time= " << swatch.getTimeInSeconds() << " secs"
		    << std::endl;
      }
      catch (const std::exception& e) 
      {
	QDP_error_exit("%s: caught exception: %s\n", name.c_str(), e.what());
      }

      pop(xml_out);

      for (auto& db : qdp_db)
	db.close();

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
