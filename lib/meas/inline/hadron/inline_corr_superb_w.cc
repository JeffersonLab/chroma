/*! \file
 * \brief Inline measurement of baryon operators via colorstd::vector matrix elements
 */

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
#include "util/ferm/superb_contractions.h"
#include "util/info/proginfo.h"

#ifdef BUILD_SB
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
	      InlineCorrSuperbEnv::Params::Param_t& param)
    {
      XMLReader paramtop(xml, path);

      param.num_vecs = 0;
      read(paramtop, "num_vecs", param.num_vecs);

      read(paramtop, "flavor_to_mass", param.flavor_to_mass);

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

      param.link_smearing = readXMLGroup(paramtop, "LinkSmearing", "LinkSmearingType");
    }

    // Writer for input parameters
    void write(XMLWriter& xml, const std::string& path,
	       const InlineCorrSuperbEnv::Params::Param_t& param)
    {
      push(xml, path);

      write(xml, "num_vecs", param.num_vecs);
      write(xml, "flavor_to_mass", param.flavor_to_mass);
      write(xml, "t_origin", param.t_origin);
      write(xml, "mesons", param.meson_files);
      write(xml, "baryons", param.baryon_files);
      write(xml, "props", param.prop_files);
      write(xml, "genprops", param.genprop_files);
      write(xml, "Nt_forward", param.Nt_forward);
      write(xml, "ensemble", param.ensemble);
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

      const bool zeroUnsmearedGraphsP = true;
      const auto& corr = Hadron::evaluate_graphs_with_superb(
	corr_graph, zeroUnsmearedGraphsP, storage_prop, storage_baryon, storage_meson,
	storage_genprop, flavor_to_mass, nev, params.param.t_origin, params.param.Nt_forward);

      std::cout << "Storing the correlation functions" << std::endl;
      const int decay_dir = 3;
      Hadron::writeCorrMap(corr, params.param.ensemble, corr_graph.layout.latt_size, decay_dir,
			   params.param.t_origin, params.named_obj.corr_file);

      // Close colorvecs storage
      SB::closeColorvecStorage(colorvecsSto);

      pop(xml_out); // ElementalOps

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

#endif // BUILD_SB
