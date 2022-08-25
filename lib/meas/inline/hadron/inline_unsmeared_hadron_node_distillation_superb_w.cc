/*! \file
 * \brief Inline measurement that construct unsmeared hadron nodes using distillation
 */

#include "meas/inline/hadron/inline_unsmeared_hadron_node_distillation_superb_w.h"

#include "qdp_map_obj_memory.h"
#include "qdp_map_obj_disk.h"
#include "qdp_map_obj_disk_multiple.h"
#include "qdp_disk_map_slice.h"
#include "handle.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"
#include "io/xml_group_reader.h"
#include "meas/glue/mesplq.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/make_xml_file.h"
#include "meas/smear/link_smearing_aggregate.h"
#include "meas/smear/link_smearing_factory.h"
#include "util/ferm/disp_soln_cache.h"
#include "util/ferm/key_timeslice_colorvec.h"
#include "util/ferm/key_val_db.h"
#include "util/ferm/subset_vectors.h"
#include "util/ferm/superb_contractions.h"
#include "util/ferm/mgproton.h"
#include "util/ferm/transf.h"
#include "util/ft/sftmom.h"
#include "util/ft/time_slice_set.h"
#include "util/info/proginfo.h"

#include "meas/inline/io/named_objmap.h"
#include <set>

#ifdef BUILD_SB

namespace Chroma 
{ 
  /*!
   * \ingroup hadron
   *
   * @{
   */
  namespace InlineUnsmearedHadronNodeDistillationSuperbEnv 
  { 
    //! Propagator input
    void read(XMLReader& xml, const std::string& path, Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_id", input.gauge_id);
      read(inputtop, "colorvec_files", input.colorvec_files);
      read(inputtop, "dist_op_file", input.dist_op_file);
    }

    //! Propagator output
    void write(XMLWriter& xml, const std::string& path, const Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "gauge_id", input.gauge_id);
      write(xml, "colorvec_files", input.colorvec_files);
      write(xml, "dist_op_file", input.dist_op_file);

      pop(xml);
    }


    //! Propagator input
    void read(XMLReader& xml, const std::string& path, Params::Param_t::DispGammaMom_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "displacement", input.displacement);
      read(inputtop, "gamma", input.gamma);
      read(inputtop, "mom", input.mom);
    }

    //! Propagator output
    void write(XMLWriter& xml, const std::string& path, const Params::Param_t::DispGammaMom_t& input)
    {
      push(xml, path);

      write(xml, "displacement", input.displacement);
      write(xml, "gamma", input.gamma);
      write(xml, "mom", input.mom);

      pop(xml);
    }


    //! Propagator input
    void read(XMLReader& xml, const std::string& path, Params::Param_t::KeySolnProp_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "cacheP", input.cacheP);
      read(inputtop, "num_vecs", input.num_vecs);
      read(inputtop, "t_source", input.t_source);
    }

    //! Propagator output
    void write(XMLWriter& xml, const std::string& path, const Params::Param_t::KeySolnProp_t& input)
    {
      push(xml, path);

      write(xml, "cacheP", input.cacheP);
      write(xml, "num_vecs", input.num_vecs);
      write(xml, "t_source", input.t_source);

      pop(xml);
    }


    //! Propagator input
    void read(XMLReader& xml, const std::string& path, Params::Param_t::SinkSource_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "t_sink", input.t_sink);
      read(inputtop, "t_source", input.t_source);
      read(inputtop, "Nt_forward", input.Nt_forward);
      read(inputtop, "Nt_backward", input.Nt_backward);
    }

    //! Propagator output
    void write(XMLWriter& xml, const std::string& path, const Params::Param_t::SinkSource_t& input)
    {
      push(xml, path);

      write(xml, "t_sink", input.t_sink);
      write(xml, "t_source", input.t_source);
      write(xml, "Nt_forward", input.Nt_forward);
      write(xml, "Nt_backward", input.Nt_backward);

      pop(xml);
    }

    //! Propagator input
    void read(XMLReader& xml, const std::string& path, Params::Param_t::Contract_t& input)
    {
      XMLReader inputtop(xml, path);

      input.alt_t_start = std::numeric_limits<int>::max();
      if (inputtop.count("t_start") == 1) {
        read(inputtop, "t_start", input.alt_t_start);
      }

      input.alt_Nt_forward = std::numeric_limits<int>::max();
      if (inputtop.count("Nt_forward") == 1) {
        read(inputtop, "Nt_forward", input.alt_Nt_forward);
      }

      input.alt_num_vecs = 0;
      if (inputtop.count("num_vecs") == 1) {
        read(inputtop, "num_vecs", input.alt_num_vecs);
      }

      read(inputtop, "use_derivP", input.use_derivP);
      read(inputtop, "decay_dir", input.decay_dir);
      read(inputtop, "displacement_length", input.displacement_length);
      read(inputtop, "mass_label", input.mass_label);
      read(inputtop, "num_tries", input.num_tries);

      input.max_rhs = 8;
      if( inputtop.count("max_rhs") == 1 ) {
        read(inputtop, "max_rhs", input.max_rhs);
      }

      input.max_tslices_in_contraction = 0;
      if( inputtop.count("max_tslices_in_contraction") == 1 ) {
        read(inputtop, "max_tslices_in_contraction", input.max_tslices_in_contraction);
      }

      input.max_moms_in_contraction = 0;
      if( inputtop.count("max_moms_in_contraction") == 1 ) {
        read(inputtop, "max_moms_in_contraction", input.max_moms_in_contraction);
      }

      input.use_device_for_contractions = true;
      if( inputtop.count("use_device_for_contractions") == 1 ) {
        read(inputtop, "use_device_for_contractions", input.use_device_for_contractions);
      }

      input.use_genprop4_format = false;
      if( inputtop.count("use_genprop4_format") == 1 ) {
        read(inputtop, "use_genprop4_format", input.use_genprop4_format);
      }

      input.use_genprop5_format = false;
      if( inputtop.count("use_genprop5_format") == 1 ) {
        read(inputtop, "use_genprop5_format", input.use_genprop5_format);
      }

      input.use_multiple_writers = false;
      if( inputtop.count("use_multiple_writers") == 1 ) {
        read(inputtop, "use_multiple_writers", input.use_multiple_writers);
      }

      input.phase.resize(Nd - 1);
      for (int i = 0; i < Nd - 1; ++i)
	input.phase[i] = 0;
      if( inputtop.count("phase") == 1 ) {
        read(inputtop, "phase", input.phase);
      }
     }

    //! Propagator output
    void write(XMLWriter& xml, const std::string& path, const Params::Param_t::Contract_t& input)
    {
      push(xml, path);

      write(xml, "t_start", input.alt_t_start);
      write(xml, "Nt_forward", input.alt_Nt_forward);
      write(xml, "num_vecs", input.alt_num_vecs);
      write(xml, "use_derivP", input.use_derivP);
      write(xml, "decay_dir", input.decay_dir);
      write(xml, "displacement_length", input.displacement_length);
      write(xml, "mass_label", input.mass_label);
      write(xml, "num_tries", input.num_tries);
      write(xml, "max_rhs", input.max_rhs);
      write(xml, "max_tslices_in_contraction", input.max_tslices_in_contraction);
      write(xml, "max_moms_in_contraction", input.max_moms_in_contraction);
      write(xml, "use_genprop4_format", input.use_genprop4_format);
      write(xml, "use_genprop5_format", input.use_genprop5_format);
      write(xml, "use_device_for_contractions", input.use_device_for_contractions);
      write(xml, "use_multiple_writers", input.use_multiple_writers);
      write(xml, "phase", input.phase);

      pop(xml);
    }


    //! Propagator input
    void read(XMLReader& xml, const std::string& path, Params::Param_t& input)
    {
      XMLReader inputtop(xml, path);

      if (inputtop.count("DispGammaMomList") == 1)
        read(inputtop, "DispGammaMomList", input.disp_gamma_mom_list);
      if (inputtop.count("Propagator") == 1)
        read(inputtop, "Propagator", input.prop);
      if (inputtop.count("PropSources") == 1)
        read(inputtop, "PropSources", input.prop_sources);
      if (inputtop.count("SinkSourcePairs") == 1)
        read(inputtop, "SinkSourcePairs", input.sink_source_pairs);

      if (inputtop.count("Displacements") == 1)
        read(inputtop, "Displacements", input.alt_displacements);
      if (inputtop.count("Moms") == 1)
        read(inputtop, "Moms", input.alt_moms);
      if (inputtop.count("SinkSources") == 1)
        read(inputtop, "SinkSources", input.alt_sink_sources);

      read(inputtop, "Contractions", input.contract);

      input.link_smearing  = readXMLGroup(inputtop, "LinkSmearing", "LinkSmearingType");
    }

    //! Propagator output
    void write(XMLWriter& xml, const std::string& path, const Params::Param_t& input)
    {
      push(xml, path);

      write(xml, "DispGammaMomList", input.disp_gamma_mom_list);
      write(xml, "Propagator", input.prop);
      write(xml, "PropSources", input.prop_sources);
      write(xml, "SinkSourcePairs", input.sink_source_pairs);
      write(xml, "Contractions", input.contract);

      write(xml, "Displacements", input.alt_displacements);
      write(xml, "Moms", input.alt_moms);
      write(xml, "SinkSources", input.alt_sink_sources);
      xml << input.link_smearing.xml;

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
  } // end namespace



  //-------------------------------------------------------------------------------
  namespace InlineUnsmearedHadronNodeDistillationSuperbEnv
  { 
    // Anonymous namespace for registration
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

    const std::string name = "UNSMEARED_HADRON_NODE_DISTILLATION_SUPERB";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
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

	  for(int i=1; i < d.size(); ++i)
	    os << " " << d[i];
	}

	return os;
      }

      StandardOutputStream& operator<<(StandardOutputStream& os, const std::vector<int>& d)
      {
	if (d.size() > 0)
	{
	  os << d[0];

	  for(int i=1; i < d.size(); ++i)
	    os << " " << d[i];
	}

	return os;
      }
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
    //! KeySolnProp reader
    void read(BinaryReader& bin, Params::Param_t::KeySolnProp_t& param)
    {
      read(bin, param.cacheP);
      read(bin, param.num_vecs);
      read(bin, param.t_source);
    }

    //! KeySolnProp write
    void write(BinaryWriter& bin, const Params::Param_t::KeySolnProp_t& param)
    {
      write(bin, param.cacheP);
      write(bin, param.num_vecs);
      write(bin, param.t_source);
    }

    //----------------------------------------------------------------------------
    //! Unsmeared meson operator
    struct KeyUnsmearedMesonElementalOperator_t
    {
      int                t_sink;       /*!< Time sink */
      int                t_slice;      /*!< Meson operator time slice */
      int                t_source;     /*!< Time source */
      int                colorvec_src; /*!< Colorvec source component */
      int                gamma;        /*!< The "gamma" in Gamma(gamma) */
      bool               derivP;       /*!< Uses derivatives */
      std::vector<int>   displacement; /*!< Displacement dirs of right colorstd::vector */
      multi1d<int>       mom;          /*!< D-1 momentum of this operator */
      std::string        mass;         /*!< Some kind of mass label */
    };

    //! Meson operator: spin source, spin sink, colorvector
    struct ValUnsmearedMesonElementalOperator_t : public SB::Tensor<3, SB::ComplexD>
    {
      ValUnsmearedMesonElementalOperator_t(int num_vecs = 0)
	: SB::Tensor<3, SB::ComplexD>("Nqs", {num_vecs, Ns, Ns}, SB::OnHost, SB::Local)
      {
      }
    };


    //----------------------------------------------------------------------------
    //! KeyUnsmearedMesonElementalOperator reader
    void read(BinaryReader& bin, KeyUnsmearedMesonElementalOperator_t& param)
    {
      read(bin, param.t_sink);
      read(bin, param.t_slice);
      read(bin, param.t_source);
      read(bin, param.colorvec_src);
      read(bin, param.gamma);
      read(bin, param.derivP);
      read(bin, param.displacement);
      read(bin, param.mom);
      readDesc(bin, param.mass);
    }

    //! UnsmearedMesonElementalOperator write
    void write(BinaryWriter& bin, const KeyUnsmearedMesonElementalOperator_t& param)
    {
      write(bin, param.t_sink);
      write(bin, param.t_slice);
      write(bin, param.t_source);
      write(bin, param.colorvec_src);
      write(bin, param.gamma);
      write(bin, param.derivP);
      write(bin, param.displacement);
      write(bin, param.mom);
      writeDesc(bin, param.mass);
    }

    //! UnsmearedMesonElementalOperator reader
    void read(XMLReader& xml, const std::string& path, KeyUnsmearedMesonElementalOperator_t& param)
    {
      XMLReader paramtop(xml, path);
    
      read(paramtop, "t_sink", param.t_sink);
      read(paramtop, "t_slice", param.t_slice);
      read(paramtop, "t_source", param.t_source);
      read(paramtop, "colorvec_src", param.colorvec_src);
      read(paramtop, "gamma", param.gamma);
      read(paramtop, "derivP", param.derivP);
      read(paramtop, "displacement", param.displacement);
      read(paramtop, "mom", param.mom);
      read(paramtop, "mass", param.mass);
    }

    //! UnsmearedMesonElementalOperator writer
    void write(XMLWriter& xml, const std::string& path, const KeyUnsmearedMesonElementalOperator_t& param)
    {
      push(xml, path);

      write(xml, "t_sink", param.t_sink);
      write(xml, "t_slice", param.t_slice);
      write(xml, "t_source", param.t_source);
      write(xml, "colorvec_src", param.colorvec_src);
      write(xml, "gamma", param.gamma);
      write(xml, "derivP", param.derivP);
      write(xml, "displacement", param.displacement);
      write(xml, "mom", param.mom);
      write(xml, "mass", param.mass);

      pop(xml);
    }


    //----------------------------------------------------------------------------
    //! UnsmearedMesonElementalOperator reader
    void read(BinaryReader& bin, ValUnsmearedMesonElementalOperator_t& param)
    {
	int num_vecs, q, s;
	read(bin, num_vecs);
	read(bin, q);
	read(bin, s);
	assert(s == Ns && q == Ns);
	param = ValUnsmearedMesonElementalOperator_t(num_vecs);
	SB::Tensor<3, SB::ComplexD> &t = param;
	read(bin, t);
    }

    //! UnsmearedMesonElementalOperator write
    void write(BinaryWriter& bin, const ValUnsmearedMesonElementalOperator_t& param)
    {
	auto kvdim = param.kvdim();
	int num_vecs = kvdim['N'], q = kvdim['q'], s = kvdim['s'];
	assert(s == Ns && q == Ns);
	write(bin, num_vecs);
	write(bin, q);
	write(bin, s);
	SB::Tensor<3, SB::ComplexD> t = param.reorder("Nqs");
	write(bin, t);
    }

    //----------------------------------------------------------------------------
    //! Generalized propagator operator (new UnsmearedMesonElementalOperator)
    struct KeyGenProp4ElementalOperator_t {
      int t_sink;		     /*!< Source time slice */
      int t_slice;		     /*!< Propagator time slice */
      int t_source;		     /*!< Source time slice */
      int g;			     /*!< DR gamma */
      std::vector<int> displacement; /*!< Displacement dirs of right colorstd::vector */
      multi1d<int> mom;		     /*!< D-1 momentum of this operator */
      std::string mass;		     /*!< A mass label */
    };

    //! Generalized propagator operator
    //! Meson operator: spin source, spin sink, colorvector source, colorvector sink
    struct ValGenProp4ElementalOperator_t : public SB::Tensor<4, SB::ComplexD>
    {
      ValGenProp4ElementalOperator_t(int num_vecs_source = 0, int num_vecs_sink = 0)
	: SB::Tensor<4, SB::ComplexD>("Nnqs", {num_vecs_sink, num_vecs_source, Ns, Ns}, SB::OnHost,
				      SB::Local)
      {
      }
    };

    //----------------------------------------------------------------------------
    StandardOutputStream& operator<<(StandardOutputStream& os,
				     const KeyGenProp4ElementalOperator_t& d)
    {
      os << "GenProp4:"
	 << " t_sink= " << d.t_sink << " t_slice= " << d.t_slice << " t_source= " << d.t_source
	 << " g= " << d.g << " displacement= " << d.displacement << " mom= " << d.mom
	 << " mass= " << d.mass;

      return os;
    }

    //----------------------------------------------------------------------------
    //! KeyGenProp4ElementalOperator reader
    void read(BinaryReader& bin, KeyGenProp4ElementalOperator_t& param)
    {
      read(bin, param.t_sink);
      read(bin, param.t_slice);
      read(bin, param.t_source);
      read(bin, param.g);
      read(bin, param.displacement);
      read(bin, param.mom);
      readDesc(bin, param.mass);
    }

    //! GenProp4ElementalOperator write
    void write(BinaryWriter& bin, const KeyGenProp4ElementalOperator_t& param)
    {
      write(bin, param.t_sink);
      write(bin, param.t_slice);
      write(bin, param.t_source);
      write(bin, param.g);
      write(bin, param.displacement);
      write(bin, param.mom);
      writeDesc(bin, param.mass);
    }

    //----------------------------------------------------------------------------
    //! PropElementalOperator reader
    void read(BinaryReader& bin, ValGenProp4ElementalOperator_t& param)
    {
      int num_vecs_source, num_vecs_sink, q, s;
      read(bin, num_vecs_sink);
      read(bin, num_vecs_source);
      read(bin, q);
      read(bin, s);
      assert(s == Ns && q == Ns);
      param = ValGenProp4ElementalOperator_t(num_vecs_source, num_vecs_sink);
      SB::Tensor<4, SB::ComplexD>& t = param;
      read(bin, t);
    }

    //! GenProp4ElementalOperator write
    void write(BinaryWriter& bin, const ValGenProp4ElementalOperator_t& param)
    {
      auto kvdim = param.kvdim();
      int num_vecs_source = kvdim['n'], num_vecs_sink = kvdim['N'], q = kvdim['q'], s = kvdim['s'];
      assert(s == Ns && q == Ns);
      write(bin, num_vecs_sink);
      write(bin, num_vecs_source);
      write(bin, q);
      write(bin, s);
      SB::Tensor<4, SB::ComplexD> t = param.reorder("Nnqs");
      write(bin, t);
    }

    //-------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------

    //----------------------------------------------------------------------------
    //! Normalize just one displacement array
    std::vector<int> normDisp(const std::vector<int>& orig)
    {
      START_CODE();

      std::vector<int> disp;
      std::vector<int> empty; 
      std::vector<int> no_disp(1); no_disp[0] = 0;

      // NOTE: a no-displacement is recorded as a zero-length array
      // Convert a length one array with no displacement into a no-displacement array
      if (orig.size() == 1)
      {
	if (orig == no_disp)
	  disp = empty;
	else
	  disp = orig;
      }
      else
      {
	disp = orig;
      }

      END_CODE();

      return disp;
    } // void normDisp



    //-------------------------------------------------------------------------------
    // Function call
    void 
    InlineMeas::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
      // If xml file not empty, then use alternate
      if (params.xml_file != "")
      {
	std::string xml_file = makeXMLFileName(params.xml_file, update_no);

	push(xml_out, "HadronNode");
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


    //-------------------------------------------------------------------------------
    // Function call
    void 
    InlineMeas::func(unsigned long update_no,
		     XMLWriter& xml_out) 
    {
      START_CODE();

      StopWatch snoop;
      snoop.reset();
      snoop.start();

      QDPIO::cout << name << ": construct unsmeared hadron nodes via distillation" << std::endl;
      
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

      push(xml_out, "UnsmearedHadronNode");
      write(xml_out, "update_no", update_no);

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

      if (params.param.contract.use_derivP && params.param.contract.use_genprop4_format)
	throw std::runtime_error("`use_genprop4_format` does not support `use_derivP` for now");

      //
      // Read in the source along with relevant information.
      // 

      SB::ColorvecsStorage colorvecsSto = SB::openColorvecStorage(params.named_obj.colorvec_files);

      //
      // Read through the momentum list and find all the unique phases
      //
      QDPIO::cout << "Parse momentum list" << std::endl;
      
      // Possible momenta, gammas, and displacements
      multi2d<int> moms;
      std::vector<int> gammas;
      std::vector<std::vector<int>> disps;

      {
        std::set<SB::Coor<Nd - 1>> moms_set;
        std::set<std::vector<int>> disps_set;
        std::set<int> gammas_set;

        for (const auto &ins : params.param.disp_gamma_mom_list) {
          SB::Coor<Nd - 1> c;
          for (unsigned int i = 0; i < ins.mom.size() && i < Nd - 1; ++i)
            c[i] = ins.mom[i];
          moms_set.insert(c);
          disps_set.insert(normDisp(ins.displacement));
          gammas_set.insert(ins.gamma);
        }

        if (params.param.disp_gamma_mom_list.size() == 0) {
          for (int g = 0; g < Nd * Nd; g++)
            gammas_set.insert(g);
        }

        for (const auto &it : params.param.alt_displacements)
          disps_set.insert(it);

        for (const auto &it : params.param.alt_moms)
        {
          SB::Coor<Nd - 1> c;
          for (unsigned int i = 0; i < it.size() && i < Nd - 1; ++i)
            c[i] = it[i];
          moms_set.insert(c);
        }

        if (moms_set.size() == 0 || disps_set.size() == 0 || gammas_set.size() == 0) {
          QDPIO::cerr << name
                      << ": warning - no moms nor displacement nor gammas; nothing to do"
                      << std::endl;
          QDP_abort(1);
        }

        int num_mom = moms_set.size();
        int mom_size = Nd - 1;
        QDPIO::cout << name << ": num_mom= " << num_mom
                    << "  mom_size= " << mom_size << std::endl;
        moms.resize(num_mom, mom_size);
        int i = 0;
        for (const auto &it : moms_set) {
          for (unsigned int j = 0; j < Nd - 1; ++j)
            moms[i][j] = it[j];
          i++;
        }

        disps.resize(disps_set.size());
        std::copy(disps_set.begin(), disps_set.end(), disps.begin());
        gammas.resize(gammas_set.size());
        std::copy(gammas_set.begin(), gammas_set.end(), gammas.begin());
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

      //
      // Initialize the slow Fourier transform phases
      //
      SftMom phases(moms, params.param.contract.decay_dir);

      //
      // Capture maximum number of vecs
      //
      int num_vecs = params.param.contract.alt_num_vecs;
      for (const auto& it : params.param.prop_sources)
	num_vecs = std::max(num_vecs, it.num_vecs);

      //
      // Stores the range of time-slices used for each sink/source
      //

      std::vector<bool> cache_tslice(Lt);
      for (const auto& it : params.param.alt_sink_sources)
      {
	for (const auto& snk : it.second)
	{
	  if (it.first < 0 || snk < 0)
	    throw std::runtime_error("Invalid source or sink on SinkSources");
	  Params::Param_t::SinkSource_t ss;
	  ss.t_sink = snk % Lt;
	  ss.t_source = it.first % Lt;
	  ss.Nt_backward = it.first - params.param.contract.alt_t_start;
	  ss.Nt_forward = params.param.contract.alt_Nt_forward - ss.Nt_backward;
	  params.param.sink_source_pairs.push_back(ss);
	  cache_tslice[ss.t_source] = cache_tslice[ss.t_sink] = true;
	}
      }

      struct FromSize {
	int from;
	int size;
      };
      std::vector<FromSize> active_tslices(Lt);
      for (const auto& it : params.param.sink_source_pairs)
      {
	// Check t_source and t_sink
	if (it.t_source < 0 || it.t_sink < 0)
	  throw std::runtime_error("Invalid source or sink on SinkSourcePairs");

	int num_tslices_active = it.Nt_backward + it.Nt_forward + 1;
	// Make the number of time-slices even; required by SB::doMomGammaDisp_contractions
	num_tslices_active = std::min(num_tslices_active + num_tslices_active % 2, Lt);

	FromSize fs = active_tslices[it.t_source % Lt];
	SB::union_interval(fs.from, fs.size, it.t_source - it.Nt_backward, num_tslices_active, Lt,
			   fs.from, fs.size);
	active_tslices[it.t_source % Lt] = fs;
	fs = active_tslices[it.t_sink % Lt];
	SB::union_interval(fs.from, fs.size, it.t_source - it.Nt_backward, num_tslices_active, Lt,
			   fs.from, fs.size);
	active_tslices[it.t_sink % Lt] = fs;
      }

      //
      // Store how many times a sink/source is call
      //
      std::vector<unsigned int> edges_on_tslice(Lt);
      for (const auto& it : params.param.sink_source_pairs)
      {
	edges_on_tslice[it.t_source % Lt]++;
	edges_on_tslice[it.t_sink % Lt]++;
      }

      //
      // Store what tslices the user suggest to cache
      //
      for (const auto &it : params.param.prop_sources)
      {
	// Check t_source
	if (it.t_source < 0)
	  throw std::runtime_error("Invalid source on PropSources");

	if (it.cacheP)
	  cache_tslice[it.t_source % Lt] = true;
      }

      //
      // Premultiply the gammas by g5
      //
      // NOTE: ultimately, we are using gamma5 hermiticity to change the propagator from source at time
      // t_sink to, instead, the t_slice
      //
      // Also NOTE: the gamma5 is hermitian. We will put the gamma_5 into the insertion, but since the need for the G5
      // is a part of the sink solution vector, we will multiply here.
      //
      // NOTE: with suitable signs, the gamma_5 could be merged with the gamma

      const int g5 = Ns * Ns - 1;
      std::vector<SB::Tensor<2, SB::Complex>> gamma_mats;
      {
	for (const int g : gammas)
	{
	  SpinMatrix gmat = Gamma(g5) * (Gamma(g) * SB::SpinMatrixIdentity());
	  gamma_mats.push_back(SB::asTensorView(gmat).cloneOn<SB::Complex>(SB::OnDefaultDevice));
	}
      }

      // Set how many processes are going to write elementals; each process is going to write in a
      // independent file
      bool use_multiple_writers = params.param.contract.use_multiple_writers;
      if (params.param.contract.use_genprop5_format)
	use_multiple_writers = true;

      //
      // DB storage
      //
      std::vector<LocalBinaryStoreDB<LocalSerialDBKey<KeyUnsmearedMesonElementalOperator_t>,
				     LocalSerialDBData<ValUnsmearedMesonElementalOperator_t>>>
	qdp_db;
      std::vector<LocalBinaryStoreDB<LocalSerialDBKey<KeyGenProp4ElementalOperator_t>,
				     LocalSerialDBData<ValGenProp4ElementalOperator_t>>>
	qdp4_db;
      SB::StorageTensor<10, SB::ComplexD> qdp5_db; // nNsSgdmtpP

      // Estimate the number of keys
      std::size_t max_tslices = 0;
      for (const auto& sink_source : params.param.sink_source_pairs)
	max_tslices =
	  std::max(max_tslices, (std::size_t)sink_source.Nt_backward + sink_source.Nt_forward + 1);
      std::size_t num_keys_gp4 = phases.numMom() * gammas.size() * disps.size() * max_tslices *
				 params.param.sink_source_pairs.size();
      for (auto& db : qdp_db)
	db.setNumberBuckets(num_keys_gp4 * num_vecs * 2);
      for (auto& db : qdp4_db)
	db.setNumberBuckets(num_keys_gp4 * 2);

      // The final elementals are going to be distributed along the lattice `t`
      // dimension, with no support on the lattice spatial dimension.  Because
      // of this, not all processes are going to have support on the final
      // elementals. The processes that have are going to write them on disk

      bool db_is_open = false;  //< whether qdp_db/qdp4_db has been opened
      int this_proc_id_t = -1;	//< the process id on the tensor holding the elementals
      // This function open the output file, after changing the name with the process id if multiple writers is used
      // \param proc_id_t: process rank on the tensor
      // \param numprocs_t: number of processes with support on the tensor
      std::function<void(int, int)> open_db = [&](int proc_id_t, int numprocs_t) {
	if (params.param.contract.use_genprop5_format)
	  return;

	// If this process has not support on the tensor, do nothing
	if (proc_id_t < 0)
	  return;

	if (db_is_open)
	{
	  assert(proc_id_t == this_proc_id_t);
	  assert((!params.param.contract.use_genprop4_format && qdp_db.size() == 1) ||
		 (params.param.contract.use_genprop4_format && qdp4_db.size() == 1));
	  return;
	}
	this_proc_id_t = proc_id_t;
	db_is_open = true;

	// If the final elementals are going to be spread among several processes, append the index
	// of the current process on the `t` dimension to the filename
	if (!params.param.contract.use_genprop4_format)
	  qdp_db.resize(1);
	else
	  qdp4_db.resize(1);
	std::string filename = params.named_obj.dist_op_file;
	if (use_multiple_writers)
	  filename += "." + std::to_string(proc_id_t + 1) + "_outof_" + std::to_string(numprocs_t);

	// Open the file, and write the meta-data and the binary for this operator
	if (!params.param.contract.use_genprop4_format)
	{
	  if (!qdp_db[0].fileExists(filename))
	  {
	    XMLBufferWriter file_xml;

	    push(file_xml, "DBMetaData");
	    write(file_xml, "id", std::string("unsmearedMesonElemOp"));
	    write(file_xml, "lattSize", QDP::Layout::lattSize());
	    write(file_xml, "decay_dir", decay_dir);
	    proginfo(file_xml); // Print out basic program info
	    write(file_xml, "Config_info", gauge_xml);
	    pop(file_xml);

	    std::string file_str(file_xml.str());
	    qdp_db[0].setMaxUserInfoLen(file_str.size());

	    qdp_db[0].open(filename, O_RDWR | O_CREAT, 0664);

	    qdp_db[0].insertUserdata(file_str);
	  }
	  else
	  {
	    qdp_db[0].open(filename, O_RDWR, 0664);
	  }
	}
	else
	{
	  if (!qdp4_db[0].fileExists(filename))
	  {
	    XMLBufferWriter file_xml;

	    push(file_xml, "DBMetaData");
	    write(file_xml, "id", std::string("genprop4ElemOp"));
	    write(file_xml, "lattSize", QDP::Layout::lattSize());
	    write(file_xml, "decay_dir", decay_dir);
	    proginfo(file_xml); // Print out basic program info
	    write(file_xml, "Config_info", gauge_xml);
	    pop(file_xml);

	    std::string file_str(file_xml.str());
	    qdp4_db[0].setMaxUserInfoLen(file_str.size());

	    qdp4_db[0].open(filename, O_RDWR | O_CREAT, 0664);

	    qdp4_db[0].insertUserdata(file_str);
	  }
	  else
	  {
	    qdp4_db[0].open(filename, O_RDWR, 0664);
	  }
	}
	QDPIO::cout << "Distillation file(s) opened" << std::endl;
      };

      if (params.param.contract.use_genprop5_format)
      {
	const char* qdp5_order = "nNsqgdmtpP";
	XMLBufferWriter metadata_xml;
	push(metadata_xml, "DBMetaData");
	write(metadata_xml, "id", std::string("genprop5ElemOp"));
	write(metadata_xml, "lattSize", QDP::Layout::lattSize());
	write(metadata_xml, "decay_dir", decay_dir);
	proginfo(metadata_xml); // Print out basic program info
	write(metadata_xml, "Config_info", gauge_xml);
	write(metadata_xml, "tensorOrder", qdp5_order);
	write(metadata_xml, "displacements", disps);
	std::vector<multi1d<int>>  moms;
	for (int i = 0; i < phases.numMom(); ++i)
	  moms.push_back(phases.numToMom(i));
	write(metadata_xml, "moms", moms);
	write(metadata_xml, "mass_label", params.param.contract.mass_label);
	write(metadata_xml, "gammas", gammas);
	write(metadata_xml, "eigen_phase", params.param.contract.phase);
	pop(metadata_xml);

	// NOTE: metadata_xml only has a valid value on Master node; so do a broadcast
	std::string metadata = SB::broadcast(metadata_xml.str());

	qdp5_db =
	  SB::StorageTensor<10, SB::ComplexD>(params.named_obj.dist_op_file, metadata, qdp5_order,
					      SB::kvcoors<10>(qdp5_order, {{'n', num_vecs},
									   {'N', num_vecs},
									   {'s', Ns},
									   {'q', Ns},
									   {'g', gammas.size()},
									   {'d', disps.size()},
									   {'m', moms.size()},
									   {'t', Lt},
									   {'p', Lt},
									   {'P', Lt}}),
					      SB::Sparse, SB::checksum_type::BlockChecksum);
	qdp5_db.preallocate(num_keys_gp4 * num_vecs * num_vecs * gammas.size() * sizeof(SB::ComplexD));
      }

      //
      // Try the factories
      //
      try
      {
	StopWatch swatch;
	swatch.reset();

	// Typedefs to save typing
	typedef LatticeFermion               T;
	typedef multi1d<LatticeColorMatrix>  P;
	typedef multi1d<LatticeColorMatrix>  Q;

	//
	// Initialize fermion action
	//
	SB::ChimeraSolver PP{params.param.prop.fermact, params.param.prop.invParam, u};

	//
	// Loop over the source color and spin, creating the source
	// and calling the relevant propagator routines.
	//

	std::vector<SB::Tensor<Nd + 5, SB::Complex>> invCache(Lt); // cache inversions

	// Maximum number of linear system RHS solved at once 
	const int max_rhs = params.param.contract.max_rhs;

	// Maximum number of tslices contracted at once (it has to be even)
	int max_tslices_in_contraction = params.param.contract.max_tslices_in_contraction;
	if (max_tslices_in_contraction <= 0)
	  max_tslices_in_contraction = Lt;
	max_tslices_in_contraction = max_tslices_in_contraction + (max_tslices_in_contraction % 2);
	max_tslices_in_contraction = std::min(Lt, max_tslices_in_contraction);

	// Maximum number of momenta contracted at once
	int max_moms_in_contraction = params.param.contract.max_moms_in_contraction;
	if (max_moms_in_contraction <= 0)
	  max_moms_in_contraction = phases.numMom();

	// Set place for doing the contractions
	SB::DeviceHost dev =
	  params.param.contract.use_device_for_contractions ? SB::OnDefaultDevice : SB::OnHost;

	for (const auto& sink_source : params.param.sink_source_pairs)
	{
	  int t_sink         = sink_source.t_sink % Lt;
	  int t_source       = sink_source.t_source % Lt;

	  QDPIO::cout << "\n\n--------------------------\nSink-Source pair: t_sink = " << t_sink << " t_source = " << t_source << std::endl; 
	  swatch.reset();
	  swatch.start();

	  int first_tslice_active = t_source - sink_source.Nt_backward;
	  int num_tslices_active =
	    std::min(sink_source.Nt_backward + std::max(sink_source.Nt_forward, 1), Lt);
	  // Make the number of time-slices even; required by SB::doMomGammaDisp_contractions
	  num_tslices_active = std::min(num_tslices_active + num_tslices_active % 2, Lt);

	  if (!invCache[t_source])
	  {
	    // If this inversion is not going to be cache, just store tslices for this source-sink pair
	    if (!cache_tslice[t_source])
	    {
	      active_tslices[t_source].from = first_tslice_active;
	      active_tslices[t_source].size = num_tslices_active;
	    }

	    // Get num_vecs colorvecs on time-slice t_source
	    SB::Tensor<Nd + 3, SB::ComplexF> source_colorvec = SB::getColorvecs(
	      colorvecsSto, u, decay_dir, t_source, 1, num_vecs, SB::none, phase, dev);

	    // Invert the source for all spins and retrieve num_tslices_active
	    // time-slices starting from time-slice first_tslice_active
	    invCache[t_source] = SB::doInversion<SB::ComplexF, SB::Complex>(
	      PP, std::move(source_colorvec), t_source, active_tslices[t_source].from,
	      active_tslices[t_source].size, {0, 1, 2, 3}, max_rhs, "cxyzXnSst");
	  }

	  if (!invCache[t_sink])
	  {
	    // If this inversion is not going to be cache, just store tslices for this source-sink pair
	    if (!cache_tslice[t_sink])
	    {
	      active_tslices[t_sink].from = first_tslice_active;
	      active_tslices[t_sink].size = num_tslices_active;
	    }

	    // Get num_vecs colorvecs on time-slice t_sink
	    SB::Tensor<Nd + 3, SB::ComplexF> sink_colorvec = SB::getColorvecs(
	      colorvecsSto, u, decay_dir, t_sink, 1, num_vecs, SB::none, phase, dev);

	    // Invert the sink for all spins and retrieve num_tslices_active time-slices starting from
	    // time-slice first_tslice_active
	    invCache[t_sink] = SB::doInversion<SB::ComplexF, SB::Complex>(
	      PP, std::move(sink_colorvec), t_sink, active_tslices[t_sink].from,
	      active_tslices[t_sink].size, {0, 1, 2, 3}, max_rhs, "ScnsxyzXt");
	  }

	  // The cache may have more tslices than need it; restrict to the ones required for this source-sink pair
	  SB::Tensor<Nd + 5, SB::Complex> invSource = invCache[t_source].kvslice_from_size(
	    {{'t', first_tslice_active - active_tslices[t_source].from}},
	    {{'t', num_tslices_active}});
	  SB::Tensor<Nd + 5, SB::Complex> invSink = invCache[t_sink].kvslice_from_size(
	    {{'t', first_tslice_active - active_tslices[t_sink].from}},
	    {{'t', num_tslices_active}});

	  // Remove from cache the source/sink inversions if the user suggests it or they are not going to be used anymore
	  edges_on_tslice[t_source]--;
	  edges_on_tslice[t_sink]--;
	  if (edges_on_tslice[t_source] == 0 || !cache_tslice[t_source])
	    invCache[t_source].release();
	  if (edges_on_tslice[t_sink] == 0 || !cache_tslice[t_sink])
	    invCache[t_sink].release();

	  // Contract the spatial components of sink and source together with
	  // several momenta, gammas and displacements; but contract not more than
	  // max_tslices_in_contraction at once!

	  invSink = invSink.rename_dims({{'n', 'N'}, {'s', 'q'}, {'S', 'Q'}});
	  for (int tfrom = 0, tsize = std::min(max_tslices_in_contraction, num_tslices_active);
	       tfrom < num_tslices_active; tfrom += tsize,
		   tsize = std::min(max_tslices_in_contraction, num_tslices_active - tfrom))
	  {
	    for (int mfrom = 0, msize = std::min(max_moms_in_contraction, phases.numMom());
		 mfrom < phases.numMom();
		 mfrom += msize, msize = std::min(max_moms_in_contraction, phases.numMom() - mfrom))
	    {

	      StopWatch snarss1;
	      snarss1.reset();
	      snarss1.start();

	      const char order_out[] = "qgmNndst";
	      SB::Tensor<Nd + 5, SB::Complex> this_invSource =
		invSource.kvslice_from_size({{'t', tfrom}}, {{'t', tsize}});
	      SB::Tensor<Nd + 5, SB::Complex> this_invSink =
		invSink.kvslice_from_size({{'t', tfrom}}, {{'t', tsize}});
	      if (tfrom + tsize >= num_tslices_active && mfrom + msize >= phases.numMom())
	      {
		invSource.release();
		invSink.release();
	      }
	      std::pair<SB::Tensor<8, SB::Complex>, std::vector<int>> r =
		SB::doMomGammaDisp_contractions<8>(
		  u, std::move(this_invSink), std::move(this_invSource),
		  first_tslice_active + tfrom, phases, mfrom, msize, gamma_mats, disps,
		  params.param.contract.use_derivP, order_out, SB::none, dev);

	      // Premultiply by g5, again; see above commit about this
	      SB::Tensor<8, SB::Complex> g5_con = r.first.like_this(
		"qgmNndst", {}, SB::OnHost, use_multiple_writers ? SB::OnEveryone : SB::OnMaster);
	      g5_con.contract(SB::Gamma<SB::Complex>(g5, SB::OnDefaultDevice), {}, SB::NotConjugate,
			      std::move(r.first), {{'q', 'j'}}, SB::NotConjugate, {{'q', 'i'}});
	      const std::vector<int> disps_perm = r.second;

	      snarss1.stop();
	      QDPIO::cout << "Time to compute contractions for " << tsize
			  << " tslices from t= " << (first_tslice_active + tfrom) % Lt << " and "
			  << msize << " momenta from momentum " << mfrom << " : "
			  << snarss1.getTimeInSeconds() << " secs" << std::endl;

	      //
	      // Write the elementals
	      //

	      snarss1.reset();
	      snarss1.start();

	      if (params.param.contract.use_genprop5_format)
	      {
		auto g5_con_rearrange_d = g5_con.like_this();
		for (int d = 0; d < disps_perm.size(); ++d)
		{
		  g5_con.kvslice_from_size({{'d', d}}, {{'d', 1}})
		    .copyTo(
		      g5_con_rearrange_d.kvslice_from_size({{'d', disps_perm[d]}}, {{'d', 1}}));
		}

		qdp5_db
		  .kvslice_from_size({{'m', mfrom},
				      {'t', (tfrom + first_tslice_active) % Lt},
				      {'p', t_source},
				      {'P', t_sink}},
				     {{'p', 1}, {'P', 1}})
		  .copyFrom(g5_con_rearrange_d);
	      }
	      else
	      {
		// Open DB if they are not opened already
		open_db(g5_con.p->procRank(), g5_con.p->numProcs());

		// Store the tensor
		if (!params.param.contract.use_genprop4_format)
		{
		  // Store
		  LocalSerialDBKey<KeyUnsmearedMesonElementalOperator_t> key;
		  LocalSerialDBData<ValUnsmearedMesonElementalOperator_t> val;
		  val.data() = ValUnsmearedMesonElementalOperator_t(num_vecs);

		  for (int t = 0; t < tsize; ++t)
		  {
		    for (int g = 0; g < gammas.size(); ++g)
		    {
		      for (int mom = 0; mom < msize; ++mom)
		      {
			for (int d = 0; d < disps_perm.size(); ++d)
			{
			  for (int n = 0; n < num_vecs; ++n)
			  {
			    auto g5_con_t =
			      g5_con
				.kvslice_from_size(
				  {{'g', g}, {'m', mom}, {'n', n}, {'d', d}, {'t', t}},
				  {{'g', 1}, {'m', 1}, {'n', 1}, {'d', 1}, {'t', 1}})
				.getLocal();
			    if (g5_con_t)
			    {
			      g5_con_t.copyTo(val.data());

			      key.key().derivP = params.param.contract.use_derivP;
			      key.key().t_sink = t_sink;
			      key.key().t_slice = (t + tfrom + first_tslice_active) % Lt;
			      key.key().t_source = t_source;
			      key.key().colorvec_src = n;
			      key.key().gamma = gammas[g];
			      key.key().displacement = disps[disps_perm[d]];
			      key.key().mom = phases.numToMom(mfrom + mom);
			      key.key().mass = params.param.contract.mass_label;

			      qdp_db[use_multiple_writers ? mfrom + mom : 0].insert(key, val);
			    }
			  }
			}
		      }
		    }
		  }
		}
		else
		{
		  // Store
		  LocalSerialDBKey<KeyGenProp4ElementalOperator_t> key;
		  LocalSerialDBData<ValGenProp4ElementalOperator_t> val;
		  val.data() = ValGenProp4ElementalOperator_t(num_vecs, num_vecs);

		  for (int t = 0; t < tsize; ++t)
		  {
		    for (int g = 0; g < gammas.size(); ++g)
		    {
		      for (int mom = 0; mom < msize; ++mom)
		      {
			for (int d = 0; d < disps_perm.size(); ++d)
			{
			  auto g5_con_t =
			    g5_con
			      .kvslice_from_size({{'g', g}, {'m', mom}, {'d', d}, {'t', t}},
						 {{'g', 1}, {'m', 1}, {'d', 1}, {'t', 1}})
			      .getLocal();

			  if (g5_con_t)
			  {
			    g5_con_t.copyTo(val.data());

			    key.key().t_sink = t_sink;
			    key.key().t_slice = (t + tfrom + first_tslice_active) % Lt;
			    key.key().t_source = t_source;
			    key.key().g = gammas[g];
			    key.key().displacement = disps[disps_perm[d]];
			    key.key().mom = phases.numToMom(mfrom + mom);
			    key.key().mass = params.param.contract.mass_label;

			    qdp4_db[use_multiple_writers ? mfrom + mom : 0].insert(key, val);
			  }
			}
		      }
		    }
		  }
		}
	      }
	      snarss1.stop();
	      QDPIO::cout << "Time to store " << tsize
			  << " tslices : " << snarss1.getTimeInSeconds() << " secs" << std::endl;
	    }
	  }
	  swatch.stop();
	  QDPIO::cout << "SINK-SOURCE: time to compute all source solution vectors and insertions "
			 "for t_sink= "
		      << t_sink << "  t_source= " << t_source
		      << "  time= " << swatch.getTimeInSeconds() << " secs" << std::endl;
	} // for sink_source
      }
      catch (const std::exception& e) 
      {
	QDP_error_exit("%s: caught exception: %s\n", name.c_str(), e.what());
      }

      // Close db
      if (db_is_open)
      {
	for (auto& db : qdp_db)
	  db.close();
	for (auto& db : qdp4_db)
	  db.close();
      }

      // Close colorvecs storage
      SB::closeColorvecStorage(colorvecsSto);

      // Close the xml output file
      pop(xml_out);     // UnsmearedHadronNode

      snoop.stop();
      QDPIO::cout << name << ": total time = " 
		  << snoop.getTimeInSeconds() 
		  << " secs" << std::endl;

      QDPIO::cout << name << ": ran successfully" << std::endl;

      END_CODE();
    } // func
  } // namespace InlineUnsmearedHadronNodeDistillationSuperbEnv

  /*! @} */  // end of group hadron

} // namespace Chroma

#endif // BUILD_SB
