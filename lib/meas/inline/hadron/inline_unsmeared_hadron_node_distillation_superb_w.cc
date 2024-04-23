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

      input.alt_t_start = 0;
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
      read(inputtop, "mass_label", input.mass_label);

      input.do_summation = false;
      if (inputtop.count("do_summation") == 1)
      {
	read(inputtop, "do_summation", input.do_summation);
      }

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

      input.output_file_is_local = false;
      if( inputtop.count("output_file_is_local") == 1 ) {
        read(inputtop, "output_file_is_local", input.output_file_is_local);
      }

      if (inputtop.count("phase") == 1)
      {
	read(inputtop, "phase", input.quarkPhase);
	if (inputtop.count("quarkPhase") == 1 || inputtop.count("quarkPhase") == 1)
	{
	  QDPIO::cerr << "Error: please don't give the tag `phase' and either `quarkPhase' or "
			 "`aQuarkPhase'"
		      << std::endl;
	  QDP_abort(1);
	}
      }
      else if (inputtop.count("quarkPhase") == 1)
      {
	read(inputtop, "quarkPhase", input.quarkPhase);
      }
      else if (inputtop.count("aQuarkPhase") == 1)
      {
	QDPIO::cerr << "Label `aQuarkPhase' without the label `quarkPhase'" << std::endl;
	QDP_abort(1);
      }
      else
      {
	input.quarkPhase.resize(Nd - 1);
      }

      if (inputtop.count("aQuarkPhase") == 1)
      {
	read(inputtop, "aQuarkPhase", input.aQuarkPhase);
      }
      else
      {
	for (float i : input.quarkPhase)
	  input.aQuarkPhase.push_back(-i);
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
      write(xml, "mass_label", input.mass_label);
      write(xml, "max_rhs", input.max_rhs);
      write(xml, "max_tslices_in_contraction", input.max_tslices_in_contraction);
      write(xml, "max_moms_in_contraction", input.max_moms_in_contraction);
      write(xml, "use_genprop4_format", input.use_genprop4_format);
      write(xml, "use_genprop5_format", input.use_genprop5_format);
      write(xml, "output_file_is_local", input.output_file_is_local);
      write(xml, "use_device_for_contractions", input.use_device_for_contractions);
      write(xml, "quarkPhase", SB::tomulti1d(input.quarkPhase));
      write(xml, "aQuarkPhase", SB::tomulti1d(input.aQuarkPhase));

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
      SB::CoorMoms moms;
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

	moms = SB::CoorMoms(moms_set.begin(), moms_set.end());
        disps.resize(disps_set.size());
        std::copy(disps_set.begin(), disps_set.end(), disps.begin());
        gammas.resize(gammas_set.size());
        std::copy(gammas_set.begin(), gammas_set.end(), gammas.begin());
      }

      // Get the maximum steps in the time direction
      int t_extra = 0;
      for (const auto& disp : disps)
      {
	int this_t_extra = 0;
	for (const auto& dir : disp)
	{
	  if (std::abs(dir) == 4)
	  {
	    this_t_extra += (dir < 0 ? -1 : 1);
	    t_extra = std::max(t_extra, std::abs(this_t_extra));
	  }
	}
      }

      //
      // Parse the phase
      //
      if (params.param.contract.quarkPhase.size() != Nd - 1 || params.param.contract.aQuarkPhase.size() != Nd - 1)
      {
	QDPIO::cerr << "`phase', `quarkPhase', and `aQuarkPhase' tags should have " << Nd - 1
		    << " components" << std::endl;
	QDP_abort(1);
      }
      SB::Coor<Nd - 1> negSinkPhase, sourcePhase;
      for (int i = 0; i < Nd - 1; ++i)
      {
	if (std::fabs((int)params.param.contract.quarkPhase[i] - params.param.contract.quarkPhase[i]) > 0 ||
	    std::fabs((int)params.param.contract.aQuarkPhase[i] - params.param.contract.aQuarkPhase[i]) > 0)
	  std::runtime_error("phase', `quarkPhase', and `aQuarkPhase' should be integer");
	sourcePhase[i] = params.param.contract.quarkPhase[i];
	negSinkPhase[i] = -params.param.contract.aQuarkPhase[i];
      }

      //
      // Capture maximum number of vecs
      //
      int num_vecs = params.param.contract.alt_num_vecs;
      for (const auto& it : params.param.prop_sources)
	num_vecs = std::max(num_vecs, it.num_vecs);

      //
      // Stores the range of time-slices used for each sink/source
      //

      std::vector<bool> cache_tslice(Lt, true);
      for (const auto& it : params.param.alt_sink_sources)
      {
	for (const auto& snk : it.second)
	{
	  if (it.first < 0 || snk < 0)
	    throw std::runtime_error("Invalid source or sink on SinkSources");
	  Params::Param_t::SinkSource_t ss;
	  ss.t_sink = snk % Lt;
	  ss.t_source = it.first % Lt;
	  if (!params.param.contract.do_summation)
	  {
	    ss.Nt_backward = it.first - params.param.contract.alt_t_start;
	    ss.Nt_forward = params.param.contract.alt_Nt_forward - ss.Nt_backward;
	  }
	  else
	  {
	    // Compute from src+1 up to snk-1
	    ss.Nt_backward = -1;
	    ss.Nt_forward = std::max(SB::normalize_coor(ss.t_sink - ss.t_source, Lt) + 1 - 2, 0);
	  }
	  params.param.sink_source_pairs.push_back(ss);
	}
      }

      struct FromSize {
	int from;
	int size;
      };
      std::vector<FromSize> active_tslices_source(Lt), active_tslices_sink0(Lt);
      std::vector<FromSize>& active_tslices_sink =
	  negSinkPhase == sourcePhase ? active_tslices_source : active_tslices_sink0;
      for (auto& it : params.param.sink_source_pairs)
      {
	// Check t_source and t_sink
	if (it.t_source < 0 || it.t_sink < 0)
	  throw std::runtime_error("Invalid source or sink on SinkSourcePairs");

	if (params.param.contract.do_summation)
	{
	  // Compute from src+1 up to snk-1
	  it.Nt_backward = -1;
	  it.Nt_forward = std::max(SB::normalize_coor(it.t_sink - it.t_source, Lt) + 1 - 2, 0);
	}

	int num_tslices_active = std::min(it.Nt_backward + it.Nt_forward + 1 + 2 * t_extra, Lt);
	// Make the number of time-slices even; required by SB::doMomGammaDisp_contractions
	num_tslices_active = std::min(num_tslices_active + num_tslices_active % 2, Lt);

	FromSize fs = active_tslices_source[it.t_source % Lt];
	SB::union_interval(fs.from, fs.size, it.t_source - it.Nt_backward - t_extra,
			   num_tslices_active, Lt, fs.from, fs.size);
	active_tslices_source[it.t_source % Lt] = fs;
	fs = active_tslices_sink[it.t_sink % Lt];
	SB::union_interval(fs.from, fs.size, it.t_source - it.Nt_backward - t_extra,
			   num_tslices_active, Lt, fs.from, fs.size);
	active_tslices_sink[it.t_sink % Lt] = fs;
      }

      //
      // Store how many times a sink/source is call
      //
      std::vector<unsigned int> edges_on_tslice_source(Lt), edges_on_tslice_sink0(Lt);
      std::vector<unsigned int>& edges_on_tslice_sink =
	negSinkPhase == sourcePhase ? edges_on_tslice_source : edges_on_tslice_sink0;
      for (const auto& it : params.param.sink_source_pairs)
      {
	edges_on_tslice_source[it.t_source % Lt]++;
	edges_on_tslice_sink[it.t_sink % Lt]++;
      }

      //
      // Store what tslices the user suggest to cache
      //
      for (const auto &it : params.param.prop_sources)
      {
	// Check t_source
	if (it.t_source < 0)
	  throw std::runtime_error("Invalid source on PropSources");

	if (!it.cacheP)
	  cache_tslice[it.t_source % Lt] = false;
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
	max_tslices = std::max(max_tslices,
			       (std::size_t)(sink_source.Nt_backward + sink_source.Nt_forward + 1));
      std::size_t num_keys_gp4 = moms.size() * gammas.size() * disps.size() * max_tslices *
				 params.param.sink_source_pairs.size();

      bool db_is_open = false;  //< whether qdp_db/qdp4_db has been opened

      // This function open the output file when using filehash
      auto open_db = [&]() {
	if (db_is_open)
	  return;
	db_is_open = true;

	std::string filename = params.named_obj.dist_op_file;

	// Open the file, and write the meta-data and the binary for this operator
	if (!params.param.contract.use_genprop4_format)
	{
	  qdp_db.resize(1);
	  qdp_db[0].setNumberBuckets(num_keys_gp4 * num_vecs * 2);
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
	  qdp4_db.resize(1);
	  qdp4_db[0].setNumberBuckets(num_keys_gp4 * 2);
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
	std::vector<multi1d<int>>  moms0;
	for (int i = 0; i < moms.size(); ++i)
	  moms0.push_back(SB::tomulti1d(moms[i]));
	write(metadata_xml, "moms", moms0);
	write(metadata_xml, "mass_label", params.param.contract.mass_label);
	write(metadata_xml, "gammas", gammas);
	write(metadata_xml, "quarkPhase", SB::tomulti1d(params.param.contract.quarkPhase));
	write(metadata_xml, "aQuarkPhase", SB::tomulti1d(params.param.contract.aQuarkPhase));
	pop(metadata_xml);

	// NOTE: metadata_xml only has a valid value on Master node; so do a broadcast
	std::string metadata = SB::broadcast(metadata_xml.str());

	qdp5_db = SB::StorageTensor<10, SB::ComplexD>(
	  params.named_obj.dist_op_file, metadata, qdp5_order,
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
	  SB::Sparse, SB::checksum_type::BlockChecksum,
	  params.param.contract.output_file_is_local ? SB::LocalFSFile : SB::SharedFSFile);
	qdp5_db.preallocate(num_keys_gp4 * num_vecs * num_vecs * gammas.size() *
			    sizeof(SB::ComplexD) /
			    (params.param.contract.output_file_is_local ? Layout::numNodes() : 1));
      }

      // Initialize fermion action
      // NOTE: this gets out the following try-block because QUDA and MGPROTO solvers may
      // hang when an exception is thrown, preventing the report of the exception message
      SB::ChimeraSolver PP{params.param.prop.fermact, params.param.prop.invParam, u};

      // NOTE: qdp5_db needs MPI synchronization when closing, so capture exception and abort in that case
      //       to avoid hangs
      try
      {
	StopWatch swatch;

	//
	// Loop over the source color and spin, creating the source
	// and calling the relevant propagator routines.
	//

	std::vector<SB::Tensor<Nd + 5, SB::Complex>> invCacheSource(Lt), invCacheSink0(Lt); // cache inversions
	std::vector<SB::Tensor<Nd + 5, SB::Complex>>& invCacheSink =
	  negSinkPhase == sourcePhase ? invCacheSource : invCacheSink0;

	// Maximum number of linear system RHS solved at once 
	const int max_rhs = params.param.contract.max_rhs;

	// Maximum number of tslices contracted at once (it has to be even)
	int max_tslices_in_contraction = params.param.contract.max_tslices_in_contraction;
	if (!params.param.contract.do_summation)
	{
	  if (max_tslices_in_contraction <= 0)
	    max_tslices_in_contraction = Lt;
	  max_tslices_in_contraction = std::min(Lt, max_tslices_in_contraction);
	}
	else
	{
	  // When doing summation, compute all middle time slices at once
	  max_tslices_in_contraction = Lt;
	}

	// Maximum number of momenta contracted at once
	int max_moms_in_contraction = params.param.contract.max_moms_in_contraction;
	if (max_moms_in_contraction <= 0)
	  max_moms_in_contraction = moms.size();

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

	  int first_tslice_active; // first middle time-slice to compute
	  int num_tslices_active; // number of middle time-slices to compute
	  first_tslice_active =
	    SB::normalize_coor(t_source - sink_source.Nt_backward - t_extra, Lt);
	  num_tslices_active =
	    std::min(sink_source.Nt_backward +
		       (sink_source.Nt_forward == 0 ? 1 : sink_source.Nt_forward) + 2 * t_extra,
		     Lt);

	  // Make the number of time-slices even; required by SB::doMomGammaDisp_contractions
	  num_tslices_active = std::min(num_tslices_active + num_tslices_active % 2, Lt);

	  if (!invCacheSource[t_source])
	  {
	    // If this inversion is not going to be cache, just store tslices for this source-sink pair
	    if (!cache_tslice[t_source])
	    {
	      active_tslices_source[t_source].from = first_tslice_active;
	      active_tslices_source[t_source].size = num_tslices_active;
	    }

	    // Get num_vecs colorvecs on time-slice t_source
	    SB::Tensor<Nd + 3, SB::ComplexF> source_colorvec = SB::getColorvecs(
	      colorvecsSto, u, decay_dir, t_source, 1, num_vecs, SB::none, sourcePhase, dev);

	    // Invert the source for all spins and retrieve num_tslices_active
	    // time-slices starting from time-slice first_tslice_active
	    invCacheSource[t_source] = SB::doInversion(
	      PP, std::move(source_colorvec), t_source, active_tslices_source[t_source].from,
	      active_tslices_source[t_source].size, {0, 1, 2, 3}, max_rhs, "cxyzXnSst");
	  }

	  if (!invCacheSink[t_sink])
	  {
	    // If this inversion is not going to be cache, just store tslices for this source-sink pair
	    if (!cache_tslice[t_sink])
	    {
	      active_tslices_sink[t_sink].from = first_tslice_active;
	      active_tslices_sink[t_sink].size = num_tslices_active;
	    }

	    // Get num_vecs colorvecs on time-slice t_sink
	    SB::Tensor<Nd + 3, SB::ComplexF> sink_colorvec = SB::getColorvecs(
	      colorvecsSto, u, decay_dir, t_sink, 1, num_vecs, SB::none, negSinkPhase, dev);

	    // Invert the sink for all spins and retrieve num_tslices_active time-slices starting from
	    // time-slice first_tslice_active
	    invCacheSink[t_sink] = SB::doInversion(
	      PP, std::move(sink_colorvec), t_sink, active_tslices_sink[t_sink].from,
	      active_tslices_sink[t_sink].size, {0, 1, 2, 3}, max_rhs, "ScnsxyzXt");
	  }

	  // The cache may have more tslices than need it; restrict to the ones required for this source-sink pair
	  SB::Tensor<Nd + 5, SB::Complex> invSource = invCacheSource[t_source].kvslice_from_size(
	    {{'t', first_tslice_active - active_tslices_source[t_source].from}},
	    {{'t', num_tslices_active}});
	  SB::Tensor<Nd + 5, SB::Complex> invSink = invCacheSink[t_sink].kvslice_from_size(
	    {{'t', first_tslice_active - active_tslices_sink[t_sink].from}},
	    {{'t', num_tslices_active}});
	  invSink = invSink.rename_dims({{'n', 'N'}, {'s', 'q'}, {'S', 'Q'}});

	  // Remove from cache the source/sink inversions if the user suggests it or they are not going to be used anymore
	  edges_on_tslice_source[t_source]--;
	  edges_on_tslice_sink[t_sink]--;
	  if (edges_on_tslice_source[t_source] == 0 || !cache_tslice[t_source])
	    invCacheSource[t_source].release();
	  if (edges_on_tslice_sink[t_sink] == 0 || !cache_tslice[t_sink])
	    invCacheSink[t_sink].release();

	  double time_in_writing = 0; // time in writing in genprops

	  auto call =
	    [&](SB::Tensor<7, SB::Complex> r, int disp_index, int tfrom, int mfrom) {

	      // Premultiply by g5, again; see above comment about this
	      r = SB::contract<7>(r, SB::Gamma<SB::Complex>(g5, dev).rename_dims({{'j', 'q'}}), "q")
		    .rename_dims({{'i', 'q'}});

	      //
	      // Do summation over all time slices between t_source+1 and t_sink-1
	      //
	      int tsize = r.kvdim().at('t');
	      int msize = r.kvdim().at('m');
	      if (params.param.contract.do_summation)
	      {
		auto s = r.like_this(SB::none, {{'t', 1}});
		s.set_zero();
		bool something_was_sum_up = false;
		for (int t = 0; t < tsize; ++t)
		{
		  if (SB::normalize_coor(tfrom + t - (t_source + 1), Lt) <
		      SB::normalize_coor(t_sink - 1 - (t_source + 1), Lt))
		  {
		    r.kvslice_from_size({{'t', t}}, {{'t', 1}}).addTo(s);
		    something_was_sum_up = true;
		  }
		}
		if (!something_was_sum_up)
		  return;
		r = s;
	      }

	      //
	      // Write the elementals
	      //

	      double time_writing_this = -SB::w_time();

	      if (params.param.contract.use_genprop5_format)
	      {
		qdp5_db
		  .kvslice_from_size(
		    {{'m', mfrom},
		     {'d', disp_index},
		     {'t', !params.param.contract.do_summation ? tfrom : (t_source + 1) % Lt},
		     {'p', t_source},
		     {'P', t_sink}},
		    {{'p', 1}, {'P', 1}, {'d', 1}})
		  .copyFrom(r);
	      }
	      else
	      {
		// Move the result to the master node and do the writing (only the master node)
		r = r.make_sure(SB::none, SB::OnHost, SB::OnMaster).getLocal();
		if (!r)
		  return;

		// Open DB if they are not opened already
		open_db();

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
			for (int n = 0; n < num_vecs; ++n)
			{
			  r.kvslice_from_size({{'g', g}, {'m', mom}, {'n', n}, {'t', t}},
					      {{'g', 1}, {'m', 1}, {'n', 1}, {'t', 1}})
			    .copyTo(val.data());

			  key.key().derivP = params.param.contract.use_derivP;
			  key.key().t_sink = t_sink;
			  key.key().t_slice = SB::normalize_coor(
			    !params.param.contract.do_summation ? t + tfrom : t_source + 1, Lt);
			  key.key().t_source = t_source;
			  key.key().colorvec_src = n;
			  key.key().gamma = gammas[g];
			  key.key().displacement = disps[disp_index];
			  key.key().mom = SB::tomulti1d(moms[mfrom + mom]);
			  key.key().mass = params.param.contract.mass_label;

			  qdp_db[0].insert(key, val);
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
			r.kvslice_from_size({{'g', g}, {'m', mom}, {'t', t}},
					    {{'g', 1}, {'m', 1}, {'t', 1}})
			  .copyTo(val.data());

			key.key().t_sink = t_sink;
			key.key().t_slice = SB::normalize_coor(
			  !params.param.contract.do_summation ? t + tfrom : t_source + 1, Lt);
			key.key().t_source = t_source;
			key.key().g = gammas[g];
			key.key().displacement = disps[disp_index];
			key.key().mom = SB::tomulti1d(moms[mfrom + mom]);
			key.key().mass = params.param.contract.mass_label;

			qdp4_db[0].insert(key, val);
		      }
		    }
		  }
		}
	      }

	      time_in_writing += SB::w_time() + time_writing_this;
	    };

	  // Contract the spatial components of sink and source together with
	  // several momenta, gammas and displacements; but contract not more than
	  // max_tslices_in_contraction at once!

	  double time_contracting_and_writing = -SB::w_time();
	  SB::doMomGammaDisp_contractions<7, Nd + 5, Nd + 5, SB::Complex>(
	    u, std::move(invSink), std::move(invSource), first_tslice_active, t_extra,
	    num_tslices_active - 2 * t_extra, moms, gamma_mats, disps,
	    params.param.contract.use_derivP, call, "qgmNnst", max_tslices_in_contraction,
	    max_moms_in_contraction, dev);
	  time_contracting_and_writing += SB::w_time();

	  QDPIO::cout << "Time to contract: " << time_contracting_and_writing - time_in_writing
		      << " secs" << std::endl;
	  QDPIO::cout << "Time to store: " << time_in_writing << " secs" << std::endl;

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
