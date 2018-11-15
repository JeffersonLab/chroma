/*! \file
 * \brief Inline measurement that construct unsmeared hadron nodes using distillation
 */

// Reverted version from master
// NB: The master version is what was in Boram's branch. In master it was dated Nov 25,
// whereas the previous version I integrated from feture/unsmeared-node was dated Nov 23
// Hence this is the latest version
//
#include "meas/inline/hadron/inline_unsmeared_hadron_node_distillation_w.h"

#ifndef QDP_IS_QDPJIT

#include "qdp_map_obj_memory.h"
#include "qdp_map_obj_disk.h"
#include "qdp_map_obj_disk_multiple.h"
#include "qdp_disk_map_slice.h"
#include "handle.h"

#include "meas/inline/hadron/inline_unsmeared_hadron_node_distillation_w.h"

#include "meas/glue/mesplq.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/smear/link_smearing_aggregate.h"
#include "meas/smear/link_smearing_factory.h"
#include "util/ferm/key_timeslice_colorvec.h"
#include "util/ferm/disp_soln_cache.h"
#include "util/ferm/key_val_db.h"
#include "util/info/proginfo.h"
#include "util/ft/sftmom.h"
#include "util/ft/time_slice_set.h"
#include "io/xml_group_reader.h"
#include "meas/inline/make_xml_file.h"

#include "util/ferm/key_val_db.h"
#include "util/ferm/transf.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"
#include "meas/inline/make_xml_file.h"

#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{ 
  /*!
   * \ingroup hadron
   *
   * @{
   */
  namespace InlineUnsmearedHadronNodeDistillationEnv 
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
      read(inputtop, "Nt_forward", input.Nt_forward);
      read(inputtop, "Nt_backward", input.Nt_backward);
    }

    //! Propagator output
    void write(XMLWriter& xml, const std::string& path, const Params::Param_t::KeySolnProp_t& input)
    {
      push(xml, path);

      write(xml, "cacheP", input.cacheP);
      write(xml, "num_vecs", input.num_vecs);
      write(xml, "t_source", input.t_source);
      write(xml, "Nt_forward", input.Nt_forward);
      write(xml, "Nt_backward", input.Nt_backward);

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

      read(inputtop, "use_derivP", input.use_derivP);
      read(inputtop, "decay_dir", input.decay_dir);
      read(inputtop, "displacement_length", input.displacement_length);
      read(inputtop, "mass_label", input.mass_label);
      read(inputtop, "num_tries", input.num_tries);
    }

    //! Propagator output
    void write(XMLWriter& xml, const std::string& path, const Params::Param_t::Contract_t& input)
    {
      push(xml, path);

      write(xml, "use_derivP", input.use_derivP);
      write(xml, "decay_dir", input.decay_dir);
      write(xml, "displacement_length", input.displacement_length);
      write(xml, "mass_label", input.mass_label);
      write(xml, "num_tries", input.num_tries);

      pop(xml);
    }


    //! Propagator input
    void read(XMLReader& xml, const std::string& path, Params::Param_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "DispGammaMomList", input.disp_gamma_mom_list);
      read(inputtop, "Propagator", input.prop);
      read(inputtop, "PropSources", input.prop_sources);
      read(inputtop, "SinkSourcePairs", input.sink_source_pairs);
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
  //-------------------------------------------------------------------------------
  namespace InlineUnsmearedHadronNodeDistillationEnv
  { 
    //----------------------------------------------------------------------------
    // Convenience type
    typedef QDP::MapObjectDisk<KeyTimeSliceColorVec_t, TimeSliceIO<LatticeColorVectorF> > MOD_t;

    // Convenience type
    typedef QDP::MapObjectDiskMultiple<KeyTimeSliceColorVec_t, TimeSliceIO<LatticeColorVectorF> > MODS_t;

    // Convenience type
    typedef QDP::MapObjectMemory<KeyTimeSliceColorVec_t, SubLatticeColorVectorF> SUB_MOD_t;

  } // end namespace


  //-------------------------------------------------------------------------------
  namespace InlineUnsmearedHadronNodeDistillationEnv
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

    const std::string name = "UNSMEARED_HADRON_NODE_DISTILLATION";

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
      read(bin, param.Nt_forward);
      read(bin, param.Nt_backward);
    }

    //! KeySolnProp write
    void write(BinaryWriter& bin, const Params::Param_t::KeySolnProp_t& param)
    {
      write(bin, param.cacheP);
      write(bin, param.num_vecs);
      write(bin, param.t_source);
      write(bin, param.Nt_forward);
      write(bin, param.Nt_backward);
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

    //! Meson operator
    struct ValUnsmearedMesonElementalOperator_t
    {
      multi3d<ComplexD>  op;          /*!< Colorvector source and spin sink */
    };


    //----------------------------------------------------------------------------
    //! Holds key and value as temporaries
    struct KeyValUnsmearedMesonElementalOperator_t
    {
      SerialDBKey<KeyUnsmearedMesonElementalOperator_t>  key;
      SerialDBData<ValUnsmearedMesonElementalOperator_t> val;
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
      read(bin, param.op);    // the size is always written, even if 0
    }

    //! UnsmearedMesonElementalOperator write
    void write(BinaryWriter& bin, const ValUnsmearedMesonElementalOperator_t& param)
    {
      write(bin, param.op);
    }

 
    //----------------------------------------------------------------------------------
    //----------------------------------------------------------------------------
    //! Map holding expressions. Key is the mom, val is some unique counter
    typedef MapObjectMemory< multi1d<int>, int >  MapTwoQuarkMom_t;

    //----------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------
    //! Twoquark field
    struct KeyTwoQuarkGamma_t
    {
      KeyTwoQuarkGamma_t() {}
      KeyTwoQuarkGamma_t(int g) : gamma(g) {}

      int    gamma;   /*!< The gamma in Gamma(gamma) - a number from [0 ... Ns*Ns) */
    };

    //----------------------------------------------------------------------------
    //! Reader
    void read(BinaryReader& bin, KeyTwoQuarkGamma_t& op)
    {
      read(bin, op.gamma);
    }

    //! Writer
    void write(BinaryWriter& bin, const KeyTwoQuarkGamma_t& op)
    {
      write(bin, op.gamma);
    }


    //----------------------------------------------------------------------------
    //! Map holding expressions. Key is the two-quark, val is the mom
    typedef MapObjectMemory< KeyTwoQuarkGamma_t, MapTwoQuarkMom_t >  MapTwoQuarkGammaMom_t;


    //----------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------
    //! Twoquark field
    struct KeyTwoQuarkDisp_t
    {
      KeyTwoQuarkDisp_t() {}
      KeyTwoQuarkDisp_t(const std::vector<int>& d) : deriv(d) {}

      std::vector<int>   deriv;   /*!< Left-right derivative path. This is 1-based. */
    };

    //----------------------------------------------------------------------------
    //! Reader
    void read(BinaryReader& bin, KeyTwoQuarkDisp_t& op)
    {
      read(bin, op.deriv);
    }
  
    //! Writer
    void write(BinaryWriter& bin, const KeyTwoQuarkDisp_t& op)
    {
      write(bin, op.deriv);
    }

    //----------------------------------------------------------------------------
    //! Map holding expressions. Key is the two-quark, val is the weight.
    typedef MapObjectMemory<KeyTwoQuarkDisp_t, MapTwoQuarkGammaMom_t>  MapTwoQuarkDispGammaMom_t;


    //-------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------
    //! Invert off of each source, and do all the checking
    LatticeFermion doInversion(const SystemSolver<LatticeFermion>& PP, const LatticeColorVector& vec_srce, int spin_source, int num_tries)
    {
      //
      // Loop over each spin source and invert. 
      //
      LatticeFermion chi = zero;
      CvToFerm(vec_srce, chi, spin_source);

      LatticeFermion quark_soln = zero;
      int ncg_had;

      // Do the propagator inversion
      // Check if bad things are happening
      bool badP = true;
      for(int nn = 1; nn <= num_tries; ++nn)
      {	
	// Reset
	quark_soln = zero;
	badP = false;
	      
	// Solve for the solution std::vector
	SystemSolverResults_t res = PP(quark_soln, chi);
	ncg_had = res.n_count;

	// Some sanity checks
	if (toDouble(res.resid) > 1.0e-3)
	{
	  QDPIO::cerr << name << ": have a resid > 1.0e-3. That is unacceptable" << std::endl;
	  QDP_abort(1);
	}

	// Check for finite values - neither NaN nor Inf
	if (isfinite(quark_soln))
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
	QDPIO::cerr << name << ": this is bad - did not get a finite solution vector after num_tries= " << num_tries << std::endl;
	QDP_abort(1);
      }

      return quark_soln;
    }
    

    //-------------------------------------------------------------------------------
    class SourcePropCache
    {
    public:
      SourcePropCache(const multi1d<LatticeColorMatrix>& u,
		      const ChromaProp_t& prop, MODS_t& eigs, int num_tries);

      //! New t-slice
      void newTimeSource(const Params::Param_t::KeySolnProp_t& key_);

      //! The solution
      const LatticeColorVectorSpinMatrix& getSoln(int t_source, int colorvec_src);

      //! Getters
      int getNumVecs(int t_source);

    private:
      //! Eigenvectors
      MODS_t& eigen_source;

      //! Put it here
      int num_tries;

      //! t_source -> prop keys
      std::map<int, Params::Param_t::KeySolnProp_t> keys;

      //! Cache:  {key -> {colorvec -> LCVSM}}
      MapObjectMemory<Params::Param_t::KeySolnProp_t, std::map<int, LatticeColorVectorSpinMatrix> > cache;

      //! qprop
      Handle< SystemSolver<LatticeFermion> > PP;
    };


    //-------------------------------------------------------------------------------
    SourcePropCache::SourcePropCache(const multi1d<LatticeColorMatrix>& u,
				     const ChromaProp_t& prop, MODS_t& eigs_, int num_tries_) : eigen_source(eigs_), num_tries(num_tries_)
    {
      StopWatch swatch;
      swatch.reset();
      swatch.start();
      QDPIO::cout << __func__ << ": initialize the prop cache" << std::endl;

      // Typedefs to save typing
      typedef LatticeFermion               T;
      typedef multi1d<LatticeColorMatrix>  P;
      typedef multi1d<LatticeColorMatrix>  Q;

      //
      // Initialize fermion action
      //
      std::istringstream  xml_s(prop.fermact.xml);
      XMLReader  fermacttop(xml_s);
      QDPIO::cout << "FermAct = " << prop.fermact.id << std::endl;

      // Generic Wilson-Type stuff
      Handle< FermionAction<T,P,Q> >
	S_f(TheFermionActionFactory::Instance().createObject(prop.fermact.id,
							     fermacttop,
							     prop.fermact.path));

      Handle< FermState<T,P,Q> > state(S_f->createState(u));

      PP = S_f->qprop(state, prop.invParam);
      
      swatch.stop(); 
      QDPIO::cout << __func__ << ": finished initializing the prop cache"
		    << "  time = " << swatch.getTimeInSeconds() 
		    << " secs" << std::endl;
    }

  
    //-------------------------------------------------------------------------------
    //! New t-slice
    void SourcePropCache::newTimeSource(const Params::Param_t::KeySolnProp_t& key)
    {
      if (cache.exist(key)) {
	QDPIO::cerr << __func__ << ": cache entry already exists" << std::endl;
	QDP_abort(1);
      }

      QDPIO::cout << __func__ << ": new t_source = " << key.t_source << std::endl;

      //! Insert a new entry
      keys.insert(std::make_pair(key.t_source, key));

      //! Cache size will depend on cacheP flag
      cache.insert(key, std::map<int, LatticeColorVectorSpinMatrix>());
    }


    //-------------------------------------------------------------------------------
    //! New t-slice
    const LatticeColorVectorSpinMatrix& SourcePropCache::getSoln(int t_slice, int colorvec_ind)
    {
      if (keys.find(t_slice) == keys.end()) {
	QDPIO::cerr << __func__ << ": t_ind does not exist in cache = " << t_slice << std::endl;
	QDP_abort(1);
      }

      // Key
      const Params::Param_t::KeySolnProp_t& key = keys[t_slice];

      if (cache[key].find(colorvec_ind) != cache[key].end()) {
	//QDPIO::cout << __func__ << ": FOUND KEY - t_slice = " << t_slice << "  colorvec_ind = " << colorvec_ind << std::endl; 
	return cache[key].at(colorvec_ind);
      }
      else {
	QDPIO::cout << __func__ << ": CREATING KEY: t_slice = " << t_slice << "  colorvec_ind = " << colorvec_ind << std::endl; 
	cache[key].insert(std::make_pair(colorvec_ind, LatticeColorVectorSpinMatrix(zero)));
      }

      // Get the source vector
      StopWatch swatch;
      swatch.reset();
      swatch.start();

      LatticeColorVectorF vec_srce = zero;
      KeyTimeSliceColorVec_t src_key(t_slice, colorvec_ind);
      TimeSliceIO<LatticeColorVectorF> time_slice_io(vec_srce, t_slice);
      eigen_source.get(src_key, time_slice_io);

      // Loop over each spin source
      for(int spin_ind=0; spin_ind < Ns; ++spin_ind)
      {	
	QDPIO::cout << __func__ << ": do t_slice= " << t_slice << "  spin_ind= " << spin_ind << "  colorvec_ind= " << colorvec_ind << std::endl; 
	    
	StopWatch snarss1;
	snarss1.reset();
	snarss1.start();

	//
	// Loop over each spin source and invert. 
	// Use the same color vector source. No spin dilution will be used.
	//
	LatticeFermion tmp = doInversion(*PP, vec_srce, spin_ind, num_tries);

	// shove a fermion into a colorvec-spinmatrix 
	FermToProp(tmp, cache[key][colorvec_ind], spin_ind);

	snarss1.stop();
	QDPIO::cout << __func__ << ": time to compute prop for t_slice= " << t_slice
		    << "  spin_ind= " << spin_ind
		    << "  colorvec_ind= " << colorvec_ind
		    << "  time = " << snarss1.getTimeInSeconds() 
		    << " secs" << std::endl;
      } // for spin_src

      swatch.stop(); 
      QDPIO::cout << __func__ << ": FINISHED: total time for t_slice = " << t_slice << "  colorvec_ind = " << colorvec_ind
		  << "  time= " << swatch.getTimeInSeconds() << " secs" << std::endl;
      
      return cache[key][colorvec_ind];
    }

    //-------------------------------------------------------------------------------
    //! New t-slice
    int SourcePropCache::getNumVecs(int t_source)
    {
      if (keys.find(t_source) == keys.end()) {
	QDPIO::cerr << __func__ << ": t_source does not exist in cache = " << t_source << std::endl;
	QDP_abort(1);
      }

      // Key
      const Params::Param_t::KeySolnProp_t& key = keys[t_source];

      return key.num_vecs;
    }
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


      //
      // Read in the source along with relevant information.
      // 
      QDPIO::cout << "Snarf the source from a std::map object disk file" << std::endl;

      MODS_t eigen_source;
      eigen_source.setDebug(0);

      std::string eigen_meta_data;   // holds the eigenvalues

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

#if 0
      // Sanity check
      if (params.param.contract.num_vecs > eigen_source.size())
      {
	QDPIO::cerr << name << ": number of available eigenvectors is too small\n";
	QDP_abort(1);
      }
      QDPIO::cout << "Number of vecs available is large enough" << std::endl;
#endif

      //
      // Setup the gauge smearing
      //
      QDPIO::cout << "Initalize link smearing" << std::endl;
      Handle< LinkSmearing > linkSmearing;

      try
      {
	std::istringstream  xml_l(params.param.link_smearing.xml);
	XMLReader  linktop(xml_l);
	QDPIO::cout << "Link smearing type = " << params.param.link_smearing.id << std::endl;
	
	linkSmearing = TheLinkSmearingFactory::Instance().createObject(params.param.link_smearing.id,
								       linktop, 
								       params.param.link_smearing.path);
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << name << ": Caught Exception link smearing: " << e << std::endl;
	QDP_abort(1);
      }

      QDPIO::cout << "Link smearing successfully initialized" << std::endl;

    
      //
      // Read through the momentum list and find all the unique phases
      //
      QDPIO::cout << "Parse momentum list" << std::endl;
      
      //
      // Initialize the slow Fourier transform phases
      //
      SftMom phases(0, false, params.param.contract.decay_dir);

      //
      // Check displacements
      //
      MapTwoQuarkDispGammaMom_t disp_gamma_moms;

      if (params.param.disp_gamma_mom_list.size() > 0)
      {
	int mom_size = 0;
	MapObjectMemory<multi1d<int>, int> uniq_mom;
	for(auto ins = params.param.disp_gamma_mom_list.begin(); ins != params.param.disp_gamma_mom_list.end(); ++ins)
	{	
	  // Sort out momentum
	  QDPIO::cout << name << ": mom= " << ins->mom << std::endl;
	  uniq_mom.insert(ins->mom, 1);
	  QDPIO::cout << " found mom" << std::endl;
	  mom_size = ins->mom.size();
	  QDPIO::cout << " mom_size= " << mom_size << std::endl;

	  //
	  // Build maps of displacements, gammas and their moms
	  //
	  // Make sure displacement is something sensible
	  std::vector<int> disp = normDisp(ins->displacement);

	  // Insert unique combos
	  // Drat - don't have automatic uniqueness generator
	  MapTwoQuarkMom_t george;
	  MapTwoQuarkGammaMom_t fred;
	  george.insert(ins->mom, 1);
	  fred.insert(ins->gamma, george);

	  if (disp_gamma_moms.exist(disp))
	  {
	    if (disp_gamma_moms[disp].exist(ins->gamma))
	    {
	      disp_gamma_moms[disp][ins->gamma].insert(ins->mom, 1);
	    }
	    else
	    {
	      disp_gamma_moms[disp].insert(ins->gamma, george);
	    }
	  }
	  else
	  {
	    disp_gamma_moms.insert(disp, fred);
	  }
	}

	int num_mom = uniq_mom.size();
	QDPIO::cout << name << ": num_mom= " << num_mom << "  mom_size= " << mom_size << std::endl;
	multi2d<int> moms(num_mom,mom_size);
	int i = 0;
	for(auto mom = uniq_mom.begin(); mom != uniq_mom.end(); ++mom)
	  moms[i++] = mom->first;

	SftMom temp_phases(moms, params.param.contract.decay_dir);
	phases = temp_phases;
      }
      else
      {
	QDPIO::cerr << name << ": warning - you have an empty disp_gamma_mom_list. Will allow under your insistence." << std::endl;
	QDP_abort(1);
      }


      //
      // Smear the gauge field if needed
      //
      multi1d<LatticeColorMatrix> u_smr = u;

      try
      {
	std::istringstream  xml_l(params.param.link_smearing.xml);
	XMLReader  linktop(xml_l);
	QDPIO::cout << "Link smearing type = " << params.param.link_smearing.id << std::endl;
	
	
	Handle< LinkSmearing >
	  linkSmearing(TheLinkSmearingFactory::Instance().createObject(params.param.link_smearing.id,
								       linktop, 
								       params.param.link_smearing.path));

	(*linkSmearing)(u_smr);
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << name << ": Caught Exception link smearing: " << e << std::endl;
	QDP_abort(1);
      }

      // Record the smeared observables
      MesPlq(xml_out, "Smeared_Observables", u_smr);


      //
      // DB storage
      //
      BinaryStoreDB< SerialDBKey<KeyUnsmearedMesonElementalOperator_t>, SerialDBData<ValUnsmearedMesonElementalOperator_t> > qdp_db;

      // Open the file, and write the meta-data and the binary for this operator
      if (! qdp_db.fileExists(params.named_obj.dist_op_file))
      {
	XMLBufferWriter file_xml;

	push(file_xml, "DBMetaData");
	write(file_xml, "id", std::string("unsmearedMesonElemOp"));
	write(file_xml, "lattSize", QDP::Layout::lattSize());
	write(file_xml, "decay_dir", decay_dir);
	proginfo(file_xml);    // Print out basic program info
	write(file_xml, "Config_info", gauge_xml);
	pop(file_xml);

	std::string file_str(file_xml.str());
	qdp_db.setMaxUserInfoLen(file_str.size());

	qdp_db.open(params.named_obj.dist_op_file, O_RDWR | O_CREAT, 0664);

	qdp_db.insertUserdata(file_str);
      }
      else
      {
	qdp_db.open(params.named_obj.dist_op_file, O_RDWR, 0664);
      }

      QDPIO::cout << "Finished opening distillation file" << std::endl;


      //
      // Try the factories
      //
      try
      {
	StopWatch swatch;
	swatch.reset();

	// Cache manager
	QDPIO::cout << name << ": initialize the prop cache" << std::endl;
	SourcePropCache prop_cache(u, params.param.prop, eigen_source, params.param.contract.num_tries);

	// All the desired solutions
	QDPIO::cout << name << ": initialize the time sources" << std::endl;
	for(auto key = params.param.prop_sources.begin(); key != params.param.prop_sources.end(); ++key)
	{
	  prop_cache.newTimeSource(*key);
	}

	
	//
	// Loop over the source color and spin, creating the source
	// and calling the relevant propagator routines.
	//
	const int g5 = Ns*Ns-1;

	for(auto sink_source = params.param.sink_source_pairs.begin(); sink_source != params.param.sink_source_pairs.end(); ++sink_source)
	{
	  int t_sink         = sink_source->t_sink;
	  int t_source       = sink_source->t_source;
	  int sink_num_vecs  = prop_cache.getNumVecs(t_sink);
	  int srce_num_vecs  = prop_cache.getNumVecs(t_source);

	  QDPIO::cout << "\n\n--------------------------\nSink-Source pair: t_sink = " << t_sink << " t_source = " << t_source << std::endl; 
	  swatch.reset();
	  swatch.start();

	  // 
	  // Build active t_slice list
	  //
	  std::vector<bool> active_t_slices(Lt);

	  if (1)
	  {
	    // Initialize the active time slices
	    for(int t=0; t < Lt; ++t)
	    {
	      active_t_slices[t] = false;
	    }

	    // Forward
	    for(int dt=0; dt < sink_source->Nt_forward; ++dt)
	    {
	      int t = (t_source + dt) % Lt;
	      active_t_slices[t] = true;
	    }

	    // Backward
	    for(int dt=0; dt < sink_source->Nt_backward; ++dt)
	    {
	      int t = (t_source - dt + 2*Lt) % Lt;
	      active_t_slices[t] = true;
	    }
	  }


	  //
	  // Look through sources, do the funny stuff through each soln, and stream the sink vectors past them
	  //
	  for(int colorvec_src=0; colorvec_src < srce_num_vecs; ++colorvec_src)
	  {
	    QDPIO::cout << "SOURCE: colorvec_src = " << colorvec_src << std::endl;
	    
	    StopWatch snarss1;
	    snarss1.reset();
	    snarss1.start();

	    //
	    // Loop over each spin source and invert. 
	    // Use the same color vector source. No spin dilution will be used.
	    //
	    const LatticeColorVectorSpinMatrix& soln_srce = prop_cache.getSoln(t_source, colorvec_src);

	    //
	    // Cache holding original solution vectors including displacements/derivatives
	    //
	    DispSolnCache disp_soln_cache(u_smr, soln_srce);
	    
	    //
	    // Loop over insertions for this source solutiuon vector
	    // Multiply by insertions
	    // Stream the sink solution vectors past it, and contract
	    //
	    QDPIO::cout << "SOURCE: do insertions for colorvec_src= " << colorvec_src << std::endl;
	    snarss1.reset();
	    snarss1.start();

	    for(auto dd = disp_gamma_moms.begin(); dd != disp_gamma_moms.end(); ++dd)
	    {
	      for(auto gg = dd->second.begin(); gg != dd->second.end(); ++gg)
	      {
		for(auto mm = gg->second.begin(); mm != gg->second.end(); ++mm)
		{
		  StopWatch snarss2;
		  snarss2.reset();
		  snarss2.start();

		  auto disp    = dd->first.deriv;
		  auto gamma   = gg->first.gamma;
		  auto mom     = mm->first;
		  int  mom_num = phases.momToNum(mom);

		  QDPIO::cout << "Insertion: disp= " << disp << "  gamma= " << gamma << " mom= " << mom << std::endl;

		  //
		  // Finally, the actual insertion
		  // NOTE: if we did not have the possibility for derivatives, then the displacement,
		  // and multiply by gamma could be moved outside the mom loop.
		  // The deriv/disp are cached, so the only cost is a lookup.
		  // The gamma is really just a rearrangement of spin components, but it does cost
		  // a traversal of the lattice
		  // The phases mult is a straight up cost.
		  //
		  // NOTE: ultimately, we are using gamma5 hermiticity to change the propagator from source at time
		  // t_sink to, instead, the t_slice
		  //
		  // Also NOTE: the gamma5 is hermitian. We will put the gamma_5 into the insertion, but since the need for the G5
		  // is a part of the sink solution vector, we will multiply here.
		  //
		  // NOTE: with suitable signs, the gamma_5 could be merged with the gamma
		  LatticeColorVectorSpinMatrix tmp = phases[mom_num] *
		    (Gamma(g5) * (Gamma(gamma) * disp_soln_cache.getDispVector(params.param.contract.use_derivP,
									       mom,
									       disp)));
	      
		  // Keys and stuff
		  SerialDBKey<KeyUnsmearedMesonElementalOperator_t>  key;
		  SerialDBData<ValUnsmearedMesonElementalOperator_t> val;

		  // The keys for the spin and displacements for this particular elemental operator
		  // No displacement for left colorvector, only displace right colorvector
		  // Invert the time - make it an independent key
		  multi1d<KeyValUnsmearedMesonElementalOperator_t> buf(phases.numSubsets());
		  for(int t=0; t < phases.numSubsets(); ++t)
		  {
		    if (! active_t_slices[t]) {continue;}
		
		    buf[t].key.key().derivP        = params.param.contract.use_derivP;
		    buf[t].key.key().t_sink        = t_sink;
		    buf[t].key.key().t_slice       = t;
		    buf[t].key.key().t_source      = t_source;
		    buf[t].key.key().colorvec_src  = colorvec_src;
		    buf[t].key.key().gamma         = gamma;
		    buf[t].key.key().displacement  = disp;
		    buf[t].key.key().mom           = mom;
		    buf[t].key.key().mass          = params.param.contract.mass_label;
		    buf[t].val.data().op.resize(sink_num_vecs,Ns,Ns);
		  }

		  // Stream the sink vectors past the insertion
		  // Will save a column of the genprop - corresponding to the current colorvec_src
		  for(int colorvec_snk=0; colorvec_snk < sink_num_vecs; ++colorvec_snk)
		  {
		    //QDPIO::cout << "stream: colorvec_snk = " << colorvec_snk << std::endl;
		    
		    // Slow fourier-transform
		    multi1d<SpinMatrixD> fred = sumMulti(localColorInnerProduct(prop_cache.getSoln(t_sink, colorvec_snk), tmp), phases.getSet());

		    for(int t=0; t < phases.numSubsets(); ++t)
		    {
		      if (! active_t_slices[t]) {continue;}

		      // Complete the gamma_5 hermitian adj on the sink by tacking on the gamma_5
		      SpinMatrixD gred = Gamma(g5) * fred[t];
			
		      for (int spin_snk = 0; spin_snk < Ns; ++spin_snk)
			for (int spin_src = 0; spin_src < Ns; ++spin_src)
			  buf[t].val.data().op(colorvec_snk,spin_snk,spin_src) = peekSpin(gred, spin_snk, spin_src);
		    }
		  }

		  // Insert these elementals into the db
		  for(int t=0; t < phases.numSubsets(); ++t)
		  {
		    if (! active_t_slices[t]) {continue;}

		    //QDPIO::cout << "insert key= " << buf[t].key.key() << std::endl;
		    //write(xml_out, "Insertion", buf[t].key.key());

		    qdp_db.insert(buf[t].key, buf[t].val);
		  }

		  snarss2.stop(); 
		  QDPIO::cout << " Time to build elemental: colorvec_src= " << colorvec_src
			      << "  gamma= " << gamma
			      << "  disp= " << disp
			      << "  mom= " << mom
			      << "  time = " << snarss2.getTimeInSeconds() << " secs " <<std::endl;
		} // mm
	      } // gg
	    } // dd

	    snarss1.stop(); 
	    QDPIO::cout << " Time to do all insertions: time = " << snarss1.getTimeInSeconds() << " secs " <<std::endl;
	  } // for colorvec_src
	} // for sink_source

	swatch.stop(); 
	QDPIO::cout << " SOURCE: time to compute all source solution vectors and insertions= " << swatch.getTimeInSeconds() << " secs" <<std::endl;
      }
      catch (const std::string& e) 
      {
	QDPIO::cout << name << ": caught exception around qprop: " << e << std::endl;
	QDP_abort(1);
      }

      // Close db
      qdp_db.close();

      // Close the xml output file
      pop(xml_out);     // UnsmearedHadronNode

      snoop.stop();
      QDPIO::cout << name << ": total time = " 
		  << snoop.getTimeInSeconds() 
		  << " secs" << std::endl;

      QDPIO::cout << name << ": ran successfully" << std::endl;

      END_CODE();
    } // func
  } // namespace InlineUnsmearedHadronNodeDistillationEnv

  /*! @} */  // end of group hadron

} // namespace Chroma

#endif
