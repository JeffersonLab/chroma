/*! \file
 * \brief Compute the propagator from distillution
 *
 * Propagator calculation in distillution
 */

#include "fermact.h"
#include "meas/inline/hadron/inline_prop_distillution_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "util/ferm/distillution_noise.h"
#include "util/ferm/key_timeslice_colorvec.h"
#include "qdp_map_obj.h"
#include "qdp_map_obj_disk.h"
#include "qdp_disk_map_slice.h"
#include "util/ferm/key_peram_dist.h"
#include "util/ferm/key_prop_distillution.h"
#include "util/ferm/transf.h"
#include "util/ferm/spin_rep.h"
#include "util/ferm/diractodr.h"
#include "util/ferm/twoquark_contract_ops.h"
#include "util/ft/time_slice_set.h"
#include "util/info/proginfo.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"
#include "meas/inline/make_xml_file.h"

#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{ 
  //----------------------------------------------------------------------------
  //! Abstract type for quarkline construction
  class AbsQuarkLine
  {
  public:
    //! Virtual destructor
    virtual ~AbsQuarkLine() {}

    //! Get a source
    virtual LatticeColorVector getSrc(int t_source, int dist_src) const = 0;

    //! Get number of vectors
    virtual int getNumVecs() const = 0;

    //! Get number of space dilutions
    virtual int getNumSpaceDils() const = 0;

    //! Get number of time dilutions
    virtual int getNumTimeDils() const = 0;

    //! Get quark line number
    virtual int getQuarkLine() const = 0;

    //! Get the time sources
    virtual std::vector<int> getTimeSources() const = 0;

    //! Get source keys
    virtual std::list<KeyPropDist_t> getSrcKeys(int t_source, int dist_src) const = 0;

    //! Get sink keys
    virtual std::list<KeyPropDist_t> getSnkKeys(int t_source, int dist_src) const = 0;

    //! Get perambulator keys
    virtual std::list<KeyPeramDist_t> getPeramKeys(int t_source) const = 0;
  };


  //----------------------------------------------------------------------------
  typedef QDP::MapObjectDisk< KeyTimeSliceColorVec_t,TimeSliceIO<LatticeColorVector> > MOD_t;

  //----------------------------------------------------------------------------
  //! Quark line factory (foundry)
  typedef SingletonHolder< 
    ObjectFactory<AbsQuarkLine,
		  std::string,
		  TYPELIST_7(XMLReader&, const std::string&,
			     const DistillutionNoise&, 
			     MOD_t&,
			     const TimeSliceSet&,
			     int,
			     const std::string&),
		  AbsQuarkLine* (*)(XMLReader&,
				    const std::string&,
				    const DistillutionNoise&, 
				    MOD_t&,
				    const TimeSliceSet&,
				    int,
				    const std::string&), StringFactoryError> >
  TheQuarkLineFactory;



  //----------------------------------------------------------------------------
  //! Utilities
  namespace PropDistillutionUtilEnv
  {
    // More utilities
    namespace
    {
      //! Check space dilutions
      int checkSpaceDils(int num_space_dils, int num_vecs)
      {
	int num = num_space_dils;

	// Reset/barf if bogus
	if (num_space_dils == 0)
	{
	  num = num_vecs;
	}
	else if (num_space_dils < 0 || num_space_dils > num_vecs)
	{
	  QDPIO::cerr << __func__ << ": invalid size of num_space_dils = " << num_space_dils << std::endl;
	  QDP_abort(1);
	}
	else if ((num_vecs % num_space_dils) != 0)
	{
	  QDPIO::cerr << __func__ 
		      << ": num_space_dils = " << num_space_dils 
		      << "  not a divisor of num_vecs = " << num_vecs
		      << std::endl;
	  QDP_abort(1);
	}

	return num;
      }


      //! Check time dilutions
      int checkTimeDils(int num_time_dils, int Lt)
      {
	int num = num_time_dils;

	// Reset/barf if bogus
	if (num_time_dils == 0)
	{
	  num = Lt;
	}
	else if (num_time_dils < 0 || num_time_dils > Lt)
	{
	  QDPIO::cerr << __func__ << ": invalid size of num_time_dils = " << num_time_dils << std::endl;
	  QDP_abort(1);
	}
	else if ((Lt % num_time_dils) != 0)
	{
	  QDPIO::cerr << __func__ 
		      << ": num_time_dils = " << num_time_dils 
		      << "  not a divisor of the time extent = " << Lt 
			  << std::endl;
	  QDP_abort(1);
	}

	return num;
      }
    } // anonymous namespace


    //----------------------------------------------------------------------------
    // The noise for this quark line.
    // NOTE: the noise is fixed for a type of source, but given for all time-slices
    // 
    multi2d<Complex> generateNoise(const DistillutionNoise& dist_noise_obj, 
				   int num_vecs, int quark_line, bool annih, const std::string& mass)
    {
      DistQuarkLines_t info;
      info.num_vecs   = num_vecs;
      info.quark_line = quark_line;
      info.annih      = false;
      info.mass       = mass;

      multi2d<Complex> eta(dist_noise_obj.getRNG(info));

#if 0
      // Just for debugging
      for(int i=0; i < eta.nrows(); ++i)
	for(int j=0; j < eta.ncols(); ++j)
	{
	  float re = toFloat(real(eta(i,j)));
	  float im = toFloat(imag(eta(i,j)));

	  QDPIO::cout << name << ": line= " << info.quark_line 
		      << "  eta("<<i<<","<<j<<")= ( "
		      << (fabs(re)<0.1?0:re) << " , " 
		      << (fabs(im)<0.1?0:im) << " )\n";
	}
#endif

      return eta;
    }

	  
    //----------------------------------------------------------------------------
    //----------------------------------------------------------------------------
    //----------------------------------------------------------------------------
    namespace Connected
    {
      //! Parameter structure
      struct Params 
      {
	Params();
	Params(XMLReader& xml_in, const std::string& path);
	void writeXML(XMLWriter& xml_out, const std::string& path) const;
	
	int               num_vecs;         /*!< Number of vectors */
	int               num_space_dils;   /*!< Number of eigenvector dilutions to use */

	std::vector<int>  t_sources;        /*!< Time sources */ 
	int               Nt_forward;       /*!< Time-slices in the forward direction */
	int               Nt_backward;      /*!< Time-slices in the backward direction */
      };


      //! Propagator input
      Params::Params(XMLReader& xml, const string& path)
      {
	XMLReader inputtop(xml, path);

	read(inputtop, "num_vecs", num_vecs);
	read(inputtop, "num_space_dils", num_space_dils);
	read(inputtop, "t_sources", t_sources);
	read(inputtop, "Nt_forward", Nt_forward);
	read(inputtop, "Nt_backward", Nt_backward);
      }

      //! Propagator output
      void Params::writeXML(XMLWriter& xml, const string& path) const
      {
	push(xml, path);

	write(xml, "num_vecs", num_vecs);
	write(xml, "num_space_dils", num_space_dils);
	write(xml, "t_sources", t_sources);
	write(xml, "Nt_forward", Nt_forward);
	write(xml, "Nt_backward", Nt_backward);

	pop(xml);
      }

      //! Propagator input
      void read(XMLReader& xml, const string& path, Params& input)
      {
	Params tmp(xml, path);
	input = tmp;
      }

      //! Propagator output
      void write(XMLWriter& xml, const string& path, const Params& input)
      {
	input.writeXML(xml, path);
      }


      //----------------------------------------------------------------------------
      //! Connected quark lines
      /*!< 
       * Pull out a time-slice of the color vector source, and add it in a crystal fashion
       * with other vectors
       */
      class QuarkLineFact : public AbsQuarkLine
      {
      public:
	QuarkLineFact(const Params& params_,
		      const DistillutionNoise& dist_noise_obj_, 
		      MOD_t& eigen_source_,
		      const TimeSliceSet& time_slice_set_,
		      int quark_line_,
		      const std::string& mass_)
	: params(params_), 
	  dist_noise_obj(dist_noise_obj_), eigen_source(eigen_source_), time_slice_set(time_slice_set_), 
	  quark_line(quark_line_), mass(mass_)
	  {
	    // Initialize the noise for this quark line
	    eta = generateNoise(dist_noise_obj, params.num_vecs, quark_line, false, mass);

	    // Reset/barf if bogus
	    params.num_space_dils = checkSpaceDils(params.num_space_dils, params.num_vecs);

#if 0
	    // Another sanity check
	    if (params.num_vecs > eigen_source.size())
	    {
	      QDPIO::cerr << __func__ << ": num_vecs= " << params.num_vecs
			  << " is greater than the number of available colorvectors= "
			  << eigen_source.size() << endl;
	      QDP_abort(1);
	    }
#endif
	  }

	//! Get a source
	virtual LatticeColorVector getSrc(int t_source, int dist_src) const;

	//! Get number of vectors
	virtual int getNumVecs() const {return params.num_vecs;}

	//! Get number of space dilutions
	virtual int getNumSpaceDils() const {return params.num_space_dils;}

	//! Get number of time dilutions
	virtual int getNumTimeDils() const {return Layout::lattSize()[dist_noise_obj.getDecayDir()];}

	//! Get quark line number
	virtual int getQuarkLine() const {return quark_line;}

	//! Get the time sources
	virtual std::vector<int> getTimeSources() const {return params.t_sources;}

	//! Get source keys
	virtual std::list<KeyPropDist_t> getSrcKeys(int t_source, int dist_src) const;

	//! Get sink keys
	virtual std::list<KeyPropDist_t> getSnkKeys(int t_source, int dist_src) const;

	//! Get perambulator keys
	virtual std::list<KeyPeramDist_t> getPeramKeys(int t_source) const;

      private:
	std::vector<bool> getActiveTSlices(int t_source) const;

      private:
	// Arguments
	Params                    params;
	const DistillutionNoise&  dist_noise_obj;
	MOD_t&                    eigen_source;
	const TimeSliceSet&       time_slice_set;
	int                       quark_line;
	std::string               mass;

	// Local
	multi2d<Complex>          eta;
      };



      //----------------------------------------------------------------------------
      //! Prepare a distilluted source
      LatticeColorVector QuarkLineFact::getSrc(int t_source, int dist_src) const
      {
	LatticeColorVector vec_srce = zero;

	for(int colorvec_source=dist_src; colorvec_source < params.num_vecs; colorvec_source += params.num_space_dils)
	{
	  QDPIO::cout << "colorvec_source = " << colorvec_source << endl;

	  // Get the actual time slice
	  int t_actual = dist_noise_obj.getTime(t_source);

	  KeyTimeSliceColorVec_t key_vec;
	  key_vec.t_slice = t_actual;
	  key_vec.colorvec = colorvec_source;

	  LatticeColorVector tmpvec = zero;
	  TimeSliceIO<LatticeColorVector> time_slice_io(tmpvec, t_actual);

	  eigen_source.get(key_vec, time_slice_io);

	  vec_srce[time_slice_set.getSet()[t_actual]] += eta(t_source, colorvec_source) * tmpvec;
	}

	return vec_srce;
      }


      //----------------------------------------------------------------------------
      //! Get active time-slices
      std::vector<bool> QuarkLineFact::getActiveTSlices(int t_source) const
      {
	// Initialize the active time slices
	const int decay_dir = dist_noise_obj.getDecayDir();
	const int Lt = Layout::lattSize()[decay_dir];

	std::vector<bool> active_t_slices(Lt);
	for(int t=0; t < Lt; ++t)
	{
	  active_t_slices[t] = false;
	}

	// Forward
	for(int dt=0; dt < params.Nt_forward; ++dt)
	{
	  int t = t_source + dt;
	  active_t_slices[t % Lt] = true;
	}

	// Backward
	for(int dt=0; dt < params.Nt_backward; ++dt)
	{
	  int t = t_source - dt;
	  while (t < 0) {t += Lt;} 

	  active_t_slices[t % Lt] = true;
	}

	return active_t_slices;
      }


      //----------------------------------------------------------------------------
      //! Get source keys
      std::list<KeyPropDist_t> QuarkLineFact::getSrcKeys(int t_source, int dist_src) const
      {
	std::list<KeyPropDist_t> keys;

	KeyPropDist_t key;

	key.prop_type    = "SRC";
	key.annihP       = false;
	key.t_source     = t_source;
	key.t_slice      = t_source;
	key.dist_src     = dist_src;
	key.spin_src     = -1;
	key.spin_snk     = -1;
	key.quark_line   = quark_line;
	key.mass         = mass;

	keys.push_back(key);

	return keys;
      }

	
      //----------------------------------------------------------------------------
      //! Get sink keys
      std::list<KeyPropDist_t> QuarkLineFact::getSnkKeys(int t_source, int dist_src) const
      {
	std::list<KeyPropDist_t> keys;

	std::vector<bool> active_t_slices = getActiveTSlices(t_source);
	
	const int Lt = Layout::lattSize()[dist_noise_obj.getDecayDir()];

	for(int spin_source=0; spin_source < Ns; ++spin_source)
	{
	  for(int spin_sink=0; spin_sink < Ns; ++spin_sink)
	  {
	    for(int t=0; t < Lt; ++t)
	    {
	      if (! active_t_slices[t]) {continue;}

	      KeyPropDist_t key;

	      key.prop_type    = "SNK";
	      key.annihP       = false;
	      key.t_source     = t_source;
	      key.t_slice      = t;
	      key.dist_src     = dist_src;
	      key.spin_src     = spin_source;
	      key.spin_snk     = spin_sink;
	      key.quark_line   = quark_line;
	      key.mass         = mass;

//            QDPIO::cout << key << std::flush;

	      keys.push_back(key);
	    } // for t
	  } // for spin_sink
	} // for spin_source

	return keys;
      }
      
      //----------------------------------------------------------------------------
      //! Get perambulator keys
      std::list<KeyPeramDist_t> QuarkLineFact::getPeramKeys(int t_source) const
      {
	std::list<KeyPeramDist_t> keys;

	std::vector<bool> active_t_slices = getActiveTSlices(t_source);
	
	const int Lt = Layout::lattSize()[dist_noise_obj.getDecayDir()];

	for(int spin_source=0; spin_source < Ns; ++spin_source)
	{
	  for(int spin_sink=0; spin_sink < Ns; ++spin_sink)
	  {
	    for(int t=0; t < Lt; ++t)
	    {
	      if (! active_t_slices[t]) {continue;}

	      KeyPeramDist_t key;

	      key.quark_line   = quark_line;
	      key.annihP       = false;
	      key.t_slice      = t;
	      key.t_source     = t_source;
	      key.spin_src     = spin_source;
	      key.spin_snk     = spin_sink;
	      key.mass         = mass;

//            QDPIO::cout << key << std::flush;

	      keys.push_back(key);
	    } // for t
	  } // for spin_sink
	} // for spin_source

	return keys;
      }

	
    } // namespace Connected




    //----------------------------------------------------------------------------
    //----------------------------------------------------------------------------
    //----------------------------------------------------------------------------
    namespace Annihilation
    {
      //! Parameter structure
      struct Params 
      {
	Params();
	Params(XMLReader& xml_in, const std::string& path);
	void writeXML(XMLWriter& xml_out, const std::string& path) const;
	
	int               num_vecs;         /*!< Number of vectors */
	int               num_space_dils;   /*!< Number of eigenvector dilutions to use in space direction */
	int               num_time_dils;    /*!< Number of eigenvector dilutions to use in time direction */
      };


      //! Propagator input
      Params::Params(XMLReader& xml, const string& path)
      {
	XMLReader inputtop(xml, path);

	read(inputtop, "num_vecs", num_vecs);
	read(inputtop, "num_space_dils", num_space_dils);
	read(inputtop, "num_time_dils", num_time_dils);
      }

      //! Propagator output
      void Params::writeXML(XMLWriter& xml, const string& path) const
      {
	push(xml, path);

	write(xml, "num_vecs", num_vecs);
	write(xml, "num_space_dils", num_space_dils);
	write(xml, "num_time_dils", num_time_dils);

	pop(xml);
      }

      //! Propagator input
      void read(XMLReader& xml, const string& path, Params& input)
      {
	Params tmp(xml, path);
	input = tmp;
      }

      //! Propagator output
      void write(XMLWriter& xml, const string& path, const Params& input)
      {
	input.writeXML(xml, path);
      }


      //----------------------------------------------------------------------------
      //! Annihilation quark lines
      /*!< 
       * Pull out a time-slice of the color vector source, and add it in a crystal fashion
       * with other vectors
       */
      class QuarkLineFact : public AbsQuarkLine
      {
      public:
	QuarkLineFact(const Params& params_,
		      const DistillutionNoise& dist_noise_obj_, 
		      MOD_t& eigen_source_,
		      const TimeSliceSet& time_slice_set_,
		      int quark_line_,
		      const std::string& mass_)
	: params(params_), 
	  dist_noise_obj(dist_noise_obj_), eigen_source(eigen_source_), time_slice_set(time_slice_set_), 
	  quark_line(quark_line_), mass(mass_)
	  {
	    // Initialize the noise for this quark line
	    eta = generateNoise(dist_noise_obj, params.num_vecs, quark_line, true, mass);

	    // Reset/barf if bogus
	    params.num_space_dils = checkSpaceDils(params.num_space_dils, params.num_vecs);
	    params.num_time_dils  = checkTimeDils(params.num_time_dils, Layout::lattSize()[dist_noise_obj.getDecayDir()]);

#if 0
	    // Another sanity check
	    if (params.num_vecs > eigen_source.size())
	    {
	      QDPIO::cerr << __func__ << ": num_vecs= " << params.num_vecs
			  << " is greater than the number of available colorvectors= "
			  << eigen_source.size() << endl;
	      QDP_abort(1);
	    }
#endif
	  }

	//! Get a source
	virtual LatticeColorVector getSrc(int t_source, int dist_src) const;

	//! Get number of vectors
	virtual int getNumVecs() const {return params.num_vecs;}

	//! Get number of space dilutions
	virtual int getNumSpaceDils() const {return params.num_space_dils;}

	//! Get number of time dilutions
	virtual int getNumTimeDils() const {return Layout::lattSize()[dist_noise_obj.getDecayDir()];}

	//! Get quark line number
	virtual int getQuarkLine() const {return quark_line;}

	//! Get the time sources
	virtual std::vector<int> getTimeSources() const;

	//! Get source keys
	virtual std::list<KeyPropDist_t> getSrcKeys(int t_source, int dist_src) const;

	//! Get sink keys
	virtual std::list<KeyPropDist_t> getSnkKeys(int t_source, int dist_src) const;

	//! Get perambulator keys
	virtual std::list<KeyPeramDist_t> getPeramKeys(int t_source) const;

      private:
	// Arguments
	Params                    params;
	const DistillutionNoise&  dist_noise_obj;
	MOD_t&                    eigen_source;
	const TimeSliceSet&       time_slice_set;
	int                       quark_line;
	std::string               mass;

	// Local
	multi2d<Complex>          eta;
      };



      //----------------------------------------------------------------------------
      //! Get the time sources
      std::vector<int> QuarkLineFact::getTimeSources() const
      {
	std::vector<int> t_sources;

	for(int t=0; t < params.num_time_dils; ++t)
	{
	  t_sources.push_back(t);
	}

	return t_sources;
      }


      //----------------------------------------------------------------------------
      //! Prepare a distilluted source
      LatticeColorVector QuarkLineFact::getSrc(int t_src, int dist_src) const
      {
	LatticeColorVector vec_srce = zero;

	const int Lt = Layout::lattSize()[dist_noise_obj.getDecayDir()];

	for(int time_source=t_src; time_source < Lt; time_source += params.num_time_dils)
	{
	  for(int colorvec_source=dist_src; colorvec_source < params.num_vecs; colorvec_source += params.num_space_dils)
	  {
	    // Get the actual time slice
	    int t_actual = dist_noise_obj.getTime(time_source);

	    KeyTimeSliceColorVec_t key_vec;
	    key_vec.t_slice = t_actual;
	    key_vec.colorvec = colorvec_source;

	    LatticeColorVector tmpvec = zero;
	    TimeSliceIO<LatticeColorVector> time_slice_io(tmpvec, t_actual);

	    eigen_source.get(key_vec, time_slice_io);

	    vec_srce[time_slice_set.getSet()[t_actual]] += eta(time_source, colorvec_source) * tmpvec;
	  }
	} // for time_source

	return vec_srce;
      }


      //----------------------------------------------------------------------------
      //! Get source keys
      std::list<KeyPropDist_t> QuarkLineFact::getSrcKeys(int t_source, int dist_src) const
      {
	std::list<KeyPropDist_t> keys;

	const int Lt = Layout::lattSize()[dist_noise_obj.getDecayDir()];

	for(int t=t_source; t < Lt; t += params.num_time_dils)
	{
	  KeyPropDist_t key;

	  key.prop_type    = "SRC";
	  key.annihP       = true;
	  key.t_source     = t;
	  key.t_slice      = t;
	  key.dist_src     = dist_src;
	  key.spin_src     = -1;
	  key.spin_snk     = -1;
	  key.quark_line   = quark_line;
	  key.mass         = mass;

	  keys.push_back(key);
	}

	return keys;
      }

	
      //----------------------------------------------------------------------------
      //! Get sink keys
      std::list<KeyPropDist_t> QuarkLineFact::getSnkKeys(int t_source, int dist_src) const
      {
	std::list<KeyPropDist_t> keys;

	const int Lt = Layout::lattSize()[dist_noise_obj.getDecayDir()];

	for(int spin_source=0; spin_source < Ns; ++spin_source)
	{
	  for(int spin_sink=0; spin_sink < Ns; ++spin_sink)
	  {
	    for(int t=t_source; t < Lt; t += params.num_time_dils)
	    {
	      KeyPropDist_t key;

	      key.prop_type    = "SNK";
	      key.annihP       = true;
	      key.t_source     = t;
	      key.t_slice      = t;
	      key.dist_src     = dist_src;
	      key.spin_src     = spin_source;
	      key.spin_snk     = spin_sink;
	      key.quark_line   = quark_line;
	      key.mass         = mass;

//            QDPIO::cout << key << std::flush;

	      keys.push_back(key);
	    } // for t
	  } // for spin_sink
	} // for spin_source

	return keys;
      }
      
      //----------------------------------------------------------------------------
      //! Get perambulator keys
      std::list<KeyPeramDist_t> QuarkLineFact::getPeramKeys(int t_source) const
      {
	std::list<KeyPeramDist_t> keys;

	const int Lt = Layout::lattSize()[dist_noise_obj.getDecayDir()];

	for(int spin_source=0; spin_source < Ns; ++spin_source)
	{
	  for(int spin_sink=0; spin_sink < Ns; ++spin_sink)
	  {
	    for(int t=t_source; t < Lt; t += params.num_time_dils)
	    {
	      KeyPeramDist_t key;

	      key.quark_line   = quark_line;
	      key.annihP       = true;
	      key.t_slice      = t;
	      key.t_source     = t;
	      key.spin_src     = spin_source;
	      key.spin_snk     = spin_sink;
	      key.mass         = mass;

//            QDPIO::cout << key << std::flush;

	      keys.push_back(key);
	    } // for t
	  } // for spin_sink
	} // for spin_source

	return keys;
      }

	
    } // namespace Annihilation


    //----------------------------------------------------------------------------
    namespace
    {
      AbsQuarkLine* createConn(XMLReader& xml_in, 
			       const std::string& path,
			       const DistillutionNoise& dist_noise_obj, 
			       MOD_t& eigen_source,
			       const TimeSliceSet& time_slice_set,
			       int quark_line,
			       const std::string& mass)
      {
	return new Connected::QuarkLineFact(Connected::Params(xml_in, path),
					    dist_noise_obj, eigen_source, time_slice_set, 
					    quark_line, mass);
      }

      AbsQuarkLine* createAnnih(XMLReader& xml_in, 
				const std::string& path,
				const DistillutionNoise& dist_noise_obj, 
				MOD_t& eigen_source,
				const TimeSliceSet& time_slice_set,
				int quark_line,
			       const std::string& mass)
      {
	return new Annihilation::QuarkLineFact(Annihilation::Params(xml_in, path),
					       dist_noise_obj, eigen_source, time_slice_set, 
					       quark_line, mass);
      }

      //! Local registration flag
      bool registered = false;
    }

    //----------------------------------------------------------------------------
    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= TheQuarkLineFactory::Instance().registerObject(std::string("CONN"), createConn);
	success &= TheQuarkLineFactory::Instance().registerObject(std::string("ANNIH"), createAnnih);
	registered = true;
      }
      return success;
    }
    
  } // anonymous namespace



  //----------------------------------------------------------------------------
  namespace InlinePropDistillutionEnv 
  {
    //! Propagator input
    void read(XMLReader& xml, const string& path, Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "save_peramP", input.save_peramP);
      read(inputtop, "save_srcP", input.save_srcP);
      read(inputtop, "save_solnP", input.save_solnP);
      read(inputtop, "gauge_id", input.gauge_id);
      read(inputtop, "distillution_id", input.distillution_id);
      read(inputtop, "colorvec_file", input.colorvec_file);
      read(inputtop, "soln_file", input.soln_file);
      read(inputtop, "peram_file", input.peram_file);
    }

    //! Propagator output
    void write(XMLWriter& xml, const string& path, const Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "save_peramP", input.save_peramP);
      write(xml, "save_srcP", input.save_srcP);
      write(xml, "save_solnP", input.save_solnP);
      write(xml, "gauge_id", input.gauge_id);
      write(xml, "distillution_id", input.distillution_id);
      write(xml, "colorvec_file", input.colorvec_file);
      write(xml, "soln_file", input.soln_file);
      write(xml, "peram_file", input.peram_file);

      pop(xml);
    }


    //! Propagator input
    void read(XMLReader& xml, const string& path, Params::Param_t::Contract_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "quark_lines", input.quark_lines);
      read(inputtop, "mass", input.mass);

      // Read quark-line parameters
      input.quark_line_xml = readXMLGroup(inputtop, "QuarkLine", "QuarkLineType");
    }

    //! Propagator output
    void write(XMLWriter& xml, const string& path, const Params::Param_t::Contract_t& input)
    {
      push(xml, path);

      write(xml, "quark_lines", input.quark_lines);
      write(xml, "mass", input.mass);
      xml << input.quark_line_xml.xml;

      pop(xml);
    }


    //! Propagator input
    void read(XMLReader& xml, const string& path, Params::Param_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "Propagator", input.prop);
      read(inputtop, "Contractions", input.contract);
    }

    //! Propagator output
    void write(XMLWriter& xml, const string& path, const Params::Param_t& input)
    {
      push(xml, path);

      write(xml, "Propagator", input.prop);
      write(xml, "Contractions", input.contract);

      pop(xml);
    }


    //! Propagator input
    void read(XMLReader& xml, const string& path, Params& input)
    {
      Params tmp(xml, path);
      input = tmp;
    }

    //! Propagator output
    void write(XMLWriter& xml, const string& path, const Params& input)
    {
      push(xml, path);
    
      write(xml, "Param", input.param);
      write(xml, "NamedObject", input.named_obj);

      pop(xml);
    }
  } // namespace InlinePropDistillutionEnv 


  //----------------------------------------------------------------------------
  namespace InlinePropDistillutionEnv 
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
      
    const std::string name = "PROP_DISTILLUTION";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= PropDistillutionUtilEnv::registerAll();
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
	QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << endl;
	QDP_abort(1);
      }
    }



    // Function call
    void 
    InlineMeas::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
      // If xml file not empty, then use alternate
      if (params.xml_file != "")
      {
	string xml_file = makeXMLFileName(params.xml_file, update_no);

	push(xml_out, "PropDistillution");
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
	QDPIO::cerr << name << ": caught dynamic cast error" << endl;
	QDP_abort(1);
      }
      catch (const string& e) 
      {
	QDPIO::cerr << name << ": map call failed: " << e << endl;
	QDP_abort(1);
      }

      push(xml_out, "PropDistillution");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << name << ": propagator calculation" << endl;

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

      //
      // Read in the source along with relevant information.
      // 
      QDPIO::cout << "Snarf the distillution factory from a named buffer" << endl;
      try
      {
	TheNamedObjMap::Instance().getData< Handle< DistillutionNoise > >(params.named_obj.distillution_id);
      }    
      catch (std::bad_cast) {
	QDPIO::cerr << name << ": caught dynamic cast error" << endl;
	QDP_abort(1);
      }
      catch (const string& e) {
	QDPIO::cerr << name << ": error extracting source_header: " << e << endl;
	QDP_abort(1);
      }
      catch( const char* e) {
	QDPIO::cerr << name <<": Caught some char* exception:" << endl;
	QDPIO::cerr << e << endl;
	QDPIO::cerr << "Rethrowing" << endl;
	throw;
      }

      // Cast should be valid now
      const DistillutionNoise& dist_noise_obj =
	*(TheNamedObjMap::Instance().getData< Handle<DistillutionNoise> >(params.named_obj.distillution_id));

      // Some diagnostics
      QDPIO::cout << "Distillution factory: ensemble= XX" << dist_noise_obj.getEnsemble() << "XX  "
		  << "sequence= XX" << dist_noise_obj.getSequence() << "XX\n"
		  << "t_origin= " << dist_noise_obj.getOrigin() << std::endl;

      // Will use TimeSliceSet-s a lot
      const int decay_dir = dist_noise_obj.getDecayDir();
      const int Lt        = Layout::lattSize()[decay_dir];

      // A sanity check
      if (decay_dir != Nd-1)
      {
	QDPIO::cerr << __func__ << ": TimeSliceIO only supports decay_dir= " << Nd-1 << "\n";
	QDP_abort(1);
      }

      // The time-slice set
      TimeSliceSet time_slice_set(decay_dir);


      //
      // Read in the source along with relevant information.
      // 
      QDPIO::cout << "Snarf the source from a map object disk file" << endl;

      QDP::MapObjectDisk<KeyTimeSliceColorVec_t,TimeSliceIO<LatticeColorVector> > eigen_source;
      eigen_source.setDebug(0);

      try
      {
	// Open
	QDPIO::cout << "Open file= " << params.named_obj.colorvec_file << endl;
	eigen_source.open(params.named_obj.colorvec_file);

	// Snarf the source info. 
	string user_str;
	QDPIO::cout << "Get user data" << endl;
	eigen_source.getUserdata(user_str);

	// Write it
	QDPIO::cout << "Write to an xml file" << endl;
	XMLBufferWriter xml_buf(user_str);
	write(xml_out, "Source_info", xml_buf);
      }    
      catch (std::bad_cast) {
	QDPIO::cerr << name << ": caught dynamic cast error" << endl;
	QDP_abort(1);
      }
      catch (const string& e) {
	QDPIO::cerr << name << ": error extracting source_header: " << e << endl;
	QDP_abort(1);
      }
      catch( const char* e) {
	QDPIO::cerr << name <<": Caught some char* exception:" << endl;
	QDPIO::cerr << e << endl;
	QDPIO::cerr << "Rethrowing" << endl;
	throw;
      }

      QDPIO::cout << "Set source mod file to read" << endl;

      QDPIO::cout << "Source successfully read and parsed" << endl;


      // Sanity check - write out the norm2 of the source in the decay_dir direction
      // Use this for any possible verification
      {
	QDPIO::cout << "Lookup source 0" << endl;
	LatticeColorVector tmpvec;

	for(int t=0; t < Lt; ++t)
	{
	  KeyTimeSliceColorVec_t key_vec;
	  key_vec.t_slice = t;
	  key_vec.colorvec = 0;

	  TimeSliceIO<LatticeColorVector> time_slice_io(tmpvec, t);
	  eigen_source.get(key_vec, time_slice_io);
	}

	multi1d<Double> source_corrs = sumMulti(localNorm2(tmpvec), time_slice_set.getSet());

	push(xml_out, "Source_correlators");
	write(xml_out, "source_corrs", source_corrs);
	pop(xml_out);
      }

      QDPIO::cout << "Source studied: do some other sanity checks" << endl;

      //
      // Map-object-disk storage
      //
      QDP::MapObjectDisk<KeyPropDist_t, TimeSliceIO<LatticeColorVector> > prop_obj;
      prop_obj.setDebug(0);

      if (params.named_obj.save_srcP || params.named_obj.save_solnP)
      {
	QDPIO::cout << "Open solution file" << endl;

	if (! prop_obj.fileExists(params.named_obj.soln_file))
	{
	  XMLBufferWriter file_xml;

	  push(file_xml, "MODMetaData");
	  write(file_xml, "id", string("propDist"));
	  write(file_xml, "lattSize", QDP::Layout::lattSize());
	  write(file_xml, "decay_dir", decay_dir);
	  write(file_xml, "ensemble", dist_noise_obj.getEnsemble());
	  write(file_xml, "sequence", dist_noise_obj.getSequence());
	  write(file_xml, "t_origin", dist_noise_obj.getOrigin());
	  file_xml << params.param.contract.quark_line_xml.xml;
	  write(file_xml, "quark_line", params.param.contract.quark_lines[0]);
	  proginfo(file_xml);    // Print out basic program info
	  write(file_xml, "Params", params.param);
	  write(file_xml, "Config_info", gauge_xml);
	  pop(file_xml);

	  std::string file_str(file_xml.str());
	  
	  prop_obj.insertUserdata(file_xml.str());
	  prop_obj.open(params.named_obj.soln_file, std::ios_base::in | std::ios_base::out | std::ios_base::trunc);
	}
	else
	{
	  prop_obj.open(params.named_obj.soln_file);
	}

	QDPIO::cout << "Finished opening solution file" << endl;
      }


      //
      // DB storage
      //
      QDP::MapObjectDisk<KeyPeramDist_t, ValPeramDist_t> qdp_db;
      qdp_db.setDebug(0);

      if (params.named_obj.save_peramP)
      {
	QDPIO::cout << "Open peram file" << endl;

	if (! qdp_db.fileExists(params.named_obj.peram_file))
	{
	  XMLBufferWriter file_xml;

	  push(file_xml, "MODMetaData");
	  write(file_xml, "id", string("peramDist"));
	  write(file_xml, "lattSize", QDP::Layout::lattSize());
	  write(file_xml, "decay_dir", decay_dir);
	  proginfo(file_xml);    // Print out basic program info
//	  write(file_xml, "Weights", getEigenValues(eigen_source, eigen_source.size()));
	  write(file_xml, "ensemble", dist_noise_obj.getEnsemble());
	  write(file_xml, "sequence", dist_noise_obj.getSequence());
	  file_xml << params.param.contract.quark_line_xml.xml;
	  write(file_xml, "t_origin", dist_noise_obj.getOrigin());
	  write(file_xml, "quark_line", params.param.contract.quark_lines[0]);
	  proginfo(file_xml);    // Print out basic program info
	  write(file_xml, "Params", params.param);
	  write(file_xml, "Config_info", gauge_xml);
	  pop(file_xml);

	  std::string file_str(file_xml.str());
	  
	  qdp_db.insertUserdata(file_xml.str());
	  qdp_db.open(params.named_obj.peram_file, std::ios_base::in | std::ios_base::out | std::ios_base::trunc);
	}
	else
	{
	  qdp_db.open(params.named_obj.peram_file);
	}

	QDPIO::cout << "Finished opening peram file" << endl;
      }


      // Total number of iterations
      int ncg_had = 0;


      // Rotation from DR to DP
      SpinMatrix diracToDRMat(DiracToDRMat());
      std::vector<MatrixSpinRep_t> diracToDrMatPlus = convertTwoQuarkSpin(diracToDRMat);
      std::vector<MatrixSpinRep_t> diracToDrMatMinus = convertTwoQuarkSpin(adj(diracToDRMat));


      //
      // Try the factories
      //
      try
      {
	StopWatch swatch;
	swatch.reset();
	QDPIO::cout << "Try the various factories" << endl;

	// Typedefs to save typing
	typedef LatticeFermion               T;
	typedef multi1d<LatticeColorMatrix>  P;
	typedef multi1d<LatticeColorMatrix>  Q;

	//
	// Initialize fermion action
	//
	std::istringstream  xml_s(params.param.prop.fermact.xml);
	XMLReader  fermacttop(xml_s);
	QDPIO::cout << "FermAct = " << params.param.prop.fermact.id << endl;

	// Generic Wilson-Type stuff
	Handle< FermionAction<T,P,Q> >
	  S_f(TheFermionActionFactory::Instance().createObject(params.param.prop.fermact.id,
							       fermacttop,
							       params.param.prop.fermact.path));

	Handle< FermState<T,P,Q> > state(S_f->createState(u));

	Handle< SystemSolver<LatticeFermion> > PP = S_f->qprop(state,
							       params.param.prop.invParam);
      
	QDPIO::cout << "Suitable factory found: compute all the quark props" << endl;
	swatch.start();


	//
	// Loop over the source color and spin, creating the source
	// and calling the relevant propagator routines.
	//
	// Loop over each quark-line
	for(std::vector<int>::const_iterator quark_line= params.param.contract.quark_lines.begin();
	    quark_line != params.param.contract.quark_lines.end();
	    ++quark_line)
	{
	  QDPIO::cout << "quark_line = " << *quark_line << "  type= " << params.param.contract.quark_line_xml.id << endl; 

	  //
	  // Factory for quark line
	  //
	  Handle<AbsQuarkLine> quark_line_fact;

	  try 
	  {
	    std::istringstream  xml_l(params.param.contract.quark_line_xml.xml);
	    XMLReader  linktop(xml_l);
	    QDPIO::cout << "Quark line type = " << params.param.contract.quark_line_xml.id << endl;
	    QDPIO::cout << "Quark line xml = XX" << params.param.contract.quark_line_xml.xml << "XX" << endl;
	    
	    QDPIO::cout << "Create factory" << std::endl;

	    quark_line_fact = 
	      TheQuarkLineFactory::Instance().createObject(params.param.contract.quark_line_xml.id,
							   linktop,
							   params.param.contract.quark_line_xml.path,
							   dist_noise_obj, eigen_source, 
							   time_slice_set, 
							   *quark_line,
							   params.param.contract.mass);

	    QDPIO::cout << "Factory created" << std::endl;
	  }
	  catch(const std::string& e) 
	  {
	    QDPIO::cerr << InlinePropDistillutionEnv::name << ": error creating quark-line factory: " << e << endl;
	    QDP_abort(1);
	  }
	  catch(...) 
	  {
	    QDPIO::cerr << InlinePropDistillutionEnv::name << ": generic exception in creating quark-line factory" << endl;
	    QDP_abort(1);
	  }


	  // Loop over 
	  std::vector<int> t_sources(quark_line_fact->getTimeSources());

	  for(int tt=0; tt < t_sources.size(); ++tt)
	  {
	    int t_source = t_sources[tt];  // This is the pretend time-slice. The actual value is shifted.
	    QDPIO::cout << "t_source = " << t_source << endl; 

	    //
	    // Initialize all the perambulator keys
	    //
	    std::list<KeyPeramDist_t> peram_keys;
	    multi3d<ValPeramDist_t>   peram_buf;

	    if (params.named_obj.save_peramP)
	    {
	      peram_keys = quark_line_fact->getPeramKeys(t_source);
	      peram_buf.resize(Lt,Ns,Ns);

	      for(std::list<KeyPeramDist_t>::const_iterator key= peram_keys.begin();
		  key != peram_keys.end();
		  ++key)
	      {
		peram_buf(key->t_slice,key->spin_snk,key->spin_src).mat.resize(quark_line_fact->getNumVecs(),quark_line_fact->getNumSpaceDils());
	      }
	    }

	    // The space distillution loop
	    for(int dist_src=0; dist_src < quark_line_fact->getNumSpaceDils(); ++dist_src)
	    {
	      StopWatch sniss1;
	      sniss1.reset();
	      sniss1.start();
	      QDPIO::cout << "dist_src = " << dist_src << endl; 

	      // Prepare a distilluted source
	      LatticeColorVector vec_srce = quark_line_fact->getSrc(t_source, dist_src);

	      //
	      // Loop over each spin source and invert. 
	      // Use the same colorvector source. No spin dilution will be used.
	      //
	      multi2d<LatticeColorVector> ferm_out(Ns,Ns);

	      for(int spin_source=0; spin_source < Ns; ++spin_source)
	      {
		QDPIO::cout << "spin_source = " << spin_source << endl; 

		// Insert a ColorVector into spin index spin_source
		// This only overwrites sections, so need to initialize first
		LatticeFermion chi = zero;
		CvToFerm(vec_srce, chi, spin_source);

		LatticeFermion quark_soln = zero;

		// Do the propagator inversion
		SystemSolverResults_t res = (*PP)(quark_soln, chi);
		ncg_had = res.n_count;

		// Extract into the temporary output array
		for(int spin_sink=0; spin_sink < Ns; ++spin_sink)
		{
		  ferm_out(spin_sink,spin_source) = peekSpin(quark_soln, spin_sink);
		}
	      } // for spin_source


	      // Rotate from DeGrand-Rossi (DR) to Dirac-Pauli (DP)
	      {
		multi2d<LatticeColorVector> ferm_tmp;

		multiplyRep(ferm_tmp, diracToDrMatMinus, ferm_out);
		multiplyRep(ferm_out, ferm_tmp, diracToDrMatPlus);
	      }

	      sniss1.stop();
	      QDPIO::cout << "Time to assemble and transmogrify propagators for dist_src= " << dist_src << "  time = " 
			  << sniss1.getTimeInSeconds() 
			  << " secs" << endl;


	      // Write out each time-slice chunk of a lattice colorvec soln to disk
	      QDPIO::cout << "Potentially write propagator source and solutions to disk" << std::endl;
	      StopWatch sniss2;
	      sniss2.reset();
	      sniss2.start();

	      // Write the source
	      if (params.named_obj.save_srcP)
	      {
		QDPIO::cout << "Write propagator source to disk" << std::endl;
		std::list<KeyPropDist_t> src_keys(quark_line_fact->getSrcKeys(t_source, dist_src));

		for(std::list<KeyPropDist_t>::const_iterator key= src_keys.begin();
		    key != src_keys.end();
		    ++key)
		{
		  prop_obj.insert(*key, TimeSliceIO<LatticeColorVector>(vec_srce, dist_noise_obj.getTime(key->t_slice)));
		}
	      }

	      // Write the solutions
	      if (params.named_obj.save_solnP)
	      {
		QDPIO::cout << "Write propagator solution to disk" << std::endl;
		std::list<KeyPropDist_t> snk_keys(quark_line_fact->getSnkKeys(t_source, dist_src));

		for(std::list<KeyPropDist_t>::const_iterator key= snk_keys.begin();
		    key != snk_keys.end();
		    ++key)
		{
		  prop_obj.insert(*key, TimeSliceIO<LatticeColorVector>(ferm_out(key->spin_snk,key->spin_src), 
									dist_noise_obj.getTime(key->t_slice)));
		} // for key
	      }

	      // Compute perambulators
	      if (params.named_obj.save_peramP)
	      {
		QDPIO::cout << "Computing perambulators" << std::endl;

		for(std::list<KeyPeramDist_t>::const_iterator key= peram_keys.begin();
		    key != peram_keys.end();
		    ++key)
		{
		  // Get the actual time slice
		  int t_actual = dist_noise_obj.getTime(key->t_slice);

		  for(int colorvec_sink=0; colorvec_sink < quark_line_fact->getNumVecs(); ++colorvec_sink)
		  {
		    KeyTimeSliceColorVec_t key_vec;
		    key_vec.t_slice  = t_actual;
		    key_vec.colorvec = colorvec_sink;

		    LatticeColorVector tmpvec = zero;
		    TimeSliceIO<LatticeColorVector> time_slice_io(tmpvec, t_actual);

		    eigen_source.get(key_vec, time_slice_io);

		    ComplexD hsum = sum(localInnerProduct(tmpvec, ferm_out(key->spin_snk,key->spin_src)), 
					time_slice_set.getSet()[t_actual]);

		    peram_buf(key->t_slice,key->spin_snk,key->spin_src).mat(colorvec_sink,dist_src) = hsum;

		  } // for colorvec_sink
		} // for key
	      }

	      sniss2.stop();
	      QDPIO::cout << "Time to write propagators for dist_src= " << dist_src << "  time = " 
			  << sniss2.getTimeInSeconds() 
			  << " secs" << endl;

	    } // for dist_src

	  
	    // Write perambulators
	    if (params.named_obj.save_peramP)
	    {
	      QDPIO::cout << "Write perambulators to disk" << std::endl;

	      for(std::list<KeyPeramDist_t>::const_iterator key= peram_keys.begin();
		  key != peram_keys.end();
		  ++key)
	      {
		qdp_db.insert(*key, peram_buf(key->t_slice,key->spin_snk,key->spin_src));
	      } // for key
	    }

	  } // for tt
	} // for quark_line

	swatch.stop();
	QDPIO::cout << "Propagators computed: time= " 
		    << swatch.getTimeInSeconds() 
		    << " secs" << endl;
      }
      catch (const std::string& e) 
      {
	QDPIO::cout << name << ": caught exception around qprop: " << e << endl;
	QDP_abort(1);
      }

      push(xml_out,"Relaxation_Iterations");
      write(xml_out, "ncg_had", ncg_had);
      pop(xml_out);

      pop(xml_out);  // prop_dist

      snoop.stop();
      QDPIO::cout << name << ": total time = "
		  << snoop.getTimeInSeconds() 
		  << " secs" << endl;

      QDPIO::cout << name << ": ran successfully" << endl;

      END_CODE();
    } 

  }

} // namespace Chroma
