/*! \file
 * \brief Distillution factory for producing keys * sources
 *
 * Distillution factory for producing keys * sources
 */

#include "meas/hadron/distillution_factory.h"

namespace Chroma 
{ 
  //----------------------------------------------------------------------------
  //! Get sink key
  KeyPropDistillution_t AbsQuarkLine::getSnkKey(const KeyPeramDistillution_t& peram_key, int dist_src) const
  {
    KeyPropDistillution_t snk_key;

    snk_key.prop_type  = "SNK";
    snk_key.annihP     = peram_key.annihP;
    snk_key.t_source   = peram_key.t_source;
    snk_key.t_slice    = peram_key.t_slice;
    snk_key.dist_src   = dist_src;
    snk_key.spin_src   = peram_key.spin_src;
    snk_key.spin_snk   = peram_key.spin_snk;
    snk_key.quark_line = peram_key.quark_line;
    snk_key.mass       = peram_key.mass;

    return snk_key;
  }



  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  //! Utilities
  namespace DistillutionFactoryEnv
  {
    // More utilities
    namespace
    {
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
	//! Constructor
	QuarkLineFact(const Params& params_,
		      const DistillutionNoise& dist_noise_obj_, 
		      QDP::MapObjectDisk<KeyPropDistillution_t, TimeSliceIO<LatticeColorVectorF> >& source_obj_,
		      const TimeSliceSet& time_slice_set_,
		      int quark_line_,
		      const std::string& mass_);

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

	//! Get mass
	virtual std::string getMass() const {return mass;}

	//! Get annihilation flag
	virtual bool getAnnihP() const {return false;}

	//! Get the time sources
	virtual std::vector<int> getTimeSources() const {return params.t_sources;}

	//! Get source keys
	virtual std::list<KeyPropDistillution_t> getSrcKeys(int t_source, int dist_src) const;

	//! Get sink keys
	virtual std::list<KeyPropDistillution_t> getSnkKeys(int t_source, int dist_src) const;

	//! Get perambulator keys
	virtual std::list<KeyPeramDistillution_t> getPeramKeys(int t_source) const;

	//! Get perambulator key time slices
	virtual std::list<int> getTslices(int t_source) const;

      private:
	//! The active time slices for this source
	virtual std::vector<bool> getActiveTSlices(int t_source) const;

	//! Get source key
	virtual KeyPropDistillution_t getSrcKey(int t_source, int dist_src) const;

      private:
	// Arguments
	Params                    params;
	const DistillutionNoise&  dist_noise_obj;
	QDP::MapObjectDisk<KeyPropDistillution_t, TimeSliceIO<LatticeColorVectorF> >& source_obj;
	const TimeSliceSet&       time_slice_set;
	int                       quark_line;
	std::string               mass;
      };



      //----------------------------------------------------------------------------
      //! Constructor
      QuarkLineFact::QuarkLineFact(const Params& params_,
				   const DistillutionNoise& dist_noise_obj_, 
				   QDP::MapObjectDisk<KeyPropDistillution_t, TimeSliceIO<LatticeColorVectorF> >& source_obj_,
				   const TimeSliceSet& time_slice_set_,
				   int quark_line_,
				   const std::string& mass_)
	: params(params_), 
	  dist_noise_obj(dist_noise_obj_), source_obj(source_obj_), time_slice_set(time_slice_set_),
	  quark_line(quark_line_), mass(mass_)
      {
      }

      //----------------------------------------------------------------------------
      //! Prepare a distilluted source
      LatticeColorVector QuarkLineFact::getSrc(int t_source, int dist_src) const
      {
	QDPIO::cout << "CONN: getSrc on t_source= " << t_source << "  dist_src= " << dist_src << endl;

	// Get the actual time slice
	int t_actual = dist_noise_obj.getTime(t_source);

	// Get the source vector
	KeyPropDistillution_t src_key = getSrcKey(t_source, dist_src);
	LatticeColorVectorF vec_srce = zero;

	TimeSliceIO<LatticeColorVectorF> time_slice_io(vec_srce, t_actual);
	source_obj.get(src_key, time_slice_io);

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
      //! Get source key
      KeyPropDistillution_t QuarkLineFact::getSrcKey(int t_source, int dist_src) const
      {
	KeyPropDistillution_t key;

	key.prop_type    = "SRC";
	key.annihP       = false;
	key.t_source     = t_source;
	key.t_slice      = t_source;
	key.dist_src     = dist_src;
	key.spin_src     = -1;
	key.spin_snk     = -1;
	key.quark_line   = quark_line;
	key.mass         = mass;

	return key;
      }

	
      //----------------------------------------------------------------------------
      //! Get source keys
      std::list<KeyPropDistillution_t> QuarkLineFact::getSrcKeys(int t_source, int dist_src) const
      {
	std::list<KeyPropDistillution_t> keys;

	keys.push_back(getSrcKey(t_source,dist_src));

	return keys;
      }

	
      //----------------------------------------------------------------------------
      //! Get sink keys
      std::list<KeyPropDistillution_t> QuarkLineFact::getSnkKeys(int t_source, int dist_src) const
      {
	std::list<KeyPropDistillution_t> keys;

	std::vector<bool> active_t_slices = getActiveTSlices(t_source);
	
	const int Lt = Layout::lattSize()[dist_noise_obj.getDecayDir()];

	for(int spin_source=0; spin_source < Ns; ++spin_source)
	{
	  for(int spin_sink=0; spin_sink < Ns; ++spin_sink)
	  {
	    for(int t=0; t < Lt; ++t)
	    {
	      if (! active_t_slices[t]) {continue;}

	      KeyPropDistillution_t key;

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
      std::list<KeyPeramDistillution_t> QuarkLineFact::getPeramKeys(int t_source) const
      {
	std::list<KeyPeramDistillution_t> keys;

	std::vector<bool> active_t_slices = getActiveTSlices(t_source);
	
	const int Lt = Layout::lattSize()[dist_noise_obj.getDecayDir()];

	for(int spin_source=0; spin_source < Ns; ++spin_source)
	{
	  for(int spin_sink=0; spin_sink < Ns; ++spin_sink)
	  {
	    for(int t=0; t < Lt; ++t)
	    {
	      if (! active_t_slices[t]) {continue;}

	      KeyPeramDistillution_t key;

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

	
      //----------------------------------------------------------------------------
      //! Get perambulator key time slices
      std::list<int> QuarkLineFact::getTslices(int t_source) const
      {
	std::list<int> keys;

	std::vector<bool> active_t_slices = getActiveTSlices(t_source);
	
	const int Lt = Layout::lattSize()[dist_noise_obj.getDecayDir()];

	for(int t=0; t < Lt; ++t)
	{
	  if (active_t_slices[t])
	  {
	    keys.push_back(t);
	  }
	} // for t

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
	//! Constructor
	QuarkLineFact(const Params& params_,
		      const DistillutionNoise& dist_noise_obj_, 
		      QDP::MapObjectDisk<KeyPropDistillution_t, TimeSliceIO<LatticeColorVectorF> >& source_obj_,
		      const TimeSliceSet& time_slice_set_,
		      int quark_line_,
		      const std::string& mass_);

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

	//! Get mass
	virtual std::string getMass() const {return mass;}

	//! Get annihilation flag
	virtual bool getAnnihP() const {return true;}

	//! Get the time sources
	virtual std::vector<int> getTimeSources() const;

	//! Get source keys
	virtual std::list<KeyPropDistillution_t> getSrcKeys(int t_source, int dist_src) const;

	//! Get sink keys
	virtual std::list<KeyPropDistillution_t> getSnkKeys(int t_source, int dist_src) const;

	//! Get perambulator keys
	virtual std::list<KeyPeramDistillution_t> getPeramKeys(int t_source) const;

	//! Get perambulator key time slices
	virtual std::list<int> getTslices(int t_source) const;

      private:
	//! Get source key
	virtual KeyPropDistillution_t getSrcKey(int t_source, int dist_src) const;

      private:
	// Arguments
	Params                    params;
	const DistillutionNoise&  dist_noise_obj;
	QDP::MapObjectDisk<KeyPropDistillution_t, TimeSliceIO<LatticeColorVectorF> >& source_obj;
	const TimeSliceSet&       time_slice_set;
	int                       quark_line;
	std::string               mass;
      };



      //----------------------------------------------------------------------------
      //! Constructor
      QuarkLineFact::QuarkLineFact(const Params& params_,
				   const DistillutionNoise& dist_noise_obj_, 
				   QDP::MapObjectDisk<KeyPropDistillution_t, TimeSliceIO<LatticeColorVectorF> >& source_obj_,
				   const TimeSliceSet& time_slice_set_,
				   int quark_line_,
				   const std::string& mass_)
	: params(params_), 
	  dist_noise_obj(dist_noise_obj_), source_obj(source_obj_), time_slice_set(time_slice_set_),
	  quark_line(quark_line_), mass(mass_)
      {
	// Reset/barf if bogus
	params.num_time_dils  = checkTimeDils(params.num_time_dils, Layout::lattSize()[dist_noise_obj.getDecayDir()]);
      }


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
	QDPIO::cout << "ANNIH: getSrc on t_source= " << t_src << "  dist_src= " << dist_src << endl;

	LatticeColorVector vec_srce = zero;

	const int Lt = Layout::lattSize()[dist_noise_obj.getDecayDir()];

	for(int time_source=t_src; time_source < Lt; time_source += params.num_time_dils)
	{
	  // Get the actual time slice
	  int t_actual = dist_noise_obj.getTime(time_source);

	  // Get the source vector
	  KeyPropDistillution_t src_key = getSrcKey(time_source, dist_src);
	  LatticeColorVectorF vec_tmp = zero;
	  
	  TimeSliceIO<LatticeColorVectorF> time_slice_io(vec_tmp, t_actual);
	  source_obj.get(src_key, time_slice_io);

	  vec_srce[time_slice_set.getSet()[t_actual]] += vec_tmp;
	} // for time_source

	return vec_srce;
      }


      //----------------------------------------------------------------------------
      //! Get source keys
      KeyPropDistillution_t QuarkLineFact::getSrcKey(int t_source, int dist_src) const
      {
	KeyPropDistillution_t key;

	key.prop_type    = "SRC";
	key.annihP       = true;
	key.t_source     = t_source;
	key.t_slice      = t_source;
	key.dist_src     = dist_src;
	key.spin_src     = -1;
	key.spin_snk     = -1;
	key.quark_line   = quark_line;
	key.mass         = mass;

	return key;
      }

	
      //----------------------------------------------------------------------------
      //! Get source keys
      std::list<KeyPropDistillution_t> QuarkLineFact::getSrcKeys(int t_source, int dist_src) const
      {
	std::list<KeyPropDistillution_t> keys;

	const int Lt = Layout::lattSize()[dist_noise_obj.getDecayDir()];

	for(int t=t_source; t < Lt; t += params.num_time_dils)
	{
	  keys.push_back(getSrcKey(t, dist_src));
	}

	return keys;
      }

	
      //----------------------------------------------------------------------------
      //! Get sink keys
      std::list<KeyPropDistillution_t> QuarkLineFact::getSnkKeys(int t_source, int dist_src) const
      {
	std::list<KeyPropDistillution_t> keys;

	const int Lt = Layout::lattSize()[dist_noise_obj.getDecayDir()];

	for(int spin_source=0; spin_source < Ns; ++spin_source)
	{
	  for(int spin_sink=0; spin_sink < Ns; ++spin_sink)
	  {
	    for(int t=t_source; t < Lt; t += params.num_time_dils)
	    {
	      KeyPropDistillution_t key;

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
      std::list<KeyPeramDistillution_t> QuarkLineFact::getPeramKeys(int t_source) const
      {
	std::list<KeyPeramDistillution_t> keys;

	const int Lt = Layout::lattSize()[dist_noise_obj.getDecayDir()];

	for(int spin_source=0; spin_source < Ns; ++spin_source)
	{
	  for(int spin_sink=0; spin_sink < Ns; ++spin_sink)
	  {
	    for(int t=t_source; t < Lt; t += params.num_time_dils)
	    {
	      KeyPeramDistillution_t key;

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

      //----------------------------------------------------------------------------
      //! Get perambulator key time slices
      std::list<int> QuarkLineFact::getTslices(int t_source) const
      {
	std::list<int> keys;

	const int Lt = Layout::lattSize()[dist_noise_obj.getDecayDir()];

	for(int t=t_source; t < Lt; t += params.num_time_dils)
	{
	  keys.push_back(t);
	} // for t

	return keys;
      }

    } // namespace Annihilation


    //----------------------------------------------------------------------------
    namespace
    {
      AbsQuarkLine* createConn(XMLReader& xml_in, 
			       const std::string& path,
			       const DistillutionNoise& dist_noise_obj, 
			       QDP::MapObjectDisk<KeyPropDistillution_t, TimeSliceIO<LatticeColorVectorF> >& source_obj,
			       const TimeSliceSet& time_slice_set,
			       int quark_line,
			       const std::string& mass)
      {
	return new Connected::QuarkLineFact(Connected::Params(xml_in, path),
					    dist_noise_obj, source_obj, time_slice_set,
					    quark_line, mass);
      }

      AbsQuarkLine* createAnnih(XMLReader& xml_in, 
				const std::string& path,
				const DistillutionNoise& dist_noise_obj, 
				QDP::MapObjectDisk<KeyPropDistillution_t, TimeSliceIO<LatticeColorVectorF> >& source_obj,
				const TimeSliceSet& time_slice_set,
				int quark_line,
				const std::string& mass)
      {
	return new Annihilation::QuarkLineFact(Annihilation::Params(xml_in, path),
					       dist_noise_obj, source_obj, time_slice_set,
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
    
  } // namespace DistillutionFactoryEnv


} // namespace Chroma
