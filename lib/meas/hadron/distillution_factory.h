// -*- C++ -*-
/*! \file
 * \brief Distillution factory for producing keys * sources
 *
 * Distillution factory for producing keys * sources
 */

#ifndef __distillution_factory_w_h__
#define __distillution_factory_w_h__

#include "qdp_map_obj_disk.h"
#include "qdp_disk_map_slice.h"

#include "chromabase.h"
#include "singleton.h"
#include "objfactory.h"
#include "meas/inline/abs_inline_measurement.h"
#include "util/ferm/key_prop_distillution.h"
#include "util/ferm/key_peram_distillution.h"
#include "util/ferm/distillution_noise.h"
#include "util/ft/time_slice_set.h"
#include <list>
#include <vector>

namespace Chroma 
{ 
  //----------------------------------------------------------------------------
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

    //! Get mass
    virtual std::string getMass() const = 0;

    //! Get annihilation flag
    virtual bool getAnnihP() const = 0;

    //! Get the time sources
    virtual std::vector<int> getTimeSources() const = 0;

    //! Get source keys
    virtual std::list<KeyPropDistillution_t> getSrcKeys(int t_source, int dist_src) const = 0;

    //! Get sink keys
    virtual std::list<KeyPropDistillution_t> getSnkKeys(int t_source, int dist_src) const = 0;

    //! Get perambulator keys
    virtual std::list<KeyPeramDistillution_t> getPeramKeys(int t_source) const = 0;

    //! Get sink key
    virtual KeyPropDistillution_t getSnkKey(const KeyPeramDistillution_t& peram_key, int dist_src) const;

    //! Get perambulator key time slices
    virtual std::list<int> getTslices(int t_source) const = 0;
  };


  //----------------------------------------------------------------------------
  typedef QDP::MapObjectDisk< KeyPropDistillution_t,TimeSliceIO<LatticeColorVectorF> > MOD_t;

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
		  AbsQuarkLine* (*)(XMLReader&, const std::string&,
				    const DistillutionNoise&, 
				    MOD_t&,
				    const TimeSliceSet&,
				    int,
				    const std::string&), StringFactoryError> >
  TheQuarkLineFactory;


  //----------------------------------------------------------------------------
  namespace DistillutionFactoryEnv 
  {
    bool registerAll();
  }

}

#endif
