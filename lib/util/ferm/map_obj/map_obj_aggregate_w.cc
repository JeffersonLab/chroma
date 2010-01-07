// -*- C++ -*-
/*! \file map_obj_aggregate_w.cc 
 *  \brief Aggregate registration for MapObj types
 */

#include "util/ferm/map_obj/map_obj_aggregate_w.h"

// Individual MapObj headers
#include "util/ferm/map_obj/map_obj_memory_w.h"
#include "util/ferm/map_obj/map_obj_disk_w.h"
#include "util/ferm/map_obj/map_obj_null_w.h"


namespace Chroma {
  
  //! Registration accregator 
  namespace MapObjectWilson4DEnv 
  {

    // Anonymous namespace
    namespace
    {
      bool registered = false;
    }

    bool registerAll()
    {
      bool success = true;
      if (! registered) 
      { 
	success &= MapObjectDiskEnv::registerAll();
	success &= MapObjectMemoryEnv::registerAll();
	success &= MapObjectNullEnv::registerAll();

	registered = true;
      }
      return success;
    }
      


  } // MapObjectEnv

}
