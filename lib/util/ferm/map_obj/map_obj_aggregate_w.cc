// -*- C++ -*-
/*! \file map_obj_aggregate_w.cc 
 *  \brief Aggregate registration for MapObj types
 */

#include "util/ferm/map_obj/map_obj_aggregate_w.h"

// Individual MapObj headers
#include "util/ferm/map_obj/map_obj_memory.h"


namespace Chroma {
  
  //! Registration accregator 
  namespace MapObjectWilson4DEnv {
    static bool registered = false;

    bool registerAll() 
    {
      bool success = true;
      if (!registered) {
	success &= MapObjectMemoryEnv::registerAll();
      }
      return success;
    }



  } // MapObjectEnv

}
