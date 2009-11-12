// -*- C++ -*-
/*! \file map_obj_aggregate_w.cc 
 *  \brief Aggregate registration for MapObj types
 */

#include "util/ferm/map_obj/map_obj_aggregate_w.h"

// Individual MapObj headers
#include "util/ferm/map_obj/map_obj_memory_w.h"


namespace Chroma {
  
  //! Registration accregator 
  namespace MapObjectWilson4DEnv {

    static bool registeredKeyBlockPropLF = false;
    bool registerKeyBlockPropLFAll() 
    {
      bool success = true;
      if (!registeredKeyBlockPropLF) {
	success &= Chroma::MapObjectMemoryEnv::registerKeyBlockPropLF();

	registeredKeyBlockPropLF = true;
      }
      return success;
    }

    static bool registeredKeyGridPropLF = false;


    bool registerKeyGridPropLFAll() 
    {
      bool success = true;
      if (!registeredKeyGridPropLF) {
	success &= Chroma::MapObjectMemoryEnv::registerKeyGridPropLF();

	registeredKeyGridPropLF = true;
      }
      return success;
    }

    static bool registeredKeyPropColorVecLF = false;
    bool registerKeyPropColorVecLFAll() 
    {
      bool success = true;
      if (!registeredKeyPropColorVecLF ) {
	success &= Chroma::MapObjectMemoryEnv::registerKeyPropColorVecLF();
	
	registeredKeyPropColorVecLF = true;
      }


      return success;
    }


    static bool registered = false;
    bool registerAll()
    {
      bool success = true;
      if( !registered) { 
	success &= registerKeyBlockPropLFAll();
	success &= registerKeyGridPropLFAll();
	success &= registerKeyPropColorVecLFAll();
	
	registered = true;
      }
      return success;
    }
      


  } // MapObjectEnv

}
