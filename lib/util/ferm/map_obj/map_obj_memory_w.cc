// -*- C++ -*-
/*! \file map_obj_memory_w.cc
 *  \brief Memory based map object, factory registration
 */

#include <string>
#include "chromabase.h"
#include "util/ferm/map_obj.h"
#include "util/ferm/map_obj/map_obj_factory_w.h"
#include "util/ferm/map_obj/map_obj_memory.h"
#include "util/ferm/key_block_prop.h"
#include "util/ferm/key_grid_prop.h"
#include "util/ferm/key_prop_colorvec.h"

namespace Chroma { 
  
  namespace MapObjectMemoryEnv {
    //! Name to be used
    const std::string name("MAP_OBJ_MEMORY");

    //! Callback function
    MapObject<KeyBlockProp_t,LatticeFermion>* createMapObjKeyBlockProp(XMLReader& xml_in,
							   const std::string& path) 
    {
      // Doesn't need parameters...
      return new MapObjectMemory<KeyBlockProp_t,LatticeFermion>();
    }

    //! Callback function
    MapObject<KeyGridProp_t,LatticeFermion>* createMapObjKeyGridProp(XMLReader& xml_in,
							   const std::string& path) 
    {
      // Doesn't need parameters...
      return new MapObjectMemory<KeyGridProp_t,LatticeFermion>();
    }


    //! Callback function
    MapObject<KeyPropColorVec_t,LatticeFermion>* createMapObjKeyPropColorVec(XMLReader& xml_in,
							   const std::string& path) 
    {
      // Doesn't need parameters...
      return new MapObjectMemory<KeyPropColorVec_t,LatticeFermion>();
    }

    
    //! Registration flag
    static bool registered = false;

    //! Registration function
    bool registerAll()
    {
      bool success = true;
      if( !registered ) {
	// Factory: MapObject<KeyBlockProp_t,LatticeFermion>
	success &= Chroma::TheMapObjKeyBlockPropFactory::Instance().registerObject(name, createMapObjKeyBlockProp );

	// Factory: MapObject<KeyGridProp_t,LatticeFermion>
	success &= Chroma::TheMapObjKeyGridPropFactory::Instance().registerObject(name, createMapObjKeyGridProp );

	// Factory: MapObject<KeyPropColorVec_t,LatticeFermion>
	success &= Chroma::TheMapObjKeyPropColorVecFactory::Instance().registerObject(name, createMapObjKeyPropColorVec );
	
	registered = true;
      }
      return success;
    }
  } // Namespace MapObjectMemoryEnv
} // Chroma
