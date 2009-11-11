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
    MapObject<KeyBlockProp_t,LatticeFermion>* createMapObjKeyBlockPropLF(XMLReader& xml_in,
							   const std::string& path) 
    {
      // Doesn't need parameters...
      return new MapObjectMemory<KeyBlockProp_t,LatticeFermion>();
    }

    //! Callback function
    MapObject<KeyGridProp_t,LatticeFermion>* createMapObjKeyGridPropLF(XMLReader& xml_in,
							   const std::string& path) 
    {
      // Doesn't need parameters...
      return new MapObjectMemory<KeyGridProp_t,LatticeFermion>();
    }


    //! Callback function
    MapObject<KeyPropColorVec_t,LatticeFermion>* createMapObjKeyPropColorVecLF(XMLReader& xml_in,
							   const std::string& path) 
    {
      // Doesn't need parameters...
      return new MapObjectMemory<KeyPropColorVec_t,LatticeFermion>();
    }

    
    //! Registration flag
    static bool RegisteredKeyBlockPropLF = false;

    //! Registration function
    bool registerKeyBlockPropLF()
    {
      bool success = true;
      if( !RegisteredKeyBlockPropLF ) {
	// Factory: MapObject<KeyBlockProp_t,LatticeFermion>
	success &= Chroma::TheMapObjKeyBlockPropFactory::Instance().registerObject(name, createMapObjKeyBlockPropLF );

	RegisteredKeyBlockPropLF = true;
      }
      return success;
    }


    static bool RegisteredKeyGridPropLF = false;
    bool registerKeyGridPropLF()
    {
      bool success = true;
      if( !RegisteredKeyGridPropLF ) {
	success &= Chroma::TheMapObjKeyGridPropFactory::Instance().registerObject(name, createMapObjKeyGridPropLF );

	RegisteredKeyGridPropLF = true;
      }
      return success;
    }

    static bool RegisteredKeyPropColorVecLF = false;
    bool registerKeyPropColorVecLF()
    {
      bool success = true;
      if( !RegisteredKeyPropColorVecLF ) {
	// Factory: MapObject<KeyPropColorVec_t,LatticeFermion>
	success &= Chroma::TheMapObjKeyPropColorVecFactory::Instance().registerObject(name, createMapObjKeyPropColorVecLF );
	
	RegisteredKeyPropColorVecLF = true;
      }
      return success;
    }


  } // Namespace MapObjectMemoryEnv
} // Chroma
