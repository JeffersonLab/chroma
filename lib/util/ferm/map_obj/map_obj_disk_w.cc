// -*- C++ -*-
/*! \file
 *  \brief Disk based map object, factory registration
 */

#include <string>
#include "chromabase.h"
#include "util/ferm/map_obj.h"
#include "util/ferm/map_obj/map_obj_factory_w.h"
#include "util/ferm/map_obj/map_obj_disk.h"
#include "util/ferm/key_prop_colorvec.h"

namespace Chroma 
{ 
  
  namespace MapObjectDiskEnv 
  {

    namespace
    {
      //! Callback function
      MapObject<int,EVPair<LatticeColorVector> >* createMapObjIntKeyCV(XMLReader& xml_in,
								     const std::string& path) 
      {
	// Doesn't need parameters...
	return new MapObjectDisk<int,EVPair<LatticeColorVector> >(MapObjectDiskParams(xml_in, path));
      }

      //! Callback function
      MapObject<KeyPropColorVec_t,LatticeFermion>* createMapObjKeyPropColorVecLF(XMLReader& xml_in,
										 const std::string& path) 
      {
	// Needs parameters...
	return new MapObjectDisk<KeyPropColorVec_t,LatticeFermion>(MapObjectDiskParams(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;

      //! Name to be used
      const std::string name = "MAP_OBJECT_DISK";
    } // namespace anontmous

    std::string getName() {return name;}

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::TheMapObjIntKeyColorEigenVecFactory::Instance().registerObject(name, createMapObjIntKeyCV);
	success &= Chroma::TheMapObjKeyPropColorVecFactory::Instance().registerObject(name, createMapObjKeyPropColorVecLF);
	registered = true;
      }
      return success;
    }
  } // Namespace MapObjectDiskEnv

  void
  MapObjectDiskParams::writeXML(XMLWriter& xml_out, const std::string& path) const
  {
    push(xml_out, path);
    write(xml_out, "MapObjType", MapObjectDiskEnv::name);
    write(xml_out, "FileName", getFileName());
    pop(xml_out);
  }


} // Chroma
