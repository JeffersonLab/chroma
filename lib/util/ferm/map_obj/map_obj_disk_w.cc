// -*- C++ -*-
/*! \file
 *  \brief Disk based map object, factory registration
 */

#include <string>
#include "chromabase.h"
#include "qdp_map_obj_disk.h"
#include "util/ferm/map_obj/map_obj_factory_w.h"
#include "util/ferm/map_obj/map_obj_disk_w.h"
#include "util/ferm/key_prop_colorvec.h"

namespace Chroma 
{ 
  
  namespace MapObjectDiskEnv 
  {

    namespace
    {
      // Parameter structure
      struct Params
      {
	Params() {}
	Params(XMLReader& xml_in, const std::string& path);

	std::string   file_name;
      };

      // Reader for input parameters
      Params::Params(XMLReader& xml, const string& path)
      {
	XMLReader paramtop(xml, path);

	read(paramtop, "FileName", file_name);
      }



      //! Callback function
      QDP::MapObject<int,EVPair<LatticeColorVector> >* createMapObjIntKeyCV(XMLReader& xml_in,
									    const std::string& path) 
      {
	// Needs parameters...
	Params params(xml_in, path);
	std::string  user_data;
	
	QDP::MapObjectDisk<int,EVPair<LatticeColorVector> >* obj = new QDP::MapObjectDisk<int,EVPair<LatticeColorVector> >();
	obj->insertUserdata(user_data);
	obj->open(params.file_name, std::ios_base::in |  std::ios_base::out |  std::ios_base::trunc);

	return obj;
      }

      //! Callback function
      QDP::MapObject<KeyPropColorVec_t,LatticeFermion>* createMapObjKeyPropColorVecLF(XMLReader& xml_in,
										      const std::string& path) 
      {
	// Needs parameters...
	Params params(xml_in, path);
	std::string  user_data;

	QDP::MapObjectDisk<KeyPropColorVec_t,LatticeFermion>* obj = new QDP::MapObjectDisk<KeyPropColorVec_t,LatticeFermion>();
	obj->insertUserdata(user_data);
	obj->open(params.file_name, std::ios_base::in |  std::ios_base::out |  std::ios_base::trunc);

	return obj;
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


} // Chroma
