// $Id: inline_nersc_read_obj.cc,v 3.2 2006-09-20 20:28:03 edwards Exp $
/*! \file
 * \brief Inline task to read an object from a named buffer
 *
 * Named object writing
 */

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/io/inline_read_map_obj_disk.h"
#include "meas/inline/io/named_objmap.h"
#include "util/ferm/map_obj.h"
#include "util/ferm/map_obj/map_obj_disk.h"
#include <string>

using namespace QDP;

namespace Chroma 
{ 
  namespace InlineReadMapObjDiskEnv 
  { 
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineReadMapObjDisk(InlineReadMapObjDiskParams(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "READ_MAP_OBJECT_DISK";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
    }
  }

  InlineReadMapObjDiskParams::InlineReadMapObjDiskParams(XMLReader& reader, 
							 const std::string& path)
  {
    try 
    {
      XMLReader paramtop(reader, path);

      if (paramtop.count("Frequency") == 1)
	read(paramtop, "Frequency", frequency);
      else
	frequency = 1;

      // Parameters for source construction
      read(paramtop, "NamedObject/object_id", named_obj.object_id);
      read(paramtop, "File/file_name", file.file_name );
    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << __func__ << ": caught Exception reading XML: " << e << endl;
      QDP_abort(1);
    }
    
  }

  void read(XMLReader& xml_in, const std::string& path, InlineReadMapObjDiskParams& p) 
  {
    InlineReadMapObjDiskParams tmp(xml_in, path);
    p = tmp;
  }

  void
  write(XMLWriter& xml_out, const std::string& path, const InlineReadMapObjDiskParams& p)
  {
    push(xml_out, path);
    
    write(xml_out, "Frequency", p.frequency);

    push(xml_out, "NamedObject");
    write(xml_out, "object_id", p.named_obj.object_id);
    pop(xml_out);

    push(xml_out, "File");
    write(xml_out, "file_name", p.file.file_name);
    pop(xml_out);

    pop(xml_out);
  }

  namespace { 
    static bool registered = false;

    int createDiskMapObjKeyPropColorVecLatticeFermion(const std::string& object_id,
						       const std::string& file_name)
    {

      Handle<MapObject<KeyPropColorVec_t, LatticeFermion> > obj_handle = 
	new MapObjectDisk<KeyPropColorVec_t, LatticeFermion>(file_name);

      obj_handle->openRead();

      TheNamedObjMap::Instance().create< 
      Handle<MapObject<KeyPropColorVec_t, LatticeFermion> >,
	Handle<MapObject<KeyPropColorVec_t, LatticeFermion> > >(object_id, obj_handle);

      return obj_handle->size();

    }

    int createDiskMapObjKeyGridPropLatticeFermion(const std::string& object_id,
						   const std::string& file_name)
    {

      Handle<MapObject<KeyGridProp_t, LatticeFermion> > obj_handle = 
	new MapObjectDisk<KeyGridProp_t, LatticeFermion>(file_name);

      obj_handle->openRead();

      TheNamedObjMap::Instance().create< 
      Handle<MapObject<KeyGridProp_t, LatticeFermion> >,
	Handle<MapObject<KeyGridProp_t, LatticeFermion> > >(object_id, obj_handle);

      return obj_handle->size();
    }

    int createDiskMapObjKeyBlockPropLatticeFermion(const std::string& object_id,
						    const std::string& file_name)
    {

      Handle<MapObject<KeyBlockProp_t, LatticeFermion> > obj_handle = 
	new MapObjectDisk<KeyBlockProp_t, LatticeFermion>(file_name);

      obj_handle->openRead();

      TheNamedObjMap::Instance().create< 
      Handle<MapObject<KeyBlockProp_t, LatticeFermion> >,
	Handle<MapObject<KeyBlockProp_t, LatticeFermion> > >(object_id, obj_handle);

      return obj_handle->size();
    }
   
    int createDiskMapObjIntLatticeColorVec(const std::string& object_id,
					    const std::string& file_name)
    {

      Handle<MapObject<int, LatticeColorVector> > obj_handle = 
	new MapObjectDisk<int, LatticeColorVector>(file_name);

      obj_handle->openRead();

      TheNamedObjMap::Instance().create< 
      Handle<MapObject<int, LatticeColorVector> >,
	Handle<MapObject<int, LatticeColorVector> > >(object_id, obj_handle);

      return obj_handle->size();
    }

    int createDiskMapObjCharFloat(const std::string& object_id,
				   const std::string& file_name)
    {

      Handle<MapObject<char, float> > obj_handle = 
	new MapObjectDisk<char, float>(file_name);

      obj_handle->openRead();

      TheNamedObjMap::Instance().create< 
      Handle<MapObject<char, float> >,
	Handle<MapObject<char, float> > >(object_id, obj_handle);

      return obj_handle->size();
    }

    std::map<int, int (*)(const std::string&, const std::string& )> funcmap;

    void registerAll(void) {
      if (! registered ) { 
	funcmap.insert( make_pair( MapObjTraitsNum<KeyPropColorVec_t, LatticeFermion>::filenum,
				   createDiskMapObjKeyPropColorVecLatticeFermion ) );

	funcmap.insert( make_pair( MapObjTraitsNum<KeyGridProp_t, LatticeFermion>::filenum,
				   createDiskMapObjKeyGridPropLatticeFermion ) );

	funcmap.insert( make_pair( MapObjTraitsNum<KeyBlockProp_t, LatticeFermion>::filenum,
				   createDiskMapObjKeyBlockPropLatticeFermion ) );

	funcmap.insert( make_pair( MapObjTraitsNum<int, LatticeColorVector>::filenum,
				   createDiskMapObjIntLatticeColorVec ) );

	funcmap.insert( make_pair( MapObjTraitsNum<char, float>::filenum,
				   createDiskMapObjCharFloat ) );

	registered = true;
      }
    }
				   
  };

  void 
  InlineReadMapObjDisk::operator()(unsigned long update_no,
				   XMLWriter& xml_out) 
  {
    START_CODE();

    registerAll();

    push(xml_out, "read_map_object_disk");
    write(xml_out, "update_no", update_no);

    QDPIO::cout << InlineReadMapObjDiskEnv::name << ": object reader" << endl;
    StopWatch swatch;

    // Read the object
    // ONLY SciDAC output format is supported in this task
    // Other tasks could support other disk formats
    QDPIO::cout << "Attempt to read object name = " << params.named_obj.object_id << endl;

    write(xml_out, "object_id", params.named_obj.object_id);
    write(xml_out, "file_name", params.file.file_name);

    try
    {
      swatch.reset();
      swatch.start();

      MapObjDiskEnv::file_typenum_t  type_num=peekMapObjectDiskTypeCode(params.file.file_name);

      int size=(funcmap[type_num])(params.named_obj.object_id, params.file.file_name);

      XMLBufferWriter file_xml_buf;
      push(file_xml_buf, "FileXML");
      write(file_xml_buf,  "object_id", params.named_obj.object_id);
      write(file_xml_buf,  "file_name", params.file.file_name);
      write(file_xml_buf,  "map_size", size);
      pop(file_xml_buf);

      XMLReader file_xml(file_xml_buf);

      TheNamedObjMap::Instance().get(params.named_obj.object_id).setFileXML( file_xml );

      // No particularly good record XML -- this after all is not QIO. So just use file_xml again
      TheNamedObjMap::Instance().get(params.named_obj.object_id).setRecordXML( file_xml );

      swatch.stop();

      QDPIO::cout << "Object successfully read: time= " 
		  << swatch.getTimeInSeconds() 
		  << " secs" << endl;
    }
    catch( std::bad_cast ) 
    {
      QDPIO::cerr << InlineReadMapObjDiskEnv::name << ": cast error" 
		  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << InlineReadMapObjDiskEnv::name << ": error message: " << e 
		  << endl;
      QDP_abort(1);
    }
    
    QDPIO::cout << InlineReadMapObjDiskEnv::name << ": ran successfully" << endl;

    pop(xml_out);  // read_named_obj

    END_CODE();
  } 

};
