/*! \file
 * \brief Inline task to read an object into a named buffer
 */

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/io/inline_read_map_obj_memory.h"
#include "meas/inline/io/named_objmap.h"
#include "util/ferm/subset_ev_pair.h"
#include "util/ferm/map_obj.h"
#include "util/ferm/map_obj/map_obj_disk.h"
#include "util/ferm/map_obj/map_obj_memory.h"
#include <string>

namespace Chroma 
{ 
  namespace InlineReadMapObjMemoryEnv 
  { 
    namespace ReadMapObjCallEnv
    { 
      namespace 
      { 
	static bool registered = false;

	template<typename K, typename V>
	int createMapObj(const std::string& object_id,
			 const std::string& file_name)
	{
	  MapObjectDisk<K,V> diskMap(file_name);
	  diskMap.openRead();
	  std::vector<K> keys = diskMap.dump();

	  Handle<MapObject<K,V> > obj_handle(new MapObjectMemory<K,V>());
	  TheNamedObjMap::Instance().create< Handle<MapObject<K,V> >, Handle<MapObject<K,V> > >(object_id, obj_handle);

	  obj_handle->openWrite();

	  for(int i=0; i < keys.size(); i++) 
	  {
	    V v;
	    diskMap.lookup(keys[i], v);
	    obj_handle->insert(keys[i],v);
	  }
	  obj_handle->openRead();

	  return obj_handle->size();
	}


	std::map<std::string, int (*)(const std::string&, const std::string& )> funcmap;

	void registerAll(void) 
	{
	  if (! registered ) 
	  { 
	    funcmap.insert( make_pair( MapObjTraitsNum<KeyPropColorVec_t, LatticeFermion>::type_string,
				       createMapObj<KeyPropColorVec_t, LatticeFermion> ) );

	    funcmap.insert( make_pair( MapObjTraitsNum<int, EVPair<LatticeColorVector> >::type_string,
				       createMapObj<int, EVPair<LatticeColorVector> > ) );

	    funcmap.insert( make_pair( MapObjTraitsNum<char, float>::type_string,
				       createMapObj<char, float> ) );

	    registered = true;
	  }
	}
      }
    }


    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineMeas(Params(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;

      const std::string name = "READ_MAP_OBJECT_DISK";
    }

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


    Params::Params(XMLReader& reader, const std::string& path)
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


    void 
    InlineMeas::operator()(unsigned long update_no, XMLWriter& xml_out) 
    {
      START_CODE();

      registerAll();

      push(xml_out, "read_map_object_disk");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << InlineReadMapObjMemoryEnv::name << ": object reader" << endl;
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

	std::string  type_string=peekMapObjectDiskTypeCode(params.file.file_name);

	int size=(ReadMapObjCallEnv::funcmap[type_string])(params.named_obj.object_id, params.file.file_name);

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
	QDPIO::cerr << InlineReadMapObjMemoryEnv::name << ": cast error" 
		    << endl;
	QDP_abort(1);
      }
      catch (const string& e) 
      {
	QDPIO::cerr << InlineReadMapObjMemoryEnv::name << ": error message: " << e 
		    << endl;
	QDP_abort(1);
      }
    
      QDPIO::cout << InlineReadMapObjMemoryEnv::name << ": ran successfully" << endl;

      pop(xml_out);  // read_named_obj

      END_CODE();
    } 

  }
}
